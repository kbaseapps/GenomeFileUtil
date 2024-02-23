import copy
import datetime
import hashlib
import os
import re
import shutil
import sys
import time
import uuid
from collections import Counter, defaultdict, OrderedDict
import logging

import Bio.SeqIO
import Bio.SeqUtils
from Bio import Seq
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from installed_clients.WorkspaceClient import Workspace
from GenomeFileUtil.core.GenomeUtils import (
    is_parent, propagate_cds_props_to_gene, warnings, parse_inferences,
    load_ontology_mappings, set_taxon_data, set_default_taxon_data
)

MAX_MISC_FEATURE_SIZE = 10000
MAX_PARENT_LOOKUPS = 5

# catalog params
MAX_THREADS_DEFAULT = 10
THREADS_PER_CPU_DEFAULT = 1

_WSID = 'workspace_id'
_WSNAME = 'workspace_name'
_INPUTS = 'inputs'


def _upa(object_info):
    return f'{object_info[6]}/{object_info[0]}/{object_info[4]}'


class _Genome:
    def __init__(self):
        self.time_string = str(datetime.datetime.fromtimestamp(
            time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
        self.generate_parents = False
        self.generate_ids = False
        self.genes = OrderedDict()
        self.mrnas = OrderedDict()
        self.cdss = OrderedDict()
        self.noncoding = []
        self.ontologies_present = defaultdict(dict)
        self.ontology_events = list()
        self.skiped_features = Counter()
        self.feature_counts = Counter()
        self.orphan_types = Counter()
        self.contig_seq = {}
        self.circ_contigs = set()
        self.features_spaning_zero = set()
        self.genome_warnings = []
        self.genome_suspect = False
        self.defects = Counter()
        self.spoofed_genes = 0
        self.excluded_features = ('source', 'exon', 'fasta_record')
        self.ont_mappings = load_ontology_mappings('/kb/module/data')
        self.code_table = 11
        # dict with feature 'id's that have been used more than once.
        self.used_twice_identifiers = {}

        # related info for genome process and upload
        self.genome_name = None
        self.genome_data = None
        self.genome_meta = None
        self.genome_type = None
        self.consolidated_file = None
        self.input_directory = None
        self.handle_service_output = None
        self.extra_info = None
        self.assembly_ref = None
        self.assembly_path = None
        self.assembly_info = None


class GenbankToGenome:
    def __init__(self, config):
        self.cfg = config
        self.gi = GenomeInterface(config)
        self.dfu = DataFileUtil(config.callbackURL)
        self.aUtil = AssemblyUtil(config.callbackURL)
        self.ws = Workspace(config.workspaceURL)
        self.re_api_url = config.re_api_url
        self.version = config.version
        self.default_params = {
            'source': 'Genbank',
            'taxon_wsname': config.raw['taxon-workspace-name'],
            'taxon_lookup_obj_name': config.raw['taxon-lookup-object-name'],

            'ontology_wsname': config.raw['ontology-workspace-name'],
            'ontology_GO_obj_name': config.raw['ontology-gene-ontology-obj-name'],
            'ontology_PO_obj_name': config.raw['ontology-plant-ontology-obj-name'],

            'release': None,
            'genetic_code': 11,
            'generate_ids_if_needed': 0,
            'metadata': {}
        }

    def import_genbank(self, params):
        print('validating parameters')
        mass_params = self._set_up_single_params(params)
        return self._import_genbank_mass(mass_params)[0]

    def import_genbank_mass(self, params):
        print('validating parameters')
        self._validate_mass_params(params)
        return self._import_genbank_mass(params)

    def _set_up_single_params(self, params):
        # avoid side effects and keep variables in params unmodfied
        inputs = dict(params)
        self._validate_params(inputs)
        ws_id = self._get_int(inputs.pop(_WSID, None), _WSID)
        ws_name = inputs.pop(_WSNAME, None)
        if (bool(ws_id) == bool(ws_name)):  # xnor
            raise ValueError(f"Exactly one of a '{_WSID}' or a '{_WSNAME}' parameter must be provided")
        if not ws_id:
            print(f"Translating workspace name {ws_name} to a workspace ID. Prefer submitting "
                  + "a workspace ID over a mutable workspace name that may cause race conditions")
            ws_id = self.dfu.ws_name_to_id(ws_name)
        mass_params = {_WSID: ws_id, _INPUTS: [inputs]}
        return mass_params

    def _validate_mass_params(self, params):
        ws_id = self._get_int(params.get(_WSID), _WSID)
        if not ws_id:
            raise ValueError(f"{_WSID} is required")
        inputs = params.get(_INPUTS)
        if not inputs or type(inputs) is not list:
            raise ValueError(f"{_INPUTS} field is required and must be a non-empty list")
        for i, inp in enumerate(inputs, start=1):
            if type(inp) is not dict:
                raise ValueError(f"Entry #{i} in {_INPUTS} field is not a mapping as required")
            try:
                self._validate_params(inp)
            except Exception as e:
                raise ValueError(f"Entry #{i} in {_INPUTS} field has invalid params") from e

    def _get_int(self, putative_int, name, minimum=1):
        if putative_int is not None:
            if type(putative_int) is not int:
                raise ValueError(f"{name} must be an integer, got: {putative_int}")
            if putative_int < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        return putative_int

    def _import_genbank_mass(self, params):

        workspace_id = params[_WSID]
        inputs = params[_INPUTS]
        genome_objs = []

        for idx, input_params in enumerate(inputs):

            genome_obj = _Genome()
            genome_obj.genome_name = input_params['genome_name']

            # construct the input directory staging area
            genome_obj.input_directory = self._stage_input(input_params)

            # update default params
            input_params = {**self.default_params, **input_params}
            genome_obj.genome_meta = input_params['metadata']
            inputs[idx] = input_params

            genome_obj.genome_type = input_params.get('genome_type', 'isolate')
            genome_obj.generate_parents = input_params.get('generate_missing_genes')
            genome_obj.generate_ids = input_params.get('generate_ids_if_needed')
            genome_obj.use_existing_assembly = input_params.get('use_existing_assembly')
            if input_params.get('genetic_code'):
                genome_obj.code_table = input_params['genetic_code']

            # find genbank file
            files = self._find_input_files(genome_obj.input_directory)
            genome_obj.consolidated_file = self._join_files_skip_empty_lines(files)

            # gather all objects
            genome_objs.append(genome_obj)

        self._get_contigs_and_validate_existing_assembly(genome_objs)

        self._save_files_to_blobstore_and_set_handle_service_output(genome_objs)

        self._save_assemblies(workspace_id, genome_objs)

        # TODO parse genbank files before the assemblies are saved
        for input_params, genome_obj in zip(inputs, genome_objs):
            self._parse_genbank(input_params, genome_obj)

            # clear the temp directory
            # TODO move with parse genbank function into for loop above and save disk space
            shutil.rmtree(genome_obj.input_directory)

        # TODO make an internal mass function save_genomes
        results = self._save_genomes(workspace_id, genome_objs)

        # return the result
        details = [
            {
                'genome_ref': _upa(result["info"]),
                'genome_info': result["info"],
                'assembly_ref': genome_obj.assembly_ref,
                'assembly_path': genome_obj.assembly_path,
                'assembly_info': genome_obj.assembly_info,
            }
            for genome_obj, result in zip(genome_objs, results)
        ]

        return details

    def _save_genomes(self, workspace_id, genome_objs):
        results = [
            self.gi.save_one_genome(
                {
                    'workspace': workspace_id,
                    'name': genome_obj.genome_name,
                    'data': genome_obj.genome_data,
                    "meta": genome_obj.genome_meta,
                }
            ) for genome_obj in genome_objs
        ]

        return results

    def _validate_params(self, params):
        if 'genome_name' not in params:
            raise ValueError('required "genome_name" field was not defined')
        if 'file' not in params:
            raise ValueError('required "file" field was not defined')

        # one and only one of 'path', 'shock_id', or 'ftp_url' is required
        file = params['file']
        if not isinstance(file, dict):
            raise ValueError('required "file" field must be a map/dict')
        sources = ('path', 'shock_id', 'ftp_url')
        n_valid_fields = sum(1 for f in sources if file.get(f))
        if n_valid_fields < 1:
            raise ValueError(f'required "file" field must include one source: '
                             f'{", ".join(sources)}')
        if n_valid_fields > 1:
            raise ValueError(f'required "file" field has too many sources specified: '
                             f'{", ".join(file.keys())}')
        if params.get('genetic_code'):
            if not (isinstance(params['genetic_code'], int) and 0 < params['genetic_code'] < 32):
                raise ValueError(f"Invalid genetic code specified: {params}")

    def _stage_input(self, params):
        """ Setup the input_directory by fetching the files and uncompressing if needed. """

        # construct the input directory where we stage files
        input_directory = os.path.join(self.cfg.sharedFolder,
                                       f'genome-upload-staging-{uuid.uuid4()}')
        os.makedirs(input_directory)

        # at this point, the 'file' input is validated, so we don't have to catch any special cases
        # we expect one and only one of path, shock_id, or ftp_url

        # determine how to get the file: if it is from shock, download it.  If it
        # is just sitting there, then use it.  Move the file to the staging input directory
        file = params['file']
        genbank_file_path = None
        if file.get('path') is not None:
            # copy the local file to the input staging directory
            # (NOTE: could just move it, but then this method would have the side effect of moving your
            # file which another SDK module might have an open handle on)
            local_file_path = file['path']
            genbank_file_path = os.path.join(input_directory, os.path.basename(local_file_path))
            shutil.copy2(local_file_path, genbank_file_path)

        if 'shock_id' in file and file['shock_id'] is not None:
            # handle shock file
            sys.stdout.flush()
            file_name = self.dfu.shock_to_file({
                                    'file_path': input_directory,
                                    'shock_id': file['shock_id']
                                })['node_file_name']
            genbank_file_path = os.path.join(input_directory, file_name)

        if 'ftp_url' in file and file['ftp_url'] is not None:
            local_file_path = self.dfu.download_web_file({
                'file_url': file['ftp_url'],
                'download_type': 'FTP'
            })['copy_file_path']
            genbank_file_path = os.path.join(input_directory,
                                             os.path.basename(local_file_path))
            shutil.copy2(local_file_path, genbank_file_path)

        # extract the file if it is compressed
        if genbank_file_path is not None:
            self.dfu.unpack_file({'file_path': genbank_file_path})

        else:
            raise ValueError('No valid files could be extracted based on the input')

        return input_directory

    def _save_files_to_blobstore_and_set_handle_service_output(self, genome_objs):
        logging.info("Saving original files to shock")
        handle_service_outputs = self.dfu.file_to_shock_mass(
            [
                {
                    'file_path': genome_obj.consolidated_file,
                    'make_handle': 1,
                    'pack': 'gzip',
                } for genome_obj in genome_objs
            ]
        )
        for genome_obj, handle_service_output in zip(
            genome_objs, handle_service_outputs
        ):
            genome_obj.handle_service_output = handle_service_output

    def _get_objects_data(self, assembly_refs):
        assembly_objs_spec = [{"ref": ref} for ref in assembly_refs]
        assembly_objs_info = Workspace.get_object_info3(
            {"objects": assembly_objs_spec, "includeMetadata": 1})["infos"]
        return assembly_objs_info

    def _parse_genbank(self, params, genome_obj):
        genome = {
            "id": params['genome_name'],
            "original_source_file_name": os.path.basename(genome_obj.consolidated_file),
            "assembly_ref": genome_obj.assembly_ref,
            "gc_content": genome_obj.gc_content,
            "dna_size": genome_obj.dna_size,
            "md5": genome_obj.md5,
            "genbank_handle_ref": genome_obj.handle_service_output['handle']['hid'],
            "publications": set(),
            "contig_ids": [],
            "contig_lengths": [],
        }
        genome['source'], genome['genome_tiers'] = self.gi.determine_tier(params['source'])

        if params.get('genome_type'):
            genome['genome_type'] = params['genome_type']

        # Set taxonomy-related fields in the genome
        # Also validates the given taxon ID
        if params.get('taxon_id'):
            set_taxon_data(int(params['taxon_id']), self.re_api_url, genome)
        else:
            set_default_taxon_data(genome)

        dates = []
        # Parse data from genbank file
        contigs = Bio.SeqIO.parse(genome_obj.consolidated_file, "genbank")
        for record in contigs:
            r_annot = record.annotations
            try:
                dates.append(time.strptime(r_annot.get('date'), "%d-%b-%Y"))
            except (TypeError, ValueError):
                pass
            genome['contig_ids'].append(record.id)
            genome['contig_lengths'].append(len(record))
            genome["publications"] |= self._get_pubs(r_annot)

            # only do the following once(on the first contig)
            if "source_id" not in genome:
                genome["source_id"] = record.id.split('.')[0]
                organism = r_annot.get('organism', 'Unknown Organism')
                if params.get('scientific_name'):
                    genome['scientific_name'] = params['scientific_name']
                else:
                    genome['scientific_name'] = organism
                genome_obj.code_table = genome['genetic_code']
                genome["molecule_type"] = r_annot.get('molecule_type', 'DNA')
                genome['notes'] = r_annot.get('comment', "").replace('\\n', '\n')

            self._parse_features(record, genome['source'], genome_obj)

        genome.update(self._get_feature_lists(genome_obj))

        genome['num_contigs'] = len(genome['contig_ids'])
        # add dates
        dates.sort()
        if dates:
            genome['external_source_origination_date'] = time.strftime(
                "%d-%b-%Y", dates[0])
            if dates[0] != dates[-1]:
                genome['external_source_origination_date'] += " _ " + \
                    time.strftime("%d-%b-%Y", dates[-1])

        if genome_obj.ontologies_present:
            genome['ontologies_present'] = dict(genome_obj.ontologies_present)
            genome["ontology_events"] = genome_obj.ontology_events
        genome['feature_counts'] = dict(genome_obj.feature_counts)
        # can't serialize a set
        genome['publications'] = list(genome['publications'])

        if len(genome['cdss']) and (genome_obj.defects['cds_seq_not_matching'] /
                                    float(len(genome['cdss'])) > 0.02):
            genome_obj.genome_warnings.append(
                warnings["genome_inc_translation"].format(
                    genome_obj.defects['cds_seq_not_matching'], len(genome['cdss'])
                )
            )
            genome_obj.genome_suspect = 1

        if genome_obj.defects['bad_parent_loc']:
            genome_obj.genome_warnings.append(
                f"There were {genome_obj.defects['bad_parent_loc']} parent/child "
                "relationships that were not able to be determined. Some of "
                "these may have splice variants that may be valid relationships."
            )

        if genome_obj.defects['spoofed_genes']:
            genome_obj.genome_warnings.append(
                warnings['spoofed_genome'].format(genome_obj.defects['spoofed_genes'])
            )
            genome['suspect'] = 1

        if genome_obj.defects['not_trans_spliced']:
            genome_obj.genome_warnings.append(
                warnings['genome_not_trans_spliced'].format(
                    genome_obj.defects['not_trans_spliced']
                )
            )
            genome['suspect'] = 1

        if genome_obj.genome_warnings:
            genome['warnings'] = genome_obj.genome_warnings
        if genome_obj.genome_suspect:
            genome['suspect'] = 1

        if params.get('genetic_code'):
            genome['genetic_code'] = params['genetic_code']

        genome_obj.genome_data = genome

    def _validate_existing_assembly(self, assembly_ref, genome_obj):
        if not re.match("\d+\/\d+\/\d+", assembly_ref):
            raise ValueError(f"Assembly ref: {assembly_ref} is not a valid format. Must"
                                f" be in numerical <ws>/<object>/<version> format.")
        ret = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
        if "KBaseGenomeAnnotations.Assembly" not in ret['info'][2]:
            raise ValueError(f"{assembly_ref} is not a reference to an assembly")
        unmatched_ids = list()
        unmatched_ids_md5s = list()
        for current_contig in genome_obj.contig_seq.keys():
            current_contig_md5 = hashlib.md5(
                str(genome_obj.contig_seq[current_contig]).encode('utf8')
            ).hexdigest()
            if current_contig in ret['data']['contigs']:
                if current_contig_md5 != ret['data']['contigs'][current_contig]['md5']:
                    unmatched_ids_md5s.append(current_contig)
            else:
                unmatched_ids.append(current_contig)
        if len(unmatched_ids) > 0:
            raise ValueError(warnings['assembly_ref_extra_contigs'].format(", ".join(unmatched_ids)))
        if len(unmatched_ids_md5s) > 0:
            raise ValueError(warnings["assembly_ref_diff_seq"].format(", ".join(unmatched_ids_md5s)))

    def _get_contigs_and_extra_info(self, contigs, genome_obj):
        out_contigs = []
        extra_info = defaultdict(dict)
        for in_contig in contigs:
            if in_contig.annotations.get('topology', "") == 'circular':
                extra_info[in_contig.id]['is_circ'] = 1
                genome_obj.circ_contigs.add(in_contig.id)
            elif in_contig.annotations.get('topology', "") == 'linear':
                extra_info[in_contig.id]['is_circ'] = 0
            out_contigs.append(in_contig)
            genome_obj.contig_seq[in_contig.id] = in_contig.seq.upper()
        return out_contigs, extra_info

    def _get_contigs_and_validate_existing_assembly(self, genome_objs):
        ref2genome = {}
        for genome_obj in genome_objs:
            assembly_ref = genome_obj.use_existing_assembly
            if assembly_ref:
                contigs = Bio.SeqIO.parse(genome_obj.consolidated_file, "genbank")
                _, extra_info = self._get_contigs_and_extra_info(
                    contigs, genome_obj
                )
                self._validate_existing_assembly(
                    assembly_ref, genome_obj
                )
                ref2genome[assembly_ref] = genome_obj
                genome_obj.extra_info = extra_info
                genome_obj.assembly_ref = assembly_ref
                genome_obj.assembly_path = None

        if ref2genome:
            assembly_refs = list(ref2genome.keys())
            assembly_objs_info = self._get_objects_data(assembly_refs)
            for assembly_ref, assembly_info in zip(
                assembly_refs, assembly_objs_info
            ):
                genome = ref2genome[assembly_ref]
                self._set_gc_content_dna_size_and_md5(genome, assembly_info)

    def _set_gc_content_dna_size_and_md5(self, genome, assembly_info):
        genome.assembly_info = assembly_info
        genome.gc_content = float(assembly_info[10]["GC content"])
        genome.dna_size = int(assembly_info[10]["Size"])
        genome.md5 = assembly_info[10]["MD5"]

    def _save_assemblies(self, workspace_id, genome_objs):
        id2genome = {}
        bulk_inputs = []
        for genome_obj in genome_objs:
            if genome_obj.assembly_ref:
                continue
            contigs = Bio.SeqIO.parse(genome_obj.consolidated_file, "genbank")
            genome_obj.assembly_id = f"{genome_obj.genome_name}_assembly"
            genome_obj.assembly_path = f"{self.cfg.sharedFolder}/{genome_obj.assembly_id}.fasta"

            out_contigs, genome_obj.extra_info = self._get_contigs_and_extra_info(
                contigs, genome_obj
            )

            Bio.SeqIO.write(out_contigs, genome_obj.assembly_path, "fasta")
            bulk_inputs.append(
                {
                    'file': genome_obj.assembly_path,
                    'assembly_name': genome_obj.assembly_id,
                    'type': genome_obj.genome_type,
                    'contig_info': genome_obj.extra_info,
                }
            )
            # map id to genome_object
            id2genome[genome_obj.assembly_id] = genome_obj

        if id2genome:
            assembly_refs = self.aUtil.save_assemblies_from_fastas(
                {
                    'workspace_id': workspace_id,
                    'inputs': bulk_inputs,
                }
            )["results"]

            for result in assembly_refs:
                assembly_info = result["object_info"]
                genome = id2genome[assembly_info[1]]
                genome.assembly_ref = result["upa"]
                self._set_gc_content_dna_size_and_md5(genome, assembly_info)

            logging.info(f"{len(id2genome)} assemblies saved to {workspace_id}")

    def _find_input_files(self, input_directory):
        valid_extensions = [".gbff", ".gbk", ".gb", ".genbank", ".dat", ".gbf"]

        files = os.listdir(os.path.abspath(input_directory))
        genbank_files = [x for x in files if
                         os.path.splitext(x)[-1].lower() in valid_extensions]

        if len(genbank_files) == 0:
            raise Exception(
                f"The input directory does not have any files with one of the "
                f"following extensions {','.join(valid_extensions)}.")

        input_files = []
        for genbank_file in genbank_files:
            input_files.append(os.path.join(input_directory, genbank_file))

        return input_files

    def _join_files_skip_empty_lines(self, input_files):
        """ Applies strip to each line of each input file.
        Args:
            input_files: Paths to input files in Genbank format.
        Returns:
            Path to resulting file (currenly it's the same file as input).
        """
        if len(input_files) == 0:
            raise ValueError("NO GENBANK FILE")
        temp_dir = os.path.join(os.path.dirname(input_files[0]), "combined")
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        ret_file = os.path.join(temp_dir, os.path.basename(input_files[0]))

        # take in Genbank file and remove all empty lines from it.
        with open(ret_file, 'w', buffering=2 ** 20) as f_out:
            for input_file in input_files:
                with open(input_file, 'r') as f_in:
                    for line in f_in:
                        line = line.rstrip('\r\n')
                        if line.strip():
                            f_out.write(line + '\n')
        return ret_file

    def _get_pubs(self, r_annotations):
        """Get a contig's publications"""
        pub_list = []
        for in_pub in r_annotations.get('references', []):
            # don't add blank pubs
            if not in_pub.authors:
                continue
            out_pub = [
                0,  # pmid
                "",  # source
                in_pub.title,
                "",  # web address
                "",  # date
                in_pub.authors,
                in_pub.journal,
            ]
            date_match = re.match("\((\d{4})\)", in_pub.journal)
            if date_match:
                out_pub[4] = date_match.group(1)
            if in_pub.pubmed_id:
                out_pub[0:4] = [
                    int(in_pub.pubmed_id),
                    "PubMed",
                    in_pub.title,
                    f"http://www.ncbi.nlm.nih.gov/pubmed/{in_pub.pubmed_id}"]
            pub_list.append(tuple(out_pub))
        return set(pub_list)

    def _get_id(self, feat, genome_obj, tags=None):
        """Assign a id to a feature based on the first tag that exists"""
        _id = ""
        if not tags:
            tags = ['locus_tag', 'kbase_id']
        for t in tags:
            _id = feat.qualifiers.get(t, [""])[0]
            if _id:
                break

        if not _id:
            if feat.type == 'gene':
                if not genome_obj.generate_ids:
                    raise ValueError(f"Unable to find a valid id for gene "
                                     f"among these tags: {', '.join(tags)}. Correct the "
                                     f"file or rerun with generate_ids\n {feat}")
                genome_obj.orphan_types['gene'] += 1
                _id = f"gene_{genome_obj.orphan_types['gene']}"
            if 'rna' in feat.type.lower() or feat.type in {'CDS', 'sig_peptide',
                                                           'five_prime_UTR', 'three_prime_UTR'}:
                _id = f"gene_{genome_obj.orphan_types['gene']}"

        return _id

    def _parse_features(self, record, source, genome_obj):
        def _location(feat):
            """Convert to KBase style location objects"""
            strand_trans = ("", "+", "-")
            loc = []
            for part in feat.location.parts:
                contig_id = part.ref if part.ref else record.id
                if part.strand >= 0:
                    begin = int(part.start) + 1
                else:
                    begin = int(part.end)
                loc.append((
                        contig_id,
                        begin,
                        strand_trans[part.strand],
                        len(part)))
            return loc

        def _warn(message):
            if message not in out_feat.get('warnings', []):
                out_feat['warnings'] = out_feat.get('warnings', []) + [message]

        def _check_suspect_location(parent=None):
            if 'trans_splicing' in out_feat.get('flags', []):
                return

            if out_feat['location'] == sorted(out_feat['location'],
                    reverse=(in_feature.location.strand == -1)):
                return

            if record.id in genome_obj.circ_contigs and \
                    in_feature.location.start == 0 \
                    and in_feature.location.end == len(record):
                genome_obj.features_spaning_zero.add(out_feat['id'])
                return

            if parent and parent['id'] in genome_obj.features_spaning_zero:
                return

            _warn(warnings['not_trans_spliced'])
            genome_obj.defects['not_trans_spliced'] += 1

        for in_feature in record.features:
            if in_feature.type in genome_obj.excluded_features:
                genome_obj.skiped_features[in_feature.type] += 1
                continue
            feat_seq = self._get_seq(in_feature, record.id, genome_obj)
            if source == "Ensembl":
                _id = self._get_id(
                    in_feature, genome_obj, ['gene', 'locus_tag']
                )
            else:
                _id = self._get_id(in_feature, genome_obj)

            # The following is common to all the feature types
            out_feat = {
                "id": "_".join([_id, in_feature.type]),
                "location": _location(in_feature),
                "dna_sequence": str(feat_seq),
                "dna_sequence_length": len(feat_seq),
                "md5": hashlib.md5(str(feat_seq).encode('utf8')).hexdigest(),
            }
            if not _id:
                out_feat['id'] = in_feature.type

            # validate input feature
            # note that end is the larger number regardless of strand
            if int(in_feature.location.end) > len(record):
                genome_obj.genome_warnings.append(
                    warnings["coordinates_off_end"].format(out_feat['id']))
                genome_obj.genome_suspect = 1
                continue

            for piece in in_feature.location.parts:
                if not isinstance(piece.start, ExactPosition) \
                        or not isinstance(piece.end, ExactPosition):
                    _warn(warnings["non_exact_coordinates"])

            genome_obj.feature_counts[in_feature.type] += 1

            # add optional fields
            if 'note' in in_feature.qualifiers:
                out_feat['note'] = in_feature.qualifiers["note"][0]

            out_feat.update(self._get_aliases_flags_functions(in_feature))

            ont, db_xrefs = self._get_ontology_db_xrefs(in_feature, genome_obj)
            if ont:
                out_feat['ontology_terms'] = ont
            if db_xrefs:
                out_feat['db_xrefs'] = db_xrefs

            if 'inference' in in_feature.qualifiers:
                out_feat['inference_data'] = parse_inferences(
                    in_feature.qualifiers['inference'])

            _check_suspect_location(genome_obj.genes.get(_id))

            # add type specific features
            if in_feature.type == 'CDS':
                self._process_cds(
                    _id, feat_seq, in_feature, out_feat, genome_obj
                )

            elif in_feature.type == 'gene':
                self._process_gene(_id, out_feat, genome_obj)

            elif in_feature.type == 'mRNA':
                self._process_mrna(_id, out_feat, genome_obj)

            else:
                genome_obj.noncoding.append(
                    self._process_noncoding(
                        _id, in_feature.type, out_feat, genome_obj
                    )
                )

    def _get_feature_lists(self, genome_obj):
        """sort genes into their final arrays"""
        coding = []
        for g in genome_obj.genes.values():
            if len(g['cdss']):
                if g['mrnas'] and len(g['mrnas']) != len(g['cdss']):
                    msg = "The length of the mrna and cdss arrays are not equal"
                    g['warnings'] = g.get('warnings', []) + [msg]

                # remove duplicates that may arise from CDS info propagation
                for key in ('functions', 'aliases', 'db_xrefs'):
                    if key in g:
                        g[key] = list(set(g[key]))
                if not g['mrnas']:
                    del g['mrnas']
                del g['type']
                coding.append(g)
                genome_obj.feature_counts["protein_encoding_gene"] += 1
            else:
                del g['mrnas'], g['cdss']
                genome_obj.noncoding.append(g)
                genome_obj.feature_counts["non_coding_genes"] += 1

        genome_obj.feature_counts["non_coding_features"] = len(genome_obj.noncoding)
        return {
            'features': coding,
            'non_coding_features': genome_obj.noncoding,
            'cdss': list(genome_obj.cdss.values()),
            'mrnas': list(genome_obj.mrnas.values()),
        }

    def _get_seq(self, feat, contig, genome_obj):
        """Extract the DNA sequence for a feature"""
        seq = []
        for part in feat.location.parts:
            strand = part.strand
            # handle trans-splicing across contigs
            if part.ref:
                part_contig = part.ref
            else:
                part_contig = contig

            if strand >= 0:
                seq.append(str(genome_obj.contig_seq[part_contig]
                               [part.start:part.end]))
            else:
                seq.append(str(genome_obj.contig_seq[part_contig]
                               [part.start:part.end].reverse_complement()))
        return "".join(seq)

    def _create_ontology_event(self, ontology_type, genome_obj):
        """Creates the ontology_event if necessary
        Returns the index of the ontology event back."""
        if ontology_type not in genome_obj.ont_mappings:
            raise ValueError(f"{ontology_type} is not a supported ontology")

        if "event_index" not in genome_obj.ont_mappings[ontology_type]:
            genome_obj.ont_mappings[ontology_type]['event_index'] = len(genome_obj.ontology_events)
            if ontology_type == "GO":
                ontology_ref = "KBaseOntology/gene_ontology"
            elif ontology_type == "PO":
                ontology_ref = "KBaseOntology/plant_ontology"
            else:
                ontology_ref = f"KBaseOntology/{ontology_type.lower()}_ontology"
            genome_obj.ontology_events.append({
                "method": "GenomeFileUtils Genbank uploader from annotations",
                "method_version": self.version,
                "timestamp": genome_obj.time_string,
                "id": ontology_type,
                "ontology_ref": ontology_ref
            })

        return genome_obj.ont_mappings[ontology_type]['event_index']

    def _get_ontology_db_xrefs(self, feature, genome_obj):
        """Splits the ontology info from the other db_xrefs"""
        ontology = defaultdict(dict)
        db_xrefs = []
        for key in ("GO_process", "GO_function", "GO_component"):
            ontology_event_index = self._create_ontology_event("GO", genome_obj)
            for term in feature.qualifiers.get(key, []):
                sp = term.split(" - ")
                ontology['GO'][sp[0]] = [ontology_event_index]
                genome_obj.ontologies_present['GO'][sp[0]] = genome_obj.ont_mappings['GO'].get(sp[0], '')

        for ref in feature.qualifiers.get('db_xref', []):
            if ref.startswith('GO:'):
                ontology['GO'][ref] = [self._create_ontology_event("GO", genome_obj)]
                genome_obj.ontologies_present['GO'][ref] = genome_obj.ont_mappings['GO'].get(ref, '')
            elif ref.startswith('PO:'):
                ontology['PO'][ref] = [self._create_ontology_event("PO", genome_obj)]
                genome_obj.ontologies_present['PO'][ref] = genome_obj.ont_mappings['PO'].get(ref, '')
            elif ref.startswith('KO:'):
                ontology['KO'][ref] = [self._create_ontology_event("KO", genome_obj)]
                genome_obj.ontologies_present['KO'][ref] = genome_obj.ont_mappings['KO'].get(ref, '')
            elif ref.startswith('COG'):
                ontology['COG'][ref] = [self._create_ontology_event("COG", genome_obj)]
                genome_obj.ontologies_present['COG'][ref] = genome_obj.ont_mappings['COG'].get(ref, '')
            elif ref.startswith('PF'):
                ontology['PFAM'][ref] = [self._create_ontology_event("PFAM", genome_obj)]
                genome_obj.ontologies_present['PFAM'][ref] = genome_obj.ont_mappings['PFAM'].get(ref, '')
            elif ref.startswith('TIGR'):
                ontology['TIGRFAM'][ref] = [self._create_ontology_event("TIGRFAM", genome_obj)]
                genome_obj.ontologies_present['TIGRFAM'][ref] = genome_obj.ont_mappings['TIGRFAM'].get(ref, '')
            elif ":" not in ref:
                db_xrefs.append(tuple(["Unknown_Source", ref]))
            else:
                db_xrefs.append(tuple(ref.split(":", 1)))

        return dict(ontology), sorted(db_xrefs)

    def _get_aliases_flags_functions(self, feat):
        """Get the values for aliases flags and features from qualifiers"""
        alias_keys = {'locus_tag', 'old_locus_tag', 'protein_id',
                      'transcript_id', 'gene', 'EC_number', 'gene_synonym'}
        result = defaultdict(list)
        for key, val_list in feat.qualifiers.items():
            if key in alias_keys:
                result['aliases'].extend([(key, val) for val in val_list])
            # flags have no other information associated with them
            if val_list == ['']:
                result['flags'].append(key)
            if key == 'function':
                result['functional_descriptions'].extend(val_list[0].split('; '))
            if key == 'product':
                result['functions'] = val_list

        return result

    def _find_parent_gene(self, potential_id, feature, genome_obj):
        """Unfortunately, Genbank files don't have a parent ID and the features can be out of
        order at times. To account for this, the this function works backwards from the end of
        list of IDs and stops when if finds a parent with valid coordinates or it hits the maximum
        number of tries"""
        if potential_id in genome_obj.genes:
            lookup_attempts = 0
            while lookup_attempts < MAX_PARENT_LOOKUPS:
                if is_parent(genome_obj.genes[potential_id], feature):
                    return potential_id

                lookup_attempts += 1
                try:
                    potential_id = list(genome_obj.genes.keys())[-(lookup_attempts + 1)]
                except IndexError:
                    break  # no more genes that could match exist

            genome_obj.defects['bad_parent_loc'] += 1
        return None

    def _assign_new_id(self, _id, genome_obj):
        """given a feature id that has already been used, add a unique modifier to it"""
        _id_modifier = genome_obj.used_twice_identifiers.get(_id, 1)
        genome_obj.used_twice_identifiers[_id] = _id_modifier + 1
        return _id + "." + str(_id_modifier)

    def _process_gene(self, _id, out_feat, genome_obj):
        out_feat.update({
            "id": _id,
            "type": 'gene',
            "mrnas": [],
            'cdss': [],
        })
        if _id in genome_obj.genes:
            _id = self._assign_new_id(_id, genome_obj)
            out_feat.update({"id": _id})
            # raise ValueError(f"Duplicate gene ID: {_id}")
        genome_obj.genes[_id] = out_feat

    def _process_noncoding(self, gene_id, feat_type, out_feat, genome_obj):
        out_feat["type"] = feat_type

        # this prevents big misc_features from blowing up the genome size
        if out_feat['dna_sequence_length'] > MAX_MISC_FEATURE_SIZE:
            del out_feat['dna_sequence']

        gene_id = self._find_parent_gene(gene_id, out_feat, genome_obj)
        if gene_id:
            if 'children' not in genome_obj.genes[gene_id]:
                genome_obj.genes[gene_id]['children'] = []
            out_feat['id'] += "_" + str(len(genome_obj.genes[gene_id]['children']) + 1)
            genome_obj.genes[gene_id]['children'].append(out_feat['id'])
            out_feat['parent_gene'] = gene_id
        else:
            genome_obj.orphan_types[feat_type] += 1
            out_feat['id'] += "_" + str(genome_obj.orphan_types[feat_type])

        return out_feat

    def _process_mrna(self, gene_id, out_feat, genome_obj):
        if gene_id not in genome_obj.genes and genome_obj.generate_parents:
            self._process_gene(gene_id, copy.copy(out_feat), genome_obj)

        gene_id = self._find_parent_gene(gene_id, out_feat, genome_obj)
        if gene_id:
            out_feat['id'] = "_".join(
                (gene_id, "mRNA", str(len(genome_obj.genes[gene_id]['mrnas']) + 1))
            )
            genome_obj.genes[gene_id]['mrnas'].append(out_feat['id'])
            out_feat['parent_gene'] = gene_id
        else:
            genome_obj.orphan_types['mrna'] += 1
            out_feat['id'] = f"mRNA_{genome_obj.orphan_types['mrna']}"
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                'Unable to find parent gene for ' + str(out_feat['id'])]

        genome_obj.mrnas[out_feat['id']] = out_feat

    def _process_cds(self, gene_id, feat_seq, in_feature, out_feat, genome_obj):
        # Associate CDS with parents
        cds_warnings = out_feat.get('warnings', [])
        validated_gene_id = self._find_parent_gene(gene_id, out_feat, genome_obj)
        if validated_gene_id:
            out_feat['id'] = "_".join((validated_gene_id, "CDS",
                                       str(len(genome_obj.genes[validated_gene_id]['cdss']) + 1)))
            genome_obj.genes[validated_gene_id]['cdss'].append(out_feat['id'])
            out_feat['parent_gene'] = validated_gene_id
        elif genome_obj.generate_parents and gene_id not in genome_obj.genes:
            new_feat = copy.copy(out_feat)
            new_feat['id'] = gene_id
            new_feat['warnings'] = [warnings['spoofed_gene']]
            genome_obj.orphan_types['gene'] += 1
            genome_obj.defects['spoofed_genes'] += 1
            self._process_gene(new_feat['id'], new_feat, genome_obj)

            out_feat['id'] = "_".join((gene_id, "CDS", str(len(genome_obj.genes[gene_id]['cdss']) + 1)))
            genome_obj.genes[gene_id]['cdss'].append(out_feat['id'])
            out_feat['parent_gene'] = gene_id
        else:
            genome_obj.orphan_types['cds'] += 1
            out_feat['id'] = f"CDS_{genome_obj.orphan_types['cds']}"
            cds_warnings.append(f"Unable to find parent gene for {out_feat['id']}")

        # there is a 1 to 1 relationship of mRNA to CDS so XXX_mRNA_1 will match XXX_CDS_1
        mrna_id = out_feat["id"].replace('CDS', 'mRNA')
        if mrna_id in genome_obj.mrnas:
            if not is_parent(genome_obj.mrnas[mrna_id], out_feat):
                cds_warnings.append(warnings['cds_mrna_cds'].format(mrna_id))
                genome_obj.mrnas[mrna_id]['warnings'] = genome_obj.mrnas[mrna_id].get(
                    'warnings', []) + [warnings['cds_mrna_mrna']]
                genome_obj.defects['bad_parent_loc'] += 1
            else:
                out_feat['parent_mrna'] = mrna_id
                genome_obj.mrnas[mrna_id]['cds'] = out_feat['id']

        # process protein
        prot_seq = in_feature.qualifiers.get("translation", [""])[0]

        # allow a little slack to account for frameshift and stop codon
        if prot_seq and abs(len(prot_seq) * 3 - len(feat_seq)) > 4:
            cds_warnings.append(warnings["inconsistent_CDS_length"].format(len(feat_seq),
                                                                           len(prot_seq)))
            genome_obj.genome_warnings.append(
                warnings['genome_inc_CDS_length'].format(
                    out_feat['id'], len(feat_seq), len(prot_seq)))
            genome_obj.genome_suspect = 1

        try:
            if prot_seq and prot_seq != Seq.translate(
                    feat_seq, genome_obj.code_table, cds=True).strip("*"):
                cds_warnings.append(warnings["inconsistent_translation"])
                genome_obj.defects['cds_seq_not_matching'] += 1

        except TranslationError as e:
            cds_warnings.append("Unable to verify protein sequence:" + str(e))

        if not prot_seq:
            try:
                prot_seq = Seq.translate(
                        feat_seq, genome_obj.code_table, cds=True).strip("*")
                cds_warnings.append(warnings["no_translation_supplied"])

            except TranslationError as e:
                cds_warnings.append(warnings["no_translation_supplied"] + str(e))

        out_feat.update({
            "protein_translation": prot_seq,
            "protein_md5": hashlib.md5(prot_seq.encode('utf8')).hexdigest(),
            "protein_translation_length": len(prot_seq),
        })

        if out_feat.get('parent_gene'):
            propagate_cds_props_to_gene(
                out_feat, genome_obj.genes[out_feat['parent_gene']]
            )

        if cds_warnings:
            out_feat['warnings'] = cds_warnings

        genome_obj.cdss[out_feat['id']] = out_feat
