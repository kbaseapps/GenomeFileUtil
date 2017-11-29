import os
import sys
import shutil
import uuid
import time
import json
import gzip
import copy
import hashlib
import collections
import datetime
import re

# KBase imports
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeInterface import GenomeInterface

# 3rd party imports
from Bio.Data import CodonTable
import Bio.SeqIO

codon_table = CodonTable.ambiguous_generic_by_name["Standard"]

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class FastaGFFToGenome:

    def __init__(self, config):
        self.cfg = config
        self.au = AssemblyUtil(config.callbackURL)
        self.dfu = DataFileUtil(self.cfg.callbackURL)
        self.gi = GenomeInterface(self.cfg)
        self.taxon_wsname = self.cfg.raw['taxon-workspace-name']
        self.time_string = str(datetime.datetime.fromtimestamp(
            time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
        yml_text = open('/kb/module/kbase.yml').read()
        self.version = re.search("module-version:\n\W+(.+)\n", yml_text).group(
            1)
        self.is_phytozome = False
        self._warnings = []
        self.feature_dict = {}
        self.ontologies_present = collections.defaultdict(dict)
        self.skiped_features = collections.Counter()
        self.feature_counts = collections.Counter()

    def warn(self, message):
        self._warnings.append(message)
        print message

    def import_file(self, params):

        # 1) validate parameters
        self._validate_import_file_params(params)

        # 2) construct the input directory staging area
        input_directory = os.path.join(self.cfg.sharedFolder, 'fast_gff_upload_'+str(uuid.uuid4()))
        os.makedirs(input_directory)
        file_paths = self._stage_input(params, input_directory)

        # 3) extract out the parameters
        params = self._set_parsed_params(params)

        # 4) do the upload
        result = self.upload_genome(
            input_fasta_file=file_paths["fasta_file"],
            input_gff_file=file_paths["gff_file"],

            workspace_name=params['workspace_name'],
            core_genome_name=params['genome_name'],
            scientific_name=params['scientific_name'],
            source=params['source'],
            genome_type=params['type'],
            release=params['release'],
            metadata=params['metadata']
        )

        # 5) generate report
        output_data_ref = params['workspace_name']+"/"+params['genome_name']
        reportObj = {'objects_created': [{'ref': output_data_ref,
                                          'description': 'KBase Genome object'}],
                     'text_message': result['report_string'],
                     'warnings': self._warnings}

        reportClient = KBaseReport(os.environ['SDK_CALLBACK_URL'])
        report_info = reportClient.create({'report': reportObj,
                                           'workspace_name': params['workspace_name']})

        # 6) clear the temp directory
        shutil.rmtree(input_directory)

        # 7) return the result
        info = result['genome_info']
        details = {
            'genome_ref': str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
            'genome_info': info,
            'report_name': report_info['name'],
            'report_ref':  report_info['ref']
        }

        return details

    def upload_genome(self, input_gff_file=None, input_fasta_file=None,
                      workspace_name=None, core_genome_name=None,
                      scientific_name="unknown_taxon", source=None,
                      release=None, genome_type=None, metadata=None):

        # save assembly file
        assembly_ref = self.au.save_assembly_from_fasta(
            {'file': {'path': input_fasta_file},
             'workspace_name': workspace_name,
             'assembly_name': core_genome_name + ".assembly"})
        assembly_data = self.dfu.get_objects(
            {'object_refs': [assembly_ref],
             'ignore_errors': 0})['data'][0]['data']

        # reading in GFF file
        features_by_contig = self._retrieve_gff_file(input_gff_file)
        contig_ids = set(assembly_data['contigs'])
        for cid in set(features_by_contig.keys()) - contig_ids:
            self.warn("Sequence name {} does not match a sequence id in the "
                      "FASTA file. {} features will not be imported."
                      .format(cid, len(features_by_contig[cid])))

        # parse feature information
        fasta_contigs = Bio.SeqIO.parse(input_fasta_file, "fasta")
        for contig in fasta_contigs:
            molecule_type = str(contig.seq.alphabet).replace(
                'IUPACAmbiguous', '').strip('()')
            for feature in features_by_contig.get(contig.id, []):
                if feature['type'] == 'exon':
                    self._parse_exon(feature)
                else:
                    self._transform_feature(contig, feature)

        # generate genome info
        genome = self._gen_genome_info(core_genome_name, scientific_name,
                                       assembly_ref, source, assembly_data,
                                       input_gff_file, molecule_type)

        result = self.gi.save_one_genome({
            'workspace': workspace_name,
            'name': core_genome_name,
            'data': genome,
            "meta": metadata,
        })
        report_string = ''

        return {'genome_info': result['info'], 'report_string': report_string}

    def _parse_exon(self, feature):
        parent_id = feature.get('Parent', '')
        if not parent_id:
            self.warn("An exon has no valid parent: {}".format(feature))
            return
        parent = self.feature_dict[parent_id]
        # If the parent does not have this exon flag, the parent location
        # should be overwritten.
        if parent.get('exons', False):
            self.feature_dict[parent_id]['location'].append(
                self._location(feature))
        else:
            parent['exons'] = True
            self.feature_dict[parent_id]['location'] = [self._location(feature)]

    @staticmethod
    def _location(in_feature):
        return [
            in_feature['contig'],
            in_feature['start'],
            in_feature['strand'],
            in_feature['end'] - in_feature['start']
        ]

    def _get_ontology(self, feature):
        ontology = collections.defaultdict(dict)
        for key in ("GO_process", "GO_function", "GO_component"):
            if key in feature['attributes']:
                sp = feature['attributes'][key][0][3:].split(" - ")
                ontology['GO'][sp[0]] = [1]
                self.ontologies_present['GO'][sp[0]] = sp[1]
        # TODO: Support other ontologies
        return dict(ontology)

    @staticmethod
    def _validate_import_file_params(params):
        """
        validate_import_file_params:
                    validates params passed to FastaGFFToGenome.import_file method

        """

        # check for required parameters
        for p in ['workspace_name', 'genome_name', 'fasta_file', 'gff_file']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        # one and only one of 'path', or 'shock_id' is required
        for key in ('fasta_file', 'gff_file'):
            file = params[key]
            if not isinstance(file, dict):
                raise ValueError('Required "{}" field must be a map/dict'.format(key))
            n_valid_fields = 0
            if 'path' in file and file['path'] is not None:
                n_valid_fields += 1
            if 'shock_id' in file and file['shock_id'] is not None:
                n_valid_fields += 1
            if 'ftp_url' in file and file['ftp_url'] is not None:
                n_valid_fields += 1
                raise ValueError('FTP link is currently not supported for FastaGFFToGenome')
            if n_valid_fields < 1:
                error_msg = 'Required "{}" field must include one source: '.format(key)
                error_msg += 'path | shock_id'
                raise ValueError(error_msg)
            if n_valid_fields > 1:
                error_msg = 'Required "{}" field has too many sources specified: '.format(key)
                error_msg += str(file.keys())
                raise ValueError(error_msg)

        # check for valid type param
        valid_types = ['Reference', 'User upload', 'Representative']
        if params.get('type') and params['type'] not in valid_types:
            error_msg = 'Entered value for type is not one of the valid entries of '
            error_msg += '[' + ''.join('"' + str(e) + '", ' for e in valid_types)[0: -2] + ']'
            raise ValueError(error_msg)

    def _set_parsed_params(self, params):
        log('Setting params')

        # default params
        default_params = {
            'taxon_wsname': self.cfg.raw['taxon-workspace-name'],
            'scientific_name': 'unknown_taxon',
            'taxon_reference': None,
            'source': 'User',
            'release': None,
            'type': 'User upload',
            'metadata': {}
        }

        for field in default_params:
            if field not in params:
                params[field] = default_params[field]

        log(json.dumps(params, indent=1))

        return params

    def _stage_input(self, params, input_directory):
        """
        stage_input: Setup the input_directory by fetching the files and uncompressing if needed

        """

        file_paths = dict()
        for key in ('fasta_file', 'gff_file'):
            file = params[key]
            file_path = None
            if 'path' in file and file['path'] is not None:
                local_file_path = file['path']
                file_path = os.path.join(input_directory, os.path.basename(local_file_path))
                log('Moving file from {} to {}'.format(local_file_path, file_path))
                shutil.copy2(local_file_path, file_path)

            if 'shock_id' in file and file['shock_id'] is not None:
                # handle shock file
                log('Downloading file from SHOCK node: {}-{}'.format(
                                                        self.cfg.sharedFolder, file['shock_id']))
                sys.stdout.flush()
                file_name = self.dfu.shock_to_file({'file_path': input_directory,
                                                    'shock_id': file['shock_id']
                                                    })['node_file_name']
                file_path = os.path.join(input_directory, file_name)

            # extract the file if it is compressed
            if file_path is not None:
                print("staged input file =" + file_path)
                sys.stdout.flush()
                dfUtil_result = self.dfu.unpack_file({'file_path': file_path})
                file_paths[key] = dfUtil_result['file_path']
            else:
                raise ValueError('No valid files could be extracted based on the input')

        return file_paths

    def _retrieve_gff_file(self, input_gff_file):
        """
        _retrieve_gff_file: retrieve info from gff_file
    
        """
        log("Reading GFF file")
    
        feature_list = collections.defaultdict(list)
        is_patric = 0

        gff_file_handle = open(input_gff_file, 'rb')
        current_line = gff_file_handle.readline()
        line_count = 0

        while (current_line != ''):
            current_line = current_line.strip()

            if(current_line.isspace() or current_line == "" or current_line.startswith("#")):
                pass
            else:
                #Split line
                (contig_id, source_id, feature_type, start, end,
                 score, strand, phase, attributes) = current_line.split('\t')

                #Checking to see if Phytozome
                if "phytozome" in source_id.lower():
                    self.is_phytozome = True

                #Checking to see if Phytozome
                if "PATRIC" in source_id:
                    is_patric = True

                #PATRIC prepends their contig ids with some gibberish
                if is_patric and "|" in contig_id:
                    contig_id = contig_id.split("|", 1)[1]

                #Populating basic feature object
                ftr = {'contig': contig_id, 'source': source_id,
                       'type': feature_type, 'start': int(start),
                       'end': int(end), 'score': score, 'strand': strand,
                       'phase': phase, 'attributes': collections.defaultdict(list)}

                #Populating with attribute key-value pair
                #This is where the feature id is from
                for attribute in attributes.split(";"):
                    attribute = attribute.strip()

                    #Sometimes empty string
                    if not attribute:
                        continue

                    #Use of 1 to limit split as '=' character can also be made available later
                    #Sometimes lack of "=", assume spaces instead
                    if("=" in attribute):
                        key, value = attribute.split("=", 1)
                        ftr['attributes'][key].append(value)
                    elif(" " in attribute):
                        key, value = attribute.split(" ", 1)
                        ftr['attributes'][key].append(value)
                    else:
                        log("Warning: attribute "+attribute+" cannot be separated into key,value pair")

                ftr['attributes']['raw'] = attributes
                if "ID" in ftr['attributes']:
                    ftr['ID'] = ftr['attributes']['ID'][0]
                if "Parent" in ftr['attributes']:
                    ftr['Parent'] = ftr['attributes']['Parent'][0]

                feature_list[contig_id].append(ftr)

            current_line = gff_file_handle.readline()

        gff_file_handle.close()

        #Some GFF/GTF files don't use "ID" so we go through the possibilities        
        feature_list = self._add_missing_identifiers(feature_list)

        #Most bacterial files have only CDSs
        #In order to work with prokaryotic and eukaryotic gene structure synonymously
        #Here we add feature dictionaries representing the parent gene and mRNAs
        feature_list = self._add_missing_parents(feature_list)

        #Phytozome has the annoying habit of editing their identifiers so we fix them
        if self.is_phytozome:
            self._update_phytozome_features(feature_list)

        #All identifiers need to be checked so that they follow the same general rules
        #Rules are listed within the function itself
        feature_list = self._update_identifiers(feature_list)

        return feature_list

    @staticmethod
    def _add_missing_identifiers(feature_list):

        #General rule is to iterate through a range of possibilities if "ID" is missing
        for contig in feature_list.keys():
            for i in range(len(feature_list[contig])):
                if("ID" not in feature_list[contig][i]):
                    for key in ("transcriptId", "proteinId", "PACid", "pacid", "Parent"):
                        if(key in feature_list[contig][i]):
                            feature_list[contig][i]['ID']=feature_list[contig][i][key]
                            break

                    #If the process fails, throw an error
                    for ftr_type in ("gene", "mRNA", "CDS"):
                        if(ftr_type not in feature_list[contig][i]):
                            continue

                        if("ID" not in feature_list[contig][i]):
                            log("Error: Cannot find unique ID to utilize in GFF attributes: "+ \
                                    feature_list[contig][i]['contig']+"."+ \
                                    feature_list[contig][i]['source']+"."+ \
                                    feature_list[contig][i]['type']+": "+ \
                                    feature_list[contig][i]['attributes'])
        return feature_list

    @staticmethod
    def _add_missing_parents(feature_list):

        #General rules is if CDS or RNA missing parent, add them
        for contig in feature_list.keys():
            ftrs = feature_list[contig]
            new_ftrs = []
            for i in range(len(ftrs)):
                if("Parent" not in ftrs[i]):
                    #Assuming parent doesn't exist at all, so create de novo instead of trying to find it
                    if("RNA" in ftrs[i]["type"] or "CDS" in ftrs[i]["type"]):
                        new_gene_ftr = copy.deepcopy(ftrs[i])
                        new_gene_ftr["type"] = "gene"
                        ftrs[i]["Parent"]=new_gene_ftr["ID"]
                        new_ftrs.append(new_gene_ftr)

                    if("CDS" in ftrs[i]["type"]):
                        new_rna_ftr = copy.deepcopy(ftrs[i])
                        new_rna_ftr["type"] = "mRNA"
                        new_ftrs.append(new_rna_ftr)
                        ftrs[i]["Parent"]=new_rna_ftr["ID"]

                new_ftrs.append(ftrs[i])
            feature_list[contig]=new_ftrs
        return feature_list

    @staticmethod
    def _update_phytozome_features(feature_list):

        #General rule is to use the "Name" field where possible
        #And update parent attribute correspondingly
        for contig in feature_list.keys():
            feature_position_dict = {}
            for i in range(len(feature_list[contig])):

                #Maintain old_id for reference
                #Sometimes ID isn't available, so use PACid
                old_id = None
                for key in ("ID", "PACid", "pacid"):
                    if(key in feature_list[contig][i]):
                        old_id = feature_list[contig][i][key]
                        break
                if(old_id is None):
                    #This should be an error
                    print ("Cannot find unique ID, PACid, or pacid in GFF "
                           "attributes: " + feature_list[contig][i][contig])
                    continue

                #Retain old_id
                feature_position_dict[old_id]=i

                #In Phytozome, gene and mRNA have "Name" field, CDS do not
                if("Name" in feature_list[contig][i]):
                    feature_list[contig][i]["ID"] = feature_list[contig][i]["Name"]

                if("Parent" in feature_list[contig][i]):
                    #Update Parent to match new ID of parent ftr
                    feature_list[contig][i]["Parent"] = feature_list[contig][feature_position_dict[feature_list[contig][i]["Parent"]]]["ID"]

        return feature_list

    @staticmethod
    def _update_identifiers(feature_list):

        #General rules:
        #1) Genes keep identifier
        #2) RNAs keep identifier only if its different from gene, otherwise append ".mRNA"
        #3) CDS always uses RNA identifier with ".CDS" appended
        #4) CDS appended with an incremented digit

        CDS_count_dict = dict()
        mRNA_parent_dict = dict()

        for contig in feature_list.keys():
            for ftr in feature_list[contig]:
                if("Parent" in ftr):

                    #Retain old_id of parents
                    old_id = ftr["ID"]

                    if(ftr["ID"] == ftr["Parent"] or "CDS" in ftr["type"]):
                        ftr["ID"] = ftr["Parent"]+"."+ftr["type"]

                    #link old to new ids for mRNA to use with CDS
                    if("RNA" in ftr["type"]):
                        mRNA_parent_dict[old_id]=ftr["ID"]

                    if("CDS" in ftr["type"]):
                        #Increment CDS identifier
                        if(ftr["ID"] not in CDS_count_dict):
                            CDS_count_dict[ftr["ID"]] = 1
                        else:
                            CDS_count_dict[ftr["ID"]] += 1
                        ftr["ID"]=ftr["ID"]+"."+str(CDS_count_dict[ftr["ID"]])

                        #Recall new mRNA id for parent
                        ftr["Parent"]=mRNA_parent_dict[ftr["Parent"]]

        return feature_list

    @staticmethod
    def _print_phytozome_gff(input_gff_file, feature_list):

        #Write modified feature ids to new file
        input_gff_file = input_gff_file.replace("gene", "edited_gene")+".gz"
        try:
            print "Printing to new file: "+input_gff_file
            gff_file_handle = gzip.open(input_gff_file, 'wb')
        except:
            print "Failed to open"

        for contig in sorted(feature_list.iterkeys()):
            for ftr in feature_list[contig]:

                #Re-build attributes
                attributes_dict = {}
                for attribute in ftr["attributes"]['raw'].split(";"):
                    attribute=attribute.strip()

                    #Sometimes empty string
                    if(attribute == ""):
                        continue

                    #Use of 1 to limit split as '=' character can also be made available later
                    #Sometimes lack of "=", assume spaces instead
                    if("=" in attribute):
                        key, value = attribute.split("=", 1)
                    elif(" " in attribute):
                        key, value = attribute.split(" ", 1)
                    else:
                        log("Warning: attribute "+attribute+" cannot be separated into key,value pair")

                    if(ftr[key] != value):
                        value = ftr[key]
                    attributes_dict[key]=value

                ftr["attributes"]=";".join(key+"="+attributes_dict[key] for key in attributes_dict.keys())

                new_line = "\t".join( str(ftr[key]) for key in ['contig', 'source', 'type', 'start', 'end',
                                                                'score', 'strand', 'phase', 'attributes'])
                gff_file_handle.write(new_line)
        gff_file_handle.close()
        return

    def _transform_feature(self, contig, in_feature):
        """Converts a feature from the gff ftr format into the appropriate
        format for a genome object """
        def _aliases(feat):
            keys = ('locus_tag', 'old_locus_tag', 'protein_id',
                    'transcript_id', 'gene', 'EC_number')
            alias_list = []
            for key in keys:
                if key in feat['attributes']:
                    alias_list.extend([(key, val) for val in feat['attributes'][key]])
            return alias_list

        def _propagate_cds_props_to_gene(cds, gene):
            # Put longest protein_translation to gene
            if "protein_translation" not in gene or (
                len(gene["protein_translation"]) <
                    len(cds["protein_translation"])):
                gene["protein_translation"] = cds["protein_translation"]
                gene["protein_translation_length"] = len(
                    cds["protein_translation"])
            # Merge cds list attributes with gene
            for key in ('function', 'aliases', 'db_xref'):
                if cds.get(key, []):
                    gene[key] = cds.get(key, []) + gene.get(key, [])
            # Merge cds["ontology_terms"] -> gene["ontology_terms"]
            terms2 = cds.get("ontology_terms")
            if terms2 is not None:
                terms = gene.get("ontology_terms")
                if terms is None:
                    gene["ontology_terms"] = terms2
                else:
                    for source in terms2:
                        if source in terms:
                            terms[source].update(terms2[source])
                        else:
                            terms[source] = terms2[source]

        # UTRs don't have any information not implied by the exon
        excluded_features = ('five_prime_UTR', 'three_prime_UTR')
        if in_feature['type'] in excluded_features:
            self.skiped_features[in_feature['type']] += 1
            return

        # if the feature ID is duplicated we only need to update the location
        if in_feature['ID'] in self.feature_dict:
            self.feature_dict[in_feature['ID']]['location'].append(
                self._location(in_feature))

        self.feature_counts[in_feature['type']] += 1
        if in_feature['start'] < 1 or in_feature['end'] > len(contig):
            self.warn("Skipping feature with invalid location for specified "
                      "contig: " + str(in_feature))
            return

        feat_seq = contig.seq[in_feature['start']-1:in_feature['end']-1]
        if in_feature['strand'] == '-':
            feat_seq = feat_seq.reverse_complement()

        # The following is common to all the feature types
        out_feat = {
            "id": in_feature['ID'],
            "type": in_feature['type'],
            "location": [self._location(in_feature)],
            "dna_sequence": str(feat_seq),
            "dna_sequence_length": len(feat_seq),
            "md5": hashlib.md5(str(feat_seq)).hexdigest(),
            "function": in_feature['attributes'].get('function', []),
            'warnings': []
        }
        # TODO: Flags, inference_data
        # add optional fields
        if 'note' in in_feature['attributes']:
            out_feat['note'] = in_feature['attributes']["note"][0]
        ont = self._get_ontology(in_feature)
        if ont:
            out_feat['ontology_terms'] = ont
        aliases = _aliases(in_feature)
        if aliases:
            out_feat['aliases'] = aliases
        if 'db_xref' in in_feature['attributes']:
            out_feat['db_xref'] = [tuple(x.split(":")) for x in
                                   in_feature['attributes']['db_xref']]
        if 'product' in in_feature['attributes']:
            out_feat['function'] += ["product:" + x for x in
                                     in_feature['attributes']["product"]]
        parent_id = in_feature.get('Parent', '')
        if parent_id and parent_id not in self.feature_dict:
            raise ValueError("Parent ID: {} was not found in feature ID list.")

        # add type specific features
        if in_feature['type'] == 'gene':
            out_feat['protein_translation_length'] = 0
            out_feat['mrnas'] = []
            out_feat['cdss'] = []

        elif in_feature['type'] == 'CDS':
            prot_seq = str(feat_seq.translate()).strip("*")
            out_feat.update({
                "protein_translation": prot_seq,
                "protein_md5": hashlib.md5(prot_seq).hexdigest(),
                "protein_translation_length": len(prot_seq),
                'parent_gene': "",
            })
            if parent_id:
                # TODO: add location checks and warnings
                parent = self.feature_dict[parent_id]
                if 'cdss' in parent:  # parent must be a gene
                    parent['cdss'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent_id
                else:  # parent must be mRNA
                    parent['cds'] = in_feature['ID']
                    out_feat['parent_mrna'] = parent_id
                    parent_gene = self.feature_dict[parent['parent_gene']]
                    parent_gene['cdss'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent['parent_gene']
                    _propagate_cds_props_to_gene(out_feat, parent_gene)

        elif in_feature['type'] == 'mRNA':
            if parent_id:
                # TODO: add location checks and warnings
                parent = self.feature_dict[parent_id]
                if 'cdss' in parent:  # parent must be a gene
                    parent['mrnas'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent_id

        else:
            out_feat["type"] = in_feature['type']
            if parent_id:
                # TODO: add location checks and warnings
                parent = self.feature_dict[parent_id]
                if 'children' not in parent:
                    parent['children'] = []
                parent['children'].append(out_feat['id'])
                out_feat['parent_gene'] = parent_id

        self.feature_dict[out_feat['id']] = out_feat

    def _gen_genome_info(self, core_genome_name, scientific_name, assembly_ref,
                         source, assembly, input_gff_file, molecule_type):
        """
        _gen_genome_info: generate genome info

        """
        genome = dict()
        genome["id"] = core_genome_name
        genome["scientific_name"] = scientific_name
        genome["assembly_ref"] = assembly_ref
        genome['molecule_type'] = molecule_type
        genome["features"] = []
        genome["cdss"] = []
        genome["mrnas"] = []
        genome['non_coding_features'] = []
        genome["gc_content"] = assembly["gc_content"]
        genome["dna_size"] = assembly["dna_size"]
        genome['md5'] = assembly['md5']
        genome['contig_ids'], genome['contig_lengths'] = zip(
            *[(k, v['length']) for k, v in assembly['contigs'].items()])
        genome['num_contigs'] = len(assembly['contigs'])
        genome["ontology_events"] = [{
            "method": "GenomeFileUtils Genbank uploader from annotations",
            "method_version": self.version,
            "timestamp": self.time_string,
            # TODO: remove this hardcoding
            "id": "GO",
            "ontology_ref": "KBaseOntology/gene_ontology"
        }]
        genome['ontologies_present'] = dict(self.ontologies_present)
        genome['feature_counts'] = dict(self.feature_counts)

        genome['taxonomy'], genome['taxon_ref'], genome['domain'], \
            genome["genetic_code"] = self.gi.retrieve_taxon(self.taxon_wsname,
                                                            genome['scientific_name'])
        genome['source'], genome['genome_tiers'] = self.gi.determine_tier(
            source)

        # Phytozome gff files are not compatible with the RNASeq Pipeline
        # so it's better to build from the object than cache the file
        if not self.is_phytozome:
            gff_file_to_shock = self.dfu.file_to_shock(
                {'file_path': input_gff_file, 'make_handle': 1, 'pack': "gzip"})
            genome['gff_handle_ref'] = gff_file_to_shock['handle']['hid']

        # sort features into their respective arrays
        for feature in self.feature_dict.values():
            if feature['type'] == 'CDS':
                del feature['type']
                genome['cdss'].append(feature)
            elif feature['type'] == 'mRNA':
                if 'exons' not in feature:
                    feature['warnings'].append(
                        'This uploader requires exons to be present to '
                        'establish mRNA location.')
                else:
                    del feature['exons']
                del feature['type']
                genome['mrnas'].append(feature)
            elif feature['type'] == 'gene':
                if genome['cdss']:
                    del feature['type']
                    self.feature_counts["protein_encoding_gene"] += 1
                    genome['features'].append(feature)
                else:
                    del feature['mrnas'], feature['cdss']
                    self.feature_counts["non-protein_encoding_gene"] += 1
                    genome['non_coding_features'].append(feature)
            else:
                if 'exons' in feature:
                    del feature['exons']
                genome['non_coding_features'].append(feature)
        genome['warnings'] = self._warnings

        return genome
