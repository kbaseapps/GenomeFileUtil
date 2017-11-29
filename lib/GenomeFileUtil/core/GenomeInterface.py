
import time
import requests
import json
import re
import sys

from Workspace.WorkspaceClient import Workspace as Workspace
from GenomeFileUtil.authclient import KBaseAuth as _KBaseAuth
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService  # @UnresolvedImport @IgnorePep8
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblySequenceAPI.AssemblySequenceAPIServiceClient import AssemblySequenceAPI
from collections import defaultdict


def log(message, prefix_newline=False):
    time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
    print(('\n' if prefix_newline else '') + time_str + ': ' + message)


class GenomeInterface:

    def __init__(self, config):
        self.ws_url = config.workspaceURL
        self.handle_url = config.handleURL
        self.shock_url = config.shockURL
        self.sw_url = config.srvWizURL
        self.token = config.token
        self.auth_service_url = config.authServiceUrl
        self.callback_url = config.callbackURL

        self.ws = Workspace(self.ws_url, token=self.token)
        self.auth_client = _KBaseAuth(self.auth_service_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.taxon_wsname = config.raw['taxon-workspace-name']

    @staticmethod
    def _validate_save_one_genome_params(params):
        """
        _validate_save_one_genome_params:
                validates params passed to save_one_genome method
        """

        log('start validating save_one_genome params')

        # check for required parameters
        for p in ['workspace', 'name', 'data']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _check_shock_response(self, response, errtxt):
        """
        _check_shock_response: check shock node response (Copied from DataFileUtil)
        """
        log('start checking shock response')

        if not response.ok:
            try:
                err = json.loads(response.content)['error'][0]
            except:
                # this means shock is down or not responding.
                self.log("Couldn't parse response error content from Shock: " +
                         response.content)
                response.raise_for_status()
            raise ValueError(errtxt + str(err))

    def _own_handle(self, genome_data, handle_property):
        """
        _own_handle: check that handle_property point to shock nodes owned by calling user
        """

        log('start checking handle {} ownership'.format(handle_property))

        if handle_property in genome_data:
            handle_id = genome_data[handle_property]
            hs = HandleService(self.handle_url, token=self.token)
            handles = hs.hids_to_handles([handle_id])
            shock_id = handles[0]['id']

            # Copy from DataFileUtil.own_shock_node implementation:
            header = {'Authorization': 'Oauth {}'.format(self.token)}
            res = requests.get(self.shock_url + '/node/' + shock_id +
                               '/acl/?verbosity=full',
                               headers=header, allow_redirects=True)
            self._check_shock_response(
                res, 'Error getting ACLs for Shock node {}: '.format(shock_id))
            owner = res.json()['data']['owner']['username']
            user_id = self.auth_client.get_user(self.token)

            if owner != user_id:
                log('start copying node to owner: {}'.format(user_id))
                dfu_shock = self.dfu.copy_shock_node({'shock_id': shock_id,
                                                      'make_handle': True})
                handle_id = dfu_shock['handle']['hid']
                genome_data[handle_property] = handle_id

    def _check_dna_sequence_in_features(self, genome):
        """
        _check_dna_sequence_in_features: check dna sequence in each feature
        """
        log('start checking dna sequence in each feature')

        if 'features' in genome:
            features_to_work = {}
            for feature in genome['features']:
                if not ('dna_sequence' in feature and feature['dna_sequence']):
                    features_to_work[feature['id']] = feature['location']

            if len(features_to_work) > 0:
                aseq = AssemblySequenceAPI(self.sw_url, token=self.token)
                get_dna_params = {'requested_features': features_to_work}
                if 'assembly_ref' in genome:
                    get_dna_params['assembly_ref'] = genome['assembly_ref']
                elif 'contigset_ref' in genome:
                    get_dna_params['contigset_ref'] = genome['contigset_ref']
                else:
                    # Nothing to do (it may be test genome without contigs)...
                    return
                dna_sequences = aseq.get_dna_sequences(get_dna_params)['dna_sequences']
                for feature in genome['features']:
                    if feature['id'] in dna_sequences:
                        feature['dna_sequence'] = dna_sequences[feature['id']]
                        feature['dna_sequence_length'] = len(feature['dna_sequence'])
        
    def save_one_genome(self, params):
        log('start saving genome object')

        self._validate_save_one_genome_params(params)

        workspace = params['workspace']
        name = params['name']
        data = params['data']
        if 'meta' in params and params['meta']:
            meta = params['meta']
        else:
            meta = {}
        if 'feature_counts' not in data:
            data = self._update_genome(data)

        # check all handles point to shock nodes owned by calling user
        self._own_handle(data, 'genbank_handle_ref')
        self._own_handle(data, 'gff_handle_ref')

        self._check_dna_sequence_in_features(data)
        data['warnings'] = self.validate_genome(data)

        if 'hidden' in params and str(params['hidden']).lower() in ('yes', 'true', 't', '1'):
            hidden = 1
        else:
            hidden = 0

        if isinstance(workspace, int) or workspace.isdigit():
            workspace_id = workspace
        else:
            workspace_id = self.dfu.ws_name_to_id(workspace)

        dfu_save_params = {'id': workspace_id,
                           'objects': [{'type': 'NewTempGenomes.Genome',
                                        'data': data,
                                        'name': name,
                                        'meta': meta,
                                        'hidden': hidden}]}

        dfu_oi = self.dfu.save_objects(dfu_save_params)[0]

        returnVal = {'info': dfu_oi, 'warnings': data['warnings']}

        return returnVal

    def retrieve_taxon(self, taxon_wsname, scientific_name):
        """
        _retrieve_taxon: retrieve taxonomy and taxon_reference

        """
        default = ('Unconfirmed Organism: '+ scientific_name, 'ReferenceTaxons/unknown_taxon', 'Unknown', 11)
        solr_url = 'http://kbase.us/internal/solr-ci/search/'
        solr_core = 'taxonomy_ci'
        query = '/select?q=scientific_name:"{}"&fl=scientific_name%2Cscientific_lineage%2Ctaxonomy_id%2Cdomain%2Cgenetic_code&rows=5&wt=json'
        match = re.match("\S+\s?\S*", scientific_name)
        if not match:
            return default
        res = requests.get(solr_url+solr_core+query.format(match.group(0)))
        results = res.json()['response']['docs']
        if not results:
            return default
        taxonomy = results[0]['scientific_lineage']
        taxon_reference = '{}/{}_taxon'.format(
            taxon_wsname, results[0]['taxonomy_id'])
        domain = results[0]['domain']
        genetic_code = results[0]['genetic_code']

        return taxonomy, taxon_reference, domain, genetic_code

    @staticmethod
    def determine_tier(source):
        """
        Given a user provided source parameter, assign a source and genome tier
        """
        low_source = source.lower()
        if 'refseq' in low_source:
            if 'reference' in low_source:
                return "Refseq", ['Reference', 'Representative',
                                  'ExternalDB']
            if 'representative' in low_source:
                return "Refseq", ['Representative', 'ExternalDB']
            return "Refseq", ['ExternalDB']
        if 'phytozome' in low_source:
            if 'flagship' in source:
                return "Phytosome", ['Reference', 'Representative',
                                     'ExternalDB']
            return "Phytosome", ['Representative', 'ExternalDB']
        if 'ensembl' in low_source:
            return "Ensembl", ['Representative', 'ExternalDB']
        return source, ['User']

    def _update_genome(self, genome):
        """Checks for missing required fields and fixes breaking changes"""
        # do top level updates
        if 'genome_tier' not in genome:
            genome['source'], genome['genome_tiers'] = self.determine_tier(
                genome['source'])
        if 'molecule_type' not in genome:
            genome['molecule_type'] = 'Unknown'
        if 'taxon_ref' not in genome:
            genome['taxonomy'], genome['taxon_ref'], genome['domain'], \
                genome['genetic_code'] = self.retrieve_taxon(
                    self.taxon_wsname, genome['scientific_name'])

        if any([x not in genome for x in ('dna_size', 'md5', 'gc_content')]):
            assembly_data = self.dfu.get_objects(
                {'object_refs': [genome['assembly_ref']],
                 'ignore_errors': 0})['data'][0]['data']
            genome["gc_content"] = assembly_data['gc_content']
            genome["dna_size"] = assembly_data['dna_size']
            genome["md5"] = assembly_data['md5']

        if 'cdss' not in genome:
            genome['cdss'] = []
        if 'mrnas' not in genome:
            genome['mrnas'] = []

        # do feature level updates
        retained_features = []
        type_counts = defaultdict(int)
        for field in ('mrnas', 'cdss', 'features'):
            for i, feat in enumerate(genome.get(field, [])):
                if 'function' in feat and not isinstance(feat, list):
                    feat['function'] = [feat['function']]
                if 'aliases' in feat:
                    feat['aliases'] = [['db_xref', x] for x in feat['aliases']]
                if 'type' in feat:
                    type_counts[feat['type']] += 1
                # TODO: Ontologies

                # split all the stuff lumped together in old versions into the
                # right arrays
                if field == 'features':
                    if feat.get('type', 'gene') == 'gene':
                        if not feat.get('cdss', []):
                            genome['non_coding_features'].append(feat)
                        else:
                            retained_features.append(feat)
                    elif feat.get('type', 'gene') == 'CDS':
                        if 'protein_md5' not in feat:
                            feat['protein_md5'] = ''
                        if 'parent_gene' not in feat:
                            feat['parent_gene'] = ''
                        genome['cdss'].append(feat)
                    elif feat.get('type', 'gene') == 'mRNA':
                        genome['mrnas'].append(feat)

        genome['features'] = retained_features

        type_counts['mRNA'] = len(genome.get('mrnas', []))
        type_counts['CDS'] = len(genome.get('cdss', []))
        type_counts['protein_encoding_gene'] = len(genome['features'])
        type_counts['non-protein_encoding_gene'] = len(genome.get('non_coding_features', []))
        genome['feature_counts'] = type_counts

        return genome

    @staticmethod
    def validate_genome(g, print_size=True):
        """
        Run a series of checks on the genome object and return any warnings
        """
        def _get_size(obj):
            return sys.getsizeof(json.dumps(obj))
        allowed_tiers = {'Representative', 'Reference', 'ExternalDB', 'User'}

        log('Validating genome object contents')
        warnings = []
        if len(g.get('cdss', [])) < len(g['features']):
            warnings.append("CDS array should be at at least as long as the "
                            "Features array.")

        if g['domain'] == "Bacteria" and len(g.get('cdss', [])) != len(g['features']):
            warnings.append("For prokaryotes, CDS array should be the same "
                            "length as the Features array.")

        if g['domain'] == "Eukaryota" and len(g.get('mrnas', [])) and \
                len(g.get('mrnas', [])) == len(g.get('cdss', [])):
            warnings.append("For Eukaryotes, CDS array should not be the same "
                            "length as the Features array due to RNA splicing.")

        if g.get('mrnas', []) and len(g.get('mrnas', [])) != len(g.get('cdss', [])):
            warnings.append("mRNA array should be the same length as the CDS"
                            "array if present.")

        if "molecule_type" in g and g['molecule_type'] not in {"DNA", 'ds-DNA'}:
            if g.get('domain', '') not in {'Virus', 'Viroid'} and \
                            g['molecule_type'] not in {"DNA", 'ds-DNA'}:
                warnings.append("Genome molecule_type {} is not expected "
                                "for domain {}.".format(g['molecule_type'],
                                                        g.get('domain', '')))

        if "genome_tiers" in g and set(g['genome_tiers']) - allowed_tiers:
            warnings.append("Undefined terms in genome_tiers: " + ", ".join(
                set(g['genome_tiers']) - allowed_tiers))

        if print_size:
            print("Subobject Sizes:")
            for x in ('cdss', 'mrnas', 'features', 'non_coding_features',
                      'ontology_present'):
                if x in g:
                    print("{}: {} bytes".format(x, _get_size(g[x])))
            print("Total size {} bytes".format(_get_size(g)))
        return warnings
