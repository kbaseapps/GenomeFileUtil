import mock
import os
import time
import unittest
from configparser import ConfigParser

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from installed_clients.WorkspaceClient import Workspace as workspaceService

class GenomeFileUtilTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'provenance': [
                            {'service': 'GenomeFileUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.dfuClient = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=token, service_ver='dev')
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsID = cls.dfuClient.ws_name_to_id(cls.wsName)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def _prep_input(self, gbk_path, ws_obj_name, g2g_mass=False):

        base_input = {
            'file': {'path': gbk_path},
            'workspace_name': self.wsName,
            'genome_name': ws_obj_name,
            'generate_ids_if_needed': 1
        }

        if g2g_mass:

            base_input.pop('workspace_name')
            return {
                'workspace_id': self.wsID,
                'inputs': [base_input]
            }

        return base_input

    def test_full_sequence(self):
        # features should not have sequences in it. But both non_coding_features and CDSs should have sequences.
        print("test_full_sequence")
        self._test_sequence(
            "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff",
            "full_sequence",
            self._check_full_sequence
        )

    def _test_sequence(self, gbk_path, ws_obj_name, check_seq_func, g2g_mass=False):
        input_data = self._prep_input(gbk_path, ws_obj_name)
        if g2g_mass:
            result = self.serviceImpl.genbanks_to_genomes(self.ctx, input_data)[0]
        else:
            result = self.serviceImpl.genbank_to_genome(self.ctx, input_data)[0]

        check_seq_func(result)

    def _check_full_sequence(self, result):
        genome = self.dfuClient.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        count_features_without_dna_sequence = 0
        for feature in genome['features']:
            if "dna_sequence" not in feature:
                count_features_without_dna_sequence += 1
        count_non_coding_features_without_sequence = 0
        for feature in genome['non_coding_features']:
            if "dna_sequence" not in feature:
                if feature["dna_sequence_length"] <= 10000:
                    count_non_coding_features_without_sequence += 1
                    print("non_coding_feature_without_sequence: " + str(feature))
        count_cdss_without_sequence = 0
        for feature in genome['cdss']:
            if "dna_sequence" not in feature:
                count_cdss_without_sequence += 1        
        self.assertTrue(count_features_without_dna_sequence == 0,"All features should have DNA sequences.")
        self.assertTrue(count_non_coding_features_without_sequence == 0,
                        "All non_coding_features should have DNA sequences.")
        self.assertTrue(count_cdss_without_sequence == 0,"All CDSs should have DNA sequences.")

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 14000000)
    def test_partial_sequence(self):
        # features should not have sequences in it. But both non_coding_features and CDSs should have sequences.
        print("test_partial_sequence")
        self._test_sequence(
            "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff",
            "partial_sequence",
            self._check_partial_sequence
        )

    def _check_partial_sequence(self, result):
        genome = self.dfuClient.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        count_features_with_dna_sequence = 0
        for feature in genome['features']:
            if "dna_sequence" in feature:
                count_features_with_dna_sequence += 1
        count_non_coding_features_without_sequence = 0
        for feature in genome['non_coding_features']:
            if "dna_sequence" not in feature:
                if feature["dna_sequence_length"] <= 10000:
                    count_non_coding_features_without_sequence += 1
                    print("non_coding_feature_without_sequence: " + str(feature))
        count_cdss_without_sequence = 0
        for feature in genome['cdss']:
            if "dna_sequence" not in feature:
                count_cdss_without_sequence += 1        
        self.assertTrue(count_features_with_dna_sequence == 0,"All features should not have DNA sequences.")
        self.assertTrue(count_non_coding_features_without_sequence == 0,
                        "All non_coding_features should have DNA sequences.")
        self.assertTrue(count_cdss_without_sequence == 0,"All CDSs should have DNA sequences.")

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 9000000)
    def test_no_sequence_kept(self):
        # features, cds, and non_coding_features should not have sequences in it.
        print("test_no_sequence_kept")
        self._test_sequence(
            "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff",
            "no_sequence",
            self._check_no_sequence
        )

    def _check_no_sequence(self, result):
        genome = self.dfuClient.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        count_features_with_dna_sequence = 0
        for feature in genome['features']:
            if "dna_sequence" in feature:
                count_features_with_dna_sequence += 1
        count_non_coding_features_with_sequence = 0
        for feature in genome['non_coding_features']:
            if "dna_sequence" in feature:
                count_non_coding_features_with_sequence += 1
        count_cdss_with_sequence = 0
        for feature in genome['cdss']:
            if "dna_sequence" in feature:
                count_cdss_with_sequence += 1        
        self.assertTrue(count_features_with_dna_sequence == 0,"All features should not have DNA sequences.")
        self.assertTrue(count_non_coding_features_with_sequence == 0,
                        "All non_coding_features should not have DNA sequences.")
        self.assertTrue(count_cdss_with_sequence == 0,"All CDSs should not have DNA sequences.")

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 1)
    def test_max_genome_size(self):
        with self.assertRaisesRegex(ValueError, "This genome size of "):
            GenomeInterface.validate_genome({"taxon_ref": "", "domain": ""})

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 1)
    def test_max_genome_size_with_genbanks_to_genomes(self):
        print("test_max_genome_size_with_genbanks_to_genomes")
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        ws_obj_name = 'large_genome_with_g2g_mass'
        input_data = self._prep_input(gbk_path, ws_obj_name, g2g_mass=True)
        with self.assertRaisesRegex(ValueError, "This genome size of "):
            self.serviceImpl.genbanks_to_genomes(self.ctx, input_data)
