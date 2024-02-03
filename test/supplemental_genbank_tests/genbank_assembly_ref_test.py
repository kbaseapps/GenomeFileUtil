import os
import time
import unittest
from configparser import ConfigParser

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeUtils import warnings
from installed_clients.WorkspaceClient import Workspace as workspaceService


# 3 tests
# 1) Happy case, make sure it kept the reference
# 2) Extra contigs
# 3) Changed sequence

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
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        ws_obj_name = 'ecoli_genome'
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        result = cls.serviceImpl.genbank_to_genome(
            cls.ctx,
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': cls.wsName,
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1,
              'source': "RefSeq Reference"
            })[0]
#        print("HERE IS THE RESULT:")
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=cls.ctx['token'],
                                service_ver='dev')
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        cls.assembly_ref = genome["assembly_ref"]
#        print("GENE 1: ")
#        pprint(cls.genome['features'][0])
#        pprint(result)
        

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_same_genome(self):
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        ws_obj_name = 'ecoli_genome'
        existing_assembly_ref = genome = self.__class__.assembly_ref
        result = self.getImpl().genbank_to_genome(
            self.getContext(),
            {
                'file': {
                    'path': gbk_path},
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name,
                'generate_ids_if_needed': 1,
                'generate_missing_genes': 1,
                'source': 'refseq reference',
                'use_existing_assembly' : existing_assembly_ref
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        self.assertTrue(genome['assembly_ref'] == existing_assembly_ref, "Same file did not keep the same assembly ref")   

    def test_diff_contigs(self):
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_altered_genomic.gbff"
        ws_obj_name = 'ecoli_genome'
        existing_assembly_ref = self.__class__.assembly_ref
        
        with self.assertRaisesRegex(
                ValueError, 
                warnings['assembly_ref_extra_contigs'].format("NC_000913_2.3")):
                self.getImpl().genbank_to_genome(
                                self.getContext(),
                                {
                                    'file': {
                                        'path': gbk_path},
                                    'workspace_name': self.getWsName(),
                                    'genome_name': ws_obj_name,
                                    'generate_ids_if_needed': 1,
                                    'generate_missing_genes': 1,
                                    'source': 'refseq reference',
                                    'use_existing_assembly' : existing_assembly_ref
                                })

    def test_diff_sequence(self):
        gbk_path = "data/e_coli/test_base_change.gbff"
        ws_obj_name = 'ecoli_genome'
        existing_assembly_ref = self.__class__.assembly_ref    
        with self.assertRaisesRegex(
                ValueError, 
                warnings['assembly_ref_diff_seq'].format("NC_000913.3")):
                self.getImpl().genbank_to_genome(
                                self.getContext(),
                                {
                                    'file': {
                                        'path': gbk_path},
                                    'workspace_name': self.getWsName(),
                                    'genome_name': ws_obj_name,
                                    'generate_ids_if_needed': 1,
                                    'generate_missing_genes': 1,
                                    'source': 'refseq reference',
                                    'use_existing_assembly' : existing_assembly_ref
                                })
