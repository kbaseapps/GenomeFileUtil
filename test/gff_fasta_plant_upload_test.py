import unittest
import os
import time
import shutil
import json

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext

from DataFileUtil.DataFileUtilClient import DataFileUtil

from GenomeFileUtil.JsonIOHelper import (download_genome_to_json_files, 
                                         compare_genome_json_files)


class GenomeFileUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
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
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

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
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_simple_upload(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        fasta_path = "data/Test_Plant/Test_v1.0.fa.gz"
        gff_path = "data/Test_Plant/Test_v1.0.gene.gff3.gz"
        ws_obj_name = 'MyGenome'
        ws_name = self.getWsName()

        ### Test for a Local Function Call
        print('attempting upload via local function directly')

        result = genomeFileUtil.fasta_gff_to_genome(self.getContext(), 
            {
                'fasta_file' : { 'path': fasta_path },
                'gff_file' : { 'path' : gff_path },
                'workspace_name':ws_name,
                'genome_name':ws_obj_name
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        #target_dir = os.path.join("/kb/module/work/tmp", "Test_Plant")
        #download_genome_to_json_files(self.getContext()['token'], result['genome_ref'],
        #                              target_dir)
        #self.assertEqual(0, len(compare_genome_json_files(target_dir, 
        #                                                  os.path.join("/kb/module/test/data", 
        #                                                               "GCF_000005845"))))
        # todo: add test that result is correct

        ### Test for upload from SHOCK - upload the file to shock first
        #print('attempting upload through shock')
        #data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
        #                        token=self.__class__.ctx['token'],
        #                        service_ver='dev')
        #shutil.copy(gbk_path, self.__class__.cfg['scratch'])
        #shock_id = data_file_cli.file_to_shock({
        #    'file_path': os.path.join(self.__class__.cfg['scratch'], gbk_path.split("/")[-1])
        #})['shock_id']
        #ws_obj_name2 = 'MyGenome.2'
        #result2 = genomeFileUtil.genbank_to_genome(self.getContext(), 
        #    {
        #        'file': {'shock_id':shock_id},
        #        'workspace_name':self.getWsName(),
        #        'genome_name':ws_obj_name2,
        #    })[0]
        #pprint(result2)
        #self.assertIsNotNone(result['genome_ref'])
        # todo: add test that result is correct
