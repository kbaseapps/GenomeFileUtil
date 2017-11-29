# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests  # noqa: F401
import inspect
import StringIO
import urllib
import shutil


from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.authclient import KBaseAuth as _KBaseAuth

from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.core.GenomeInterface import GenomeInterface


class SaveGenomeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        cls.user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': cls.user_id,
                        'provenance': [
                            {'service': 'kb_ke_util',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.shockURL = cls.cfg['shock-url']
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.dfu = DataFileUtil(cls.callback_url)
        cls.cfg['KB_AUTH_TOKEN'] = cls.token

        # build genome interface instance
        gi_config = SDKConfig(cls.cfg)
        cls.genome_interface = GenomeInterface(gi_config)

        # second user
        test_cfg_file = '/kb/module/work/test.cfg'
        test_cfg_text = "[test]\n"
        with open(test_cfg_file, "r") as f:
            test_cfg_text += f.read()

        config = ConfigParser()
        config.readfp(StringIO.StringIO(test_cfg_text))

        test_cfg_dict = dict(config.items("test"))
        if ('test_token2' not in test_cfg_dict):
            raise ValueError("Configuration in <module>/test_local/test.cfg file should " +
                             "include second user credentials ('test_token2')")
        token2 = test_cfg_dict['test_token2']
        user2 = auth_client.get_user(token2)
        cls.ctx2 = MethodContext(None)
        cls.ctx2.update({'token': token2,
                         'user_id': user2,
                         'provenance': [
                            {'service': 'NarrativeService',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                         'authenticated': 1})

        suffix = int(time.time() * 1000)
        cls.wsName = "test_SaveGenomeTest_" + str(suffix)
        cls.wsClient.create_workspace({'workspace': cls.wsName})

        cls.nodes_to_delete = []
        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'nodes_to_delete'):
            for node in cls.nodes_to_delete:
                cls.delete_shock_node(node)

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    @classmethod
    def prepare_data(cls):
        assembly_file_path = os.path.join(cls.scratch,
                                          'e_coli_assembly.fasta')
        shutil.copy('data/e_coli/e_coli_assembly.fasta', assembly_file_path)
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        assembly_ref = au.save_assembly_from_fasta({
            'workspace_name': cls.wsName,
            'assembly_name': 'e_coli.assembly',
            'file': {'path': assembly_file_path}
        })
        cls.test_genome_data = json.load(open('data/e_coli/e_coli.json'))
        cls.test_genome_data['assembly_ref'] = assembly_ref

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def fail_save_one_genome(self, params, error, exception=ValueError, contains=False):
        with self.assertRaises(exception) as context:
            self.getImpl().save_one_genome(self.ctx, params)
        if contains:
            self.assertIn(error, str(context.exception.message))
        else:
            self.assertEqual(error, str(context.exception.message))

    def check_save_one_genome_output(self, ret, genome_name):
        self.assertTrue('info' in ret)

        genome_info = ret['info']
        self.assertEqual(genome_info[1], genome_name)
        self.assertEqual(genome_info[2].split('-')[0].split('.')[1], 'Genome')
        self.assertEqual(genome_info[5], self.user_id)

    def test_bad_one_genome_params(self):
        self.start_test()
        invalidate_params = {'missing_workspace': 'workspace',
                             'name': 'name',
                             'data': 'data'}
        error_msg = '"workspace" parameter is required, but missing'
        self.fail_save_one_genome(invalidate_params, error_msg)

    def test_one_genome(self):
        self.start_test()
        genome_name = 'test_genome'
        params = {'workspace': self.wsName,
                  'name': genome_name,
                  'data': self.test_genome_data}
        ret = self.getImpl().save_one_genome(self.ctx, params)[0]
        self.check_save_one_genome_output(ret, genome_name)

    def test_one_genome_with_hidden(self):
        self.start_test()
        genome_name = 'test_genome_hidden'
        params = {'workspace': self.wsName,
                  'name': genome_name,
                  'data': self.test_genome_data,
                  'hidden': 1}
        ret = self.getImpl().save_one_genome(self.ctx, params)[0]
        self.check_save_one_genome_output(ret, genome_name)

        params = {'workspace': self.wsName,
                  'name': genome_name,
                  'data': self.test_genome_data,
                  'hidden': True}
        ret = self.getImpl().save_one_genome(self.ctx, params)[0]
        self.check_save_one_genome_output(ret, genome_name)

    def test_GenomeInterface_check_dna_sequence_in_features(self):
        # no feature in genome
        genome = {'missing_features': 'features'}
        copied_genome = genome.copy()
        self.genome_interface._check_dna_sequence_in_features(copied_genome)
        self.assertItemsEqual(copied_genome, genome)

        # with contigs
        copied_genome = self.test_genome_data.copy()
        for feat in copied_genome['features']:
            if 'dna_sequence' in feat:
                del feat['dna_sequence']
        self.genome_interface._check_dna_sequence_in_features(copied_genome)

        feature_dna_sum = 0
        for feature in copied_genome['features']:
            if 'dna_sequence' in feature:
                feature_dna_sum += len(feature['dna_sequence'])

        origin_feature_dna_sum = 0
        for feature in self.test_genome_data['features']:
            if 'dna_sequence' in feature:
                origin_feature_dna_sum += len(feature['dna_sequence'])

        self.assertTrue(feature_dna_sum > 3000000)
        self.assertItemsEqual(copied_genome, self.test_genome_data)

    def test_GenomeInterface_own_handle(self):
        # no handle in genome
        genome = {'missing_genbank_handle_ref': 'hid'}
        origin_genome = genome.copy()
        self.genome_interface._own_handle(genome, 'genbank_handle_ref')
        self.assertItemsEqual(origin_genome, genome)

        # user unauthorized
        temp_shock_file = "/kb/module/work/tmp/shock1.txt"
        with open(temp_shock_file, "w") as f1:
            f1.write("Test Shock Handle")
        token2 = self.ctx2['token']
        dfu2 = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=token2)
        shock_ret = dfu2.file_to_shock({'file_path': temp_shock_file, 'make_handle': 1})
        self.nodes_to_delete.append(shock_ret['shock_id'])
        hid = shock_ret['handle']['hid']

        genome = {'genbank_handle_ref': hid}
        with self.assertRaisesRegexp(ValueError, 
                                     'Error getting ACLs for Shock node'):
            self.genome_interface._own_handle(genome, 'genbank_handle_ref')

        # same user
        shock_ret = self.dfu.file_to_shock({'file_path': temp_shock_file, 'make_handle': 1})
        self.nodes_to_delete.append(shock_ret['shock_id'])
        hid = shock_ret['handle']['hid']
        genome = {'genbank_handle_ref': hid}
        origin_genome = genome.copy()
        self.genome_interface._own_handle(genome, 'genbank_handle_ref')
        self.assertDictEqual(origin_genome, genome)

        # differet user
        self.wsClient.set_permissions({'workspace': self.wsName, 'new_permission': 'w',
                                       'users': [self.ctx2['user_id']]})

        token2 = self.ctx2['token']
        dfu2 = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=token2)
        shock_ret = dfu2.file_to_shock({'file_path': temp_shock_file, 'make_handle': 1})
        node = shock_ret['shock_id']
        self.nodes_to_delete.append(node)
        hid = shock_ret['handle']['hid']
    
        # grant user1 read access to node
        user1 = self.ctx['user_id']
        acl = 'read'
        url = self.shockURL + '/node/' + node + '/acl'
        url += '/' + acl + '?users=' + urllib.quote(user1)
        auth_header = {'Authorization': 'OAuth {}'.format(token2)}
        req = requests.put(url, headers=auth_header, allow_redirects=True)
        if not req.ok:
            err = json.loads(req.content)['error'][0]
            print 'response error: {}'.format(err)

        genome = {'genbank_handle_ref': hid}
        origin_genome = genome.copy()
        self.genome_interface._own_handle(genome, 'genbank_handle_ref')
        self.assertNotEqual(origin_genome['genbank_handle_ref'], 
                            genome['genbank_handle_ref'])
