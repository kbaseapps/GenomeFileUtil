import os
import pytest
import re
import shutil
import time
import unittest
from configparser import ConfigParser

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from conftest import assert_exception_correct


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
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        cls.wsID = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_genome_upload_fungi(self):
        gbk_path = "data/fungi/GCF_000002945.1_ASM294v2_genomic.gbff.gz"
        ws_obj_name = 'FungiGenome.1'
        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': self.wsName,
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        self.assertTrue(int(
            result['genome_info'][10]['Number of Protein Encoding Genes']) > 0)

    def test_genome_upload_two_files(self):
        gbk_path = "data/drosophila/small_test.tar.gz"
        ws_obj_name = 'DrosophilaGenome.1'
        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
              'file': {'path': gbk_path},
              'workspace_name': self.wsName,
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        self.assertTrue(int(
            result['genome_info'][10]['Number of Protein Encoding Genes']) > 0)

    def test_unable_to_find_locus_tag(self):
        gbk_path = "data/drosophila/small_test.tar.gz"
        ws_obj_name = 'DrosophilaGenome.1'
        with self.assertRaisesRegex(
                        ValueError, 
                        "Unable to find a valid id for gene among these tags: locus_tag," + \
                        " kbase_id. Correct the file or rerun with generate_ids"):
            self.serviceImpl.genbank_to_genome(
                self.ctx,
                {
                    'file': {'path': gbk_path},
                    'workspace_name': self.wsName,
                    'genome_name': ws_obj_name,
                    'generate_ids_if_needed': 0
            })

    def test_feature_id_duplication_bug(self):
        gbk_path = "data/duplication.gbff"
        ws_obj_name = 'BugGenome.1'
        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
                'file': {'path': gbk_path},
                'workspace_name': self.wsName,
                'genome_name': ws_obj_name,
                'generate_ids_if_needed': 1,
                'generate_missing_genes': 1
            })[0]
        self.assertEqual(result['genome_info'][10]['Number contigs'], '1')
        self.assertEqual(result['genome_info'][10]['Number of Protein Encoding Genes'], '2')

    def test_upload_prokka(self):
        gbk_path = "data/prokka/PROKKA_012345.gbk.gz"
        ws_obj_name = 'CustomGenome.1'
        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
                'file': {'path': gbk_path},
                'workspace_name': self.wsName,
                'genome_name': ws_obj_name,
                'generate_ids_if_needed': 1,
                "generate_missing_genes": 1
            })[0]
        self.assertEqual(result['genome_info'][10]['Number contigs'], '1')
        self.assertEqual(result['genome_info'][10]['Number of Protein Encoding Genes'], '4313')
        self.assertEqual(result['genome_info'][10]['Domain'], 'Unknown')

    def test_ftp_upload_bug(self):
        gbk_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/589/275/GCF_001589275.1_sce3192.1/GCF_001589275.1_sce3192.1_genomic.gbff.gz"
        ws_obj_name = 'BugGenome.2'
        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
                'file': {'ftp_url': gbk_url },
                'workspace_name': self.wsName,
                'genome_name': ws_obj_name,
                'generate_ids_if_needed': 1
            })[0]
        self.assertTrue(int(
            result['genome_info'][10]['Number of Protein Encoding Genes']) > 0)

    def test_simple_shock_upload(self):
        ### Test for upload from SHOCK - upload the file to shock first
        print('attempting upload through shock')
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        shutil.copy(gbk_path, self.cfg['scratch'])
        shock_id = data_file_cli.file_to_shock({
            'file_path': os.path.join(self.cfg['scratch'], gbk_path.split("/")[-1])
        })['shock_id']
        print("Running test")
        ws_obj_name2 = 'MyGenome.2'
        result = self.serviceImpl.genbank_to_genome(self.ctx, {
                'file': {'shock_id': shock_id},
                'workspace_name': self.wsName,
                'genome_name': ws_obj_name2,
            })[0]
        self.assertIsNotNone(result['genome_ref'])
        self.assertTrue(int(
            result['genome_info'][10]['Number of Protein Encoding Genes']) > 0)

    def _run_test_fail(self, params, error_message):
        with pytest.raises(Exception) as got:
            self.serviceImpl.genbank_to_genome(self.ctx, params)
        assert_exception_correct(got.value, ValueError(error_message))

    def _run_test_fail_mass(self, params, error_message):
        with pytest.raises(Exception) as got:
            self.serviceImpl.genbanks_to_genomes(self.ctx, params)
        assert_exception_correct(got.value, ValueError(error_message))

    def _check_result_object_info_fields(self, results, file_names, object_metas):
        object_version_pattern = re.compile(r'^[0-9]+\/1$')
        for idx, res in enumerate(results):
            assert object_version_pattern.match("/".join(res['genome_ref'].split("/")[-2:]))
            # Check relevant object info fields
            obj = self.wsClient.get_object_info3(
                {'objects': [{'ref': res['genome_ref'], 'includeMetadata': 1}]}
            )
            info = obj['infos'][0]
            assert info == res['genome_info']
            assert info[1] == file_names[idx]
            assert info[2].split('-')[0] == 'KBaseGenomes.Genome'
            assert info[6] == self.wsID
            print(info[10])
            assert info[10] == object_metas[idx]

    def test_genbank_to_genome_invalid_workspace(self):
        genome_name = "GCF_000970165.1_ASM97016v1_genomic.gbff.gz"
        params = {
            "file": {'path': f"data/gbff/{genome_name}"},
            "genome_name": genome_name,
        }
        self._run_test_fail(
            params,
            "Exactly one of a 'workspace_id' or a 'workspace_name' parameter must be provided",
        )

    def test_genbank_to_genome_fail_both_workspace_id_and_workspace_name(self):
        genome_name = "GCF_000970165.1_ASM97016v1_genomic.gbff.gz"
        params = {
            "workspace_id": self.wsID,
            "workspace_name": self.wsName,
            "file": {'path': f"data/gbff/{genome_name}"},
            "genome_name": genome_name,
        }
        self._run_test_fail(
            params,
            "Exactly one of a 'workspace_id' or a 'workspace_name' parameter must be provided",
        )

    def test_genbanks_to_genomes(self):
        genome_name1 = "GCF_000970165.1_ASM97016v1_genomic.gbff.gz"
        genome_name2 = "GCF_000970185.1_ASM97018v1_genomic.gbff.gz"
        
        file_names = [genome_name1, genome_name2]
        object_metas = {}
        
        results = self.serviceImpl.genbanks_to_genomes(
            self.ctx,
            {
                "workspace_id": self.wsID,
                "inputs": [
                    {
                        "file": {"path": f"data/gbff/{genome_name1}"},
                        "genome_name": genome_name1,
                        "metadata": {"foo": "bar"},
                    },
                    {
                        "file": {"path": f"data/gbff/{genome_name2}"},
                        "genome_name": genome_name2,
                        "metadata": {"bar": "foo"},
                    }
                ]
            }
        )[0]['results']

        self.assertEqual(len(results), 2)
        self._check_result_object_info_fields(results, file_names, object_metas)

    def test_genbanks_to_genomes_invalid_workspace_id(self):
        genome_name = "GCF_000970165.1_ASM97016v1_genomic.gbff.gz"
        params = {
            "inputs": [
                {
                    'file': {'path': f"data/gbff/{genome_name}"},
                    'genome_name': genome_name,
                }
            ]
        }
        self._run_test_fail_mass(params, "workspace_id is required")

        wsid = float(self.wsID)
        params = {
            "workspace_id": wsid,
            "inputs": [
                {
                    'file': {'path': f"data/gbff/{genome_name}"},
                    'genome_name': genome_name,
                }
            ]
        }
        self._run_test_fail_mass(
            params, f"workspace_id must be an integer, got: {wsid}"
        )

        params = {
            "workspace_id": 0,
            "inputs": [
                {
                    'file': {'path': f"data/gbff/{genome_name}"},
                    'genome_name': genome_name,
                }
            ]
        }
        self._run_test_fail_mass(params, "workspace_id must be an integer >= 1")

    def test_genbanks_to_genomes_invalid_inputs(self):
        params = {"workspace_id": self.wsID, "inputs": []}
        self._run_test_fail_mass(
            params, "inputs field is required and must be a non-empty list"
        )

        params = {
            "workspace_id": self.wsID,
            "inputs": [["genome_file_path", "genome_name"]]
        }
        self._run_test_fail_mass(
            params,
            "Entry #1 in inputs field is not a mapping as required",
        )

    def test_genbanks_to_genomes_invalid_params(self):
        e = 'required "genome_name" field was not defined'
        params = {
            "workspace_id": self.wsID,
            "inputs": [
                {"genome": "GCF_000970165.1_ASM97016v1_genomic.gbff.gz"}
            ]
        }
        self._run_test_fail_mass(
            params,
            f"Entry #1 in inputs field has invalid params: {e}",
        )
