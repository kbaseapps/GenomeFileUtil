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
from installed_clients.WorkspaceClient import Workspace as workspaceService
from conftest import assert_exception_correct


class GenomeFileUtilTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.provenance = [
            {
                'service': 'GenomeFileUtil',
                'method': 'run_local_tests',
                'method_params': [],
            }
        ]
        cls.ctx.update(
            {
                'token': token,
                'provenance': cls.provenance,
                'authenticated': 1
            }
        )
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
        cls.dfuClient = DataFileUtil(os.environ['SDK_CALLBACK_URL'])

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
        shutil.copy(gbk_path, self.cfg['scratch'])
        shock_id = self.dfuClient.file_to_shock({
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

    def _check_result_object_info_fields_and_provenance(
        self,
        results,
        file_names,
        object_metas,
        expected_provenance,
    ):
        object_version_pattern = re.compile(r'^[0-9]+\/1$')
        for idx, res in enumerate(results):
            assert object_version_pattern.match("/".join(res['genome_ref'].split("/")[-2:]))
            data = self.wsClient.get_objects2({"objects": [{'ref': res['genome_ref']}]})["data"][0]

            # check info
            info = data["info"]
            assert info == res['genome_info']
            assert info[1] == file_names[idx]
            assert info[2].split('-')[0] == 'KBaseGenomes.Genome'
            assert info[6] == self.wsID
            assert all(info[10].get(k) == v for k, v in object_metas[idx].items())

            # check provenance
            provenance = data["provenance"][0]
            retrieved_provenance = [
                {
                    'service': provenance['service'],
                    'method': provenance['method'],
                    'method_params': provenance['method_params'],
                }
            ]

            assert retrieved_provenance == expected_provenance

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

    def test_genbank_to_genome(self):
        genome_name = "GCF_000970205.1_ASM97020v1_genomic.gbff.gz"
        object_metas = [
            {
                "GC content": "0.41421",
                "Size": "4142816",
                "Number contigs": "1",
                "MD5": "cf47d74f66a16dffcbaa7a05eb9eec70",
                "temp": "curr",
            }
        ]

        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
                "workspace_id": self.wsID,
                "file": {"path": f"data/gbff/{genome_name}"},
                "genome_name": genome_name,
                "metadata": {"temp": "curr"},
            })

        self._check_result_object_info_fields_and_provenance(
            result, [genome_name], object_metas, self.provenance
        )

    def test_genbanks_to_genomes(self):
        genome_name1 = "GCF_000970165.1_ASM97016v1_genomic.gbff.gz"
        genome_name2 = "GCF_000970185.1_ASM97018v1_genomic.gbff.gz"
        
        file_names = [genome_name1, genome_name2]
        object_metas = [
            {
                "GC content": "0.41457",
                "Size": "4096482",
                "Number contigs": "1",
                "MD5": "949a0fe665048cb917c8cf74f75c74b7",
                "foo": "bar",
            },
            {
                "GC content": "0.41487",
                "Size": "4066551",
                "Number contigs": "1",
                "MD5": "d33802829ba0686714a5d74280527615",
                "bar": "foo",
            }
        ]

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
        self._check_result_object_info_fields_and_provenance(
            results, file_names, object_metas, self.provenance
        )

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

    # TODO test other params in genbanks_to_genomes function
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

    def test_genbanks_to_genomes_generate_ids_and_missing_genes(self):
        genome_name = "Cyanidioschyzon_merolae_one_locus.gbff"

        object_metas = [
            {
                "GC content": "0.27065",
                "Size": "32211",
                "Number contigs": "1",
                "MD5": "43b94ee0851f3b9e9db521167c6fcba3",
                "curr": "temp",
            }
        ]

        results = self.serviceImpl.genbanks_to_genomes(
            self.ctx,
            {
                "workspace_id": self.wsID,
                "inputs": [
                    {
                        "file": {"path": f"data/Cyanidioschyzon/{genome_name}"},
                        "genome_name": genome_name,
                        "generate_ids_if_needed": 1,
                        "generate_missing_genes": 1,
                        "metadata": {"curr": "temp"},
                    }
                ]
            }
        )[0]['results']

        # TODO maybe more test
        self._check_result_object_info_fields_and_provenance(
            results, [genome_name], object_metas, self.provenance
        )