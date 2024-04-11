import json
import os
import pytest
import re
import shutil
import time
import unittest
from configparser import ConfigParser
from datetime import datetime

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.WorkspaceClient import Workspace as workspaceService
from conftest import assert_exception_correct

_GENOMEFILEUTIL_VERSION = "local-dev"
_DATAFILEUTIL_VERSION = "0.2.1-release"
_ASSEMBLYUTIL_VERSION = "3.1.1-release"
_GENOMEANNOTATIONAPI_VERSION = "1.0.2-release"
_WSLARGEDATAIO_VERSION = "0.0.5-beta"

def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    if isinstance(obj, int):
        return str(obj)
    else:
        return obj

class GenomeFileUtilTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.provenance = [
            {
                "service": "GenomeFileUtil",
                "service_ver": "local-dev",
                "method": "run_local_tests",
                "method_params": [],
                "input_ws_objects": [],
                "resolved_ws_objects": [],
                "intermediate_incoming": [],
                "intermediate_outgoing": [],
                "external_data": [],
                "subactions": [
                    {
                        "name": "GenomeFileUtil",
                        "ver": _GENOMEFILEUTIL_VERSION,
                        "code_url": "https://localhost",
                    },
                    {
                        "name": "AssemblyUtil",
                        "ver": _ASSEMBLYUTIL_VERSION,
                        "code_url": "https://github.com/kbaseapps/AssemblyUtil",
                    },
                    {
                        "name": "DataFileUtil",
                        "ver": _DATAFILEUTIL_VERSION,
                        "code_url": "https://github.com/kbaseapps/DataFileUtil",
                    },
                    {
                        'name': 'GenomeAnnotationAPI',
                        'ver': _GENOMEANNOTATIONAPI_VERSION,
                        'code_url': 'https://github.com/kbase/genome_annotation_api',
                    },
                    {
                        "name": "WsLargeDataIO",
                        "ver": _WSLARGEDATAIO_VERSION,
                        "code_url": "https://github.com/kbaseapps/WsLargeDataIO",
                    },
                ],
                "custom": {},
                "description": "KBase SDK method run via the KBase Execution Engine",
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

    def _load_expected_data(self, json_path):
        with open(json_path, "r") as read_file:
            data = json.load(read_file)
        return data

    def _dump_retrieved_data(self, json_path, dictionary):
        with open(json_path, "w") as outfile:
            json.dump(dictionary, outfile)

    def _retrieve_provenance(self, provenance):
        # make a copy to avoid modifying the original provenance
        provs = [prov.copy() for prov in provenance]
        for prov in provs:
            for key in ["time", "epoch"]:
                prov.pop(key)
            for subaction in prov['subactions']:
                subaction.pop("commit")
        return provs

    def _retrieve_genome_metadata(self, metadata):
        # make a copy to avoid modifying the original metadata
        metadata = dict(metadata)
        metadata.pop("Assembly Object")
        return metadata

    def _retrieve_genome_data(self, data):
        # make a copy to avoid modifying the original data
        data = dict(data)
        for key in ["assembly_ref", "genbank_handle_ref"]:
            data.pop(key)
        for ontology_event in data.get("ontology_events", []):
            ontology_event.pop("timestamp")
        return data

    def _retrieve_assembly_data(self, data):
        # make a copy to avoid modifying the original data
        data = dict(data)
        data.pop("fasta_handle_ref")
        fasta_handle_info = data["fasta_handle_info"]
        fasta_handle_info.pop("shock_id")
        for key in ["hid", "id"]:
            fasta_handle_info["handle"].pop(key)
        return data

    def _check_result_assembly_info_provenance_data(
        self,
        results,
        file_names,
        expected_metadata,
        expected_provenance,
        expected_data,
    ):
        object_version_pattern = re.compile(r'^[0-9]+\/[0-9]+\/1$')
        for idx, res in enumerate(results):

            assert object_version_pattern.match(res['assembly_ref'])
            obj = self.wsClient.get_objects2({"objects": [{'ref': res['assembly_ref']}]})["data"][0]

            # check info
            info = obj["info"]
            assert info == res['assembly_info']
            assert info[1] == file_names[idx] + "_assembly"
            assert info[2].split('-')[0] == 'KBaseGenomeAnnotations.Assembly'
            # datetime.fromisoformat is not available in Python 3.6 or below
            assert datetime.strptime(info[3], '%Y-%m-%dT%H:%M:%S+%f')
            assert info[6] == self.wsID
            assert info[7] == self.wsName
            # check metadata
            retrieved_metadata = info[10]
            assert retrieved_metadata == expected_metadata[idx]

            # check provenance
            provenance = obj["provenance"]
            retrieved_provenance = self._retrieve_provenance(provenance)
            assert retrieved_provenance == expected_provenance

            # check data
            data = obj["data"]
            retrieved_data = self._retrieve_assembly_data(data)
            print("---------------")
            print(f"retrieved data is {retrieved_data}")
            print(f"expected data is {expected_data[idx]}")
            print("---------------")
            assert retrieved_data == expected_data[idx]

    def _check_result_genome_info_provenance_data(
        self,
        results,
        file_names,
        expected_metadata,
        expected_provenance,
        expected_data,
    ):
        object_version_pattern = re.compile(r'^[0-9]+\/[0-9]+\/1$')
        for idx, res in enumerate(results):

            assert object_version_pattern.match(res['genome_ref'])
            obj = self.wsClient.get_objects2({"objects": [{'ref': res['genome_ref']}]})["data"][0]

            # check info
            info = obj["info"]
            assert info == res['genome_info']
            assert info[1] == file_names[idx]
            assert info[2].split('-')[0] == 'KBaseGenomes.Genome'
            # datetime.fromisoformat is not available in Python 3.6 or below
            assert datetime.strptime(info[3], '%Y-%m-%dT%H:%M:%S+%f')
            assert info[6] == self.wsID
            assert info[7] == self.wsName
            # check metadata
            retrieved_metadata = self._retrieve_genome_metadata(info[10])
            assert retrieved_metadata == expected_metadata[idx]

            # check provenance
            provenance = obj["provenance"]
            retrieved_provenance = self._retrieve_provenance(provenance)
            assert retrieved_provenance == expected_provenance

            # check data
            data = obj["data"]
            retrieved_data = self._retrieve_genome_data(data)
            assert ordered(retrieved_data) == ordered(expected_data[idx])

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
                "temp": "curr",
                "Taxonomy": "Unconfirmed Organism",
                "Size": "4142816",
                "Source": "Genbank",
                "Name": "Methanosarcina mazei S-6",
                "GC content": "0.41421",
                "Genetic code": "11",
                "Number of Genome Level Warnings": "1",
                "Source ID": "NZ_CP009512",
                "Number of Protein Encoding Genes": "3509",
                "Number contigs": "1",
                "Domain": "Unknown",
                "Number of CDS": "3509",
                "MD5": "cf47d74f66a16dffcbaa7a05eb9eec70",
            }
        ]

        expected_data = self._load_expected_data("data/genome_curated/genome_GCF_000970205.json")

        result = self.serviceImpl.genbank_to_genome(
            self.ctx,
            {
                "workspace_id": self.wsID,
                "file": {"path": f"data/gbff/{genome_name}"},
                "genome_name": genome_name,
                "metadata": {"temp": "curr"},
            })

        self._check_result_genome_info_provenance_data(
            result, [genome_name], object_metas, self.provenance, [expected_data]
        )

    def test_genbanks_to_genomes(self):
        genome_name1 = "GCF_000970185.1_ASM97018v1_genomic.gbff.gz"
        genome_name2 = "Cyanidioschyzon_merolae_one_locus.gbff"
        genome_name3 = "mRNA_with_no_parent.gbff"
        genome_name4 = "ontology.gbff"

        file_names = [genome_name1, genome_name2, genome_name3, genome_name4]

        genome_metas = [
            {
                "Taxonomy": "Unconfirmed Organism",
                "Size": "4066551",
                "Source": "Genbank",
                "Name": "Methanosarcina mazei SarPi",
                "GC content": "0.41487",
                "Genetic code": "11",
                "bar": "foo",
                "Number of Genome Level Warnings": "1",
                "Source ID": "NZ_CP009511",
                "Number of Protein Encoding Genes": "3403",
                "Number contigs": "1",
                "Domain": "Unknown",
                "Number of CDS": "3403",
                "MD5": "d33802829ba0686714a5d74280527615",
            },
            {
                "Taxonomy": "Unconfirmed Organism",
                "Size": "32211",
                "Source": "Genbank",
                "Name": "Cyanidioschyzon merolae strain 10D",
                "GC content": "0.27065",
                "Genetic code": "11",
                "Number of Genome Level Warnings": "2",
                "Source ID": "MitoASM9120v1",
                "Number of Protein Encoding Genes": "34",
                "Number contigs": "1",
                "curr": "temp",
                "Domain": "Unknown",
                "Number of CDS": "34",
                "MD5": "43b94ee0851f3b9e9db521167c6fcba3",
            },
            {
                "temp": "curr",
                "Taxonomy": "Unconfirmed Organism",
                "Size": "32211",
                "Source": "Genbank",
                "Name": "Cyanidioschyzon merolae strain 10D",
                "GC content": "0.27065",
                "Genetic code": "11",
                "Number of Genome Level Warnings": "1",
                "Source ID": "MitoASM9120v1",
                "Number of Protein Encoding Genes": "0",
                "Number contigs": "1",
                "Domain": "Unknown",
                "Number of CDS": "0",
                "MD5": "43b94ee0851f3b9e9db521167c6fcba3",
            },
            {
                "Taxonomy": "Unconfirmed Organism",
                "Size": "32211",
                "foo": "bar",
                "Source": "Genbank",
                "Name": "Cyanidioschyzon merolae strain 10D",
                "GC content": "0.27065",
                "Genetic code": "11",
                "Number of Genome Level Warnings": "1",
                "Source ID": "MitoASM9120v1",
                "Number of Protein Encoding Genes": "1",
                "Number contigs": "1",
                "Domain": "Unknown",
                "Number of CDS": "1",
                "MD5": "43b94ee0851f3b9e9db521167c6fcba3",
            },
        ]

        assembly_metas = [
            {'GC content': '0.41487', 'Size': '4066551', 'N Contigs': '1', 'MD5': 'd33802829ba0686714a5d74280527615'},
            {'GC content': '0.27065', 'Size': '32211', 'N Contigs': '1', 'MD5': '43b94ee0851f3b9e9db521167c6fcba3'},
            {'GC content': '0.27065', 'Size': '32211', 'N Contigs': '1', 'MD5': '43b94ee0851f3b9e9db521167c6fcba3'},
            {'GC content': '0.27065', 'Size': '32211', 'N Contigs': '1', 'MD5': '43b94ee0851f3b9e9db521167c6fcba3'},
        ]

        expected_genome_data = [
            self._load_expected_data("data/genome_curated/genome_GCF_000970185.json"),
            self._load_expected_data("data/genome_curated/genome_Cyanidioschyzon_merolae_one_locus.json"),
            self._load_expected_data("data/genome_curated/genome_mRNA_with_no_parent.json"),
            self._load_expected_data("data/genome_curated/genome_ontology.json"),
        ]

        expected_assembly_data = [
            self._load_expected_data("data/genome_curated/assembly_GCF_000970185.json"),
            self._load_expected_data("data/genome_curated/assembly_Cyanidioschyzon_merolae_one_locus.json"),
            self._load_expected_data("data/genome_curated/assembly_mRNA_with_no_parent.json"),
            self._load_expected_data("data/genome_curated/assembly_ontology.json"),
        ]

        results = self.serviceImpl.genbanks_to_genomes(
            self.ctx,
            {
                "workspace_id": self.wsID,
                "inputs": [
                    {
                        "file": {"path": f"data/gbff/{genome_name1}"},
                        "genome_name": genome_name1,
                        "metadata": {"bar": "foo"},
                    },
                    {
                        "file": {"path": f"data/Cyanidioschyzon/{genome_name2}"},
                        "genome_name": genome_name2,
                        "generate_ids_if_needed": 1,
                        "generate_missing_genes": 1,
                        "metadata": {"curr": "temp"},
                    },
                    {
                        "file": {"path": f"data/genome_curated/{genome_name3}"},
                        "genome_name": genome_name3,
                        "generate_ids_if_needed": 1,
                        "generate_missing_genes": 1,
                        "metadata": {"temp": "curr"},
                    },
                    {
                        "file": {"path": f"data/genome_curated/{genome_name4}"},
                        "genome_name": genome_name4,
                        "generate_ids_if_needed": 1,
                        "generate_missing_genes": 1,
                        "metadata": {"foo": "bar"},
                    }
                ]
            }
        )[0]['results']

        self.assertEqual(len(results), 4)
        self._check_result_genome_info_provenance_data(
            results,
            file_names,
            genome_metas,
            self.provenance,
            expected_genome_data
        )
        self._check_result_assembly_info_provenance_data(
            results,
            file_names,
            assembly_metas,
            self.provenance,
            expected_assembly_data
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
