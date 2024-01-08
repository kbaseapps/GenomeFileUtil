import os
import pytest
import unittest

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from conftest import assert_exception_correct


class ImplTest(unittest.TestCase):
    def setUp(self):
        self.cfg = {}
        self.cfg["workspace-url"] = "https://ci.kbase.us/services/ws"
        self.cfg["shock-url"] = "https://ci.kbase.us/services/shock-api"
        self.cfg["handle-service-url"] = "https://ci.kbase.us/services/handle_service"
        self.cfg["scratch"] = "/kb/module/work/tmp"
        self.cfg["srv-wiz-url"] = "https://ci.kbase.us/services/service_wizard"
        self.cfg["auth-service-url"] = "https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login"
        self.cfg["re-api-url"] = "https://ci.kbase.us/services/relation_engine_api"

    def tearDown(self):
        for kbase_secure_param in (
            "KBASE_SECURE_CONFIG_PARAM_MAX_THREADS",
            "KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU",
        ):
            os.environ.pop(kbase_secure_param, None)

    def _run_test_fail(self, cfg, error_message):
        with pytest.raises(Exception) as got:
            GenomeFileUtil(cfg)
        assert_exception_correct(got.value, ValueError(error_message))

    def test_valid_catalog_param_type(self):
        os.environ["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10"
        os.environ["KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU"] = "2.5"
        gfu = GenomeFileUtil(self.cfg)
        self.assertEqual(gfu.max_threads, 10)
        self.assertEqual(gfu.threads_per_cpu, 2.5)

    def test_valid_default_catalog_param_type(self):
        gfu = GenomeFileUtil(self.cfg)
        self.assertEqual(gfu.max_threads, 10)
        self.assertEqual(gfu.threads_per_cpu, 1)

    def test_invalid_max_threads(self):
        os.environ["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10.5"
        os.environ["KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU"] = "2.5"
        self._run_test_fail(self.cfg, "MAX_THREADS must be of type int")

    def test_invalid_threads_per_cpu(self):
        os.environ["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10"
        os.environ["KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU"] = "2.8e"
        self._run_test_fail(self.cfg, "THREADS_PER_CPU must be of type float")
