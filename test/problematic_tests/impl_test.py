import os
import pytest
import traceback
import unittest

from configparser import ConfigParser

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil


def assert_exception_correct(got: Exception, expected: Exception):
    """
    Compare raised exception with expected exception.

    Args:
        got (Exception): Exception received
        expected (Exception): Exception expected
    """
    err = "".join(traceback.TracebackException.from_exception(got).format())
    assert got.args == expected.args, err
    assert type(got) == type(expected)


class ImplTest(unittest.TestCase):
    def setUp(self):
        config_file = os.environ["KB_DEPLOYMENT_CONFIG"]
        token = os.environ.get("KB_AUTH_TOKEN", None)

        self.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items("GenomeFileUtil"):
            self.cfg[nameval[0]] = nameval[1]
        self.cfg["KB_AUTH_TOKEN"] = token

    def test_invalid_catalog_param_type(self):
        MAX_THREADS = "MAX_THREADS"
        THREADS_PER_CPU = "THREADS_PER_CPU"

        # max_threads type check fails
        os.environ["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10.5"
        os.environ["KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU"] = "2.5"

        with pytest.raises(Exception) as got:
            GenomeFileUtil(self.cfg)
        assert_exception_correct(
            got.value, ValueError(f"{MAX_THREADS} must be of type {int.__name__}")
        )

        # threads_per_cpu type check fails
        os.environ["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10"
        os.environ["KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU"] = "2.8e"

        with pytest.raises(Exception) as got:
            GenomeFileUtil(self.cfg)
        assert_exception_correct(
            got.value, ValueError(f"{THREADS_PER_CPU} must be of type {float.__name__}")
        )
