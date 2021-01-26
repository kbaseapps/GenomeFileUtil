import os
import time
import unittest
from configparser import ConfigParser
from uuid import uuid4

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.WorkspaceClient import Workspace as workspaceService


class UpdateTaxonAssignmentsTest(unittest.TestCase):

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
        config_file = os.environ['KB_DEPLOYMENT_CONFIG']
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsID = ws_info[0]
        # Copy the example genome into the test workspace
        obj_info = cls.wsClient.copy_object({
            "from": {"ref": "33192/33/17"},
            "to": {"wsid": cls.wsID, "name": "test_taxonomy_assignment_genome"},
        })
        cls.example_obj_id = obj_info[0]
        cls.example_obj_wsid = obj_info[6]

    def test_update_taxon_assignments_valid(self):
        """
        Test a valid call to the update_taxon_assignments method.
        """
        taxon_key = str(uuid4())
        taxon_val = str(uuid4())
        taxon_val_new = str(uuid4())
        # Add a new assignment
        self.serviceImpl.update_taxon_assignments(self.ctx, {
            'workspace_id': self.example_obj_wsid,
            'object_id': self.example_obj_id,
            'taxon_assignments': {
                taxon_key: taxon_val
            }
        })
        # Fetch the object and check the mapping
        get_obj_params = {
            'wsid': self.example_obj_wsid,
            'objid': self.example_obj_id,
            'included': ['/taxon_assignments']
        }
        obj = self.wsClient.get_objects2({'objects': [get_obj_params]})['data'][0]['data']
        self.assertTrue(taxon_key in obj['taxon_assignments'])
        print(obj['taxon_assignments'])
        self.assertEqual(obj['taxon_assignments'][taxon_key], taxon_val)
        # Update the assignment we just added
        self.serviceImpl.update_taxon_assignments(self.ctx, {
            'workspace_id': self.example_obj_wsid,
            'object_id': self.example_obj_id,
            'taxon_assignments': {
                taxon_key: taxon_val_new
            }
        })
        # Fetch the object and check the mapping
        obj = self.wsClient.get_objects2({'objects': [get_obj_params]})['data'][0]['data']
        self.assertTrue(taxon_key in obj['taxon_assignments'])
        self.assertEqual(obj['taxon_assignments'][taxon_key], taxon_val_new)
        # Remove the assignment we just added
        self.serviceImpl.update_taxon_assignments(self.ctx, {
            'workspace_id': self.example_obj_wsid,
            'object_id': self.example_obj_id,
            'remove_assignments': [taxon_key]
        })
        # Fetch the object and check the mapping
        obj = self.wsClient.get_objects2({'objects': [get_obj_params]})['data'][0]['data']
        self.assertTrue(taxon_key not in obj['taxon_assignments'])
        self.assertEqual(obj['taxon_assignments'].get(taxon_key), None)
