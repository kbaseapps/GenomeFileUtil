import os
import time
import unittest
from configparser import ConfigParser
from uuid import uuid4

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService

# NOTE: These tests must run against https://ci.kbase.us
_WORKSPACE_NAME = 'KBaseTestData'
_OBJECT_NAME = 'KBase_derived_GCF_002287175.1.gbff_genome'


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
        cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.callbackURL = os.environ['SDK_CALLBACK_URL']
        cls.handleURL = cls.cfg['handle-service-url']

    def test_update_taxon_assignments_valid(self):
        """
        Test a valid call to the update_taxon_assignments method.
        """
        taxon_key = str(uuid4())
        taxon_val = str(uuid4())
        taxon_val_new = str(uuid4())
        # Copy the object to test workspace
        dfu = DataFileUtil(self.callbackURL)
        obj_ref = f"{_WORKSPACE_NAME}/{_OBJECT_NAME}"
        result = dfu.get_objects({'object_refs': [obj_ref]})['data'][0]
        obj_data = result['data']
        # crate user owned handle in the object and update it
        hs = HandleService(self.handleURL)
        prev_handle_id = obj_data['genbank_handle_ref']
        prev_shock_id = hs.hids_to_handles([prev_handle_id])[0]['id']
        new_handle_id = dfu.own_shock_node({'shock_id': prev_shock_id, 'make_handle': 1})['handle']['hid']
        obj_data['genbank_handle_ref'] = new_handle_id
        # Save new object in test workspace
        obj_info = result['info']
        new_obj = {
            'type': obj_info[2],
            'data': obj_data,
            'name': 'GCF_002287175.1'
        }
        test_ws_id = dfu.ws_name_to_id(self.wsName)
        infos = dfu.save_objects({
            'id': test_ws_id,
            'objects': [new_obj]
        })
        obj_ref = f"{infos[0][6]}/{infos[0][0]}/{infos[0][4]}"
        new_ws_id = infos[0][6]
        new_obj_id = infos[0][0]
        get_obj_params = {
            'wsid': new_ws_id,
            'objid': new_obj_id,
            'included': ['/taxon_assignments']
        }
        # Add a new assignment
        self.serviceImpl.update_taxon_assignments(self.ctx, {
            'workspace_id': new_ws_id,
            'object_id': new_obj_id,
            'taxon_assignments': {
                taxon_key: taxon_val
            }
        })
        # Fetch the object and check the mapping
        obj = self.wsClient.get_objects2({'objects': [get_obj_params]})['data'][0]['data']
        self.assertTrue(taxon_key in obj['taxon_assignments'])
        self.assertEqual(obj['taxon_assignments'][taxon_key], taxon_val)
        # Update the assignment we just added
        self.serviceImpl.update_taxon_assignments(self.ctx, {
            'workspace_id': new_ws_id,
            'object_id': new_obj_id,
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
            'workspace_id': new_ws_id,
            'object_id': new_obj_id,
            'remove_assignments': [taxon_key]
        })
        # Fetch the object and check the mapping
        obj = self.wsClient.get_objects2({'objects': [get_obj_params]})['data'][0]['data']
        self.assertTrue(taxon_key not in obj['taxon_assignments'])
        self.assertEqual(obj['taxon_assignments'].get(taxon_key), None)
