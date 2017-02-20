
import os
import sys
import shutil
import traceback
import uuid
import urllib2

from urlparse import urlparse
from pprint import pprint, pformat

from biokbase.workspace.client import Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil

from GenomeFileUtil.core.FastaGFFUploaderScript import upload_genome

class FastaGFFToGenome:

    def __init__(self, sdk_config):
        self.cfg = sdk_config

    def import_file(self, ctx, params):

        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) construct the input directory staging area

        input_directory =  os.path.join(self.cfg.sharedFolder, 'genome-upload-staging-'+str(uuid.uuid4()))
        os.makedirs(input_directory)

        file_paths = self.stage_input(params,input_directory)
        pprint(file_paths)

        # 3) extract out the parameters
        parsed_params = self.set_defaults()

        # required parameters must be defined
        parsed_params['workspace_name'] = params['workspace_name']
        parsed_params['genome_name'] = params['genome_name']

        # add any optional parameters
        optional_param_fields_to_check = [
                'source',
                'taxon_wsname',
                'taxon_reference',
                'release',
                'genetic_code',
#                'exclude_ontologies',
                'type',
#                'metadata'
            ]

        for field in optional_param_fields_to_check:
            if field in params:
                parsed_params[field] = params[field]

        # 4) Do the upload
        pprint(parsed_params)
        result = upload_genome(
                shock_service_url = self.cfg.shockURL,
                handle_service_url = self.cfg.handleURL,
                workspace_service_url = self.cfg.workspaceURL,
                callback_url = self.cfg.callbackURL,

                input_fasta_file=file_paths["fasta_file"],
                input_gff_file=file_paths["gff_file"],
        
                workspace_name   = parsed_params['workspace_name'],
                core_genome_name = parsed_params['genome_name'],

#                taxon_wsname     = parsed_params['taxon_wsname'],
#                taxon_lookup_obj_name = parsed_params['taxon_lookup_obj_name'],
#                taxon_reference = parsed_params['taxon_reference'],

                source           = parsed_params['source'],
                genome_type             = parsed_params['type'],
#                release          = parsed_params['release'],
#                genetic_code     = parsed_params['genetic_code'],

#                exclude_ontologies = parsed_params['exclude_ontologies'],
#                ontology_wsname = parsed_params['ontology_wsname'],
#                ontology_GO_obj_name = parsed_params['ontology_GO_obj_name'],
#                ontology_PO_obj_name = parsed_params['ontology_PO_obj_name'],

#                provenance = ctx['provenance'],
#                usermeta = parsed_params['metadata']
            )

        # 5) clear the temp directory
        shutil.rmtree(input_directory)

        # 6) return the result
        info = result['genome_info']
        details = {
            'genome_ref': str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
            'genome_info': info
#            'report_name': result['report_name'],
#            'report_ref': result['report_ref']
        }

        return details

    def validate_params(self, params):
        for key in ('workspace_name','genome_name','fasta_file','gff_file'):
            if key not in params:
                raise ValueError('required "'+key+'" field was not defined')

        # one and only one of 'path', or 'shock_id' is required
        for key in ('fasta_file','gff_file'):
            file = params[key]
            if not isinstance(file, dict):
                raise ValueError('required "'+key+'" field must be a map/dict')
            n_valid_fields = 0
            if 'path' in file and file['path'] is not None:
                n_valid_fields += 1
            if 'shock_id' in file and file['shock_id'] is not None:
                n_valid_fields += 1
            if 'ftp_url' in file and file['ftp_url'] is not None:
                n_valid_fields += 1
            if n_valid_fields < 1:
                raise ValueError('required "'+key+'" field must include one source: path | shock_id | ftp_url')
            if n_valid_fields > 1:
                raise ValueError('required "'+key+'" field has too many sources specified: ' + str(file.keys()))

        valid_types = ['Reference','User upload','Representative']
        if 'type' in params and params['type'] not in valid_types:
            raise ValueError('Entered value for type is not one of the valid entries of "Reference", "Representative" or "User upload"')

    def set_defaults(self):
        default_params = {
            'source':'Genbank',
            'taxon_wsname': self.cfg.raw['taxon-workspace-name'],
            'taxon_lookup_obj_name': self.cfg.raw['taxon-lookup-object-name'],
            'taxon_reference': None,

            'ontology_wsname': self.cfg.raw['ontology-workspace-name'],
            'ontology_GO_obj_name': self.cfg.raw['ontology-gene-ontology-obj-name'],
            'ontology_PO_obj_name': self.cfg.raw['ontology-plant-ontology-obj-name'],

            'release': None,
            'genetic_code': None,
            'generate_ids_if_needed': 0,
            'exclude_ontologies': 0,
            'type': 'User upload',
            'metadata':{}
        }
        return default_params

    def stage_input(self, params, input_directory):
        ''' Setup the input_directory by fetching the files and uncompressing if needed. '''

        # at this point, the 'fasta_file' and 'gff_file' input is validated, so we don't have to catch any special cases
        # we expect one and only one of path, or shock_id

        # determine how to get each file: if it is from shock, download it.  If it
        # is just sitting there, then use it.  Move the files to the staging input directory
        file_paths  = dict()
        for key in ('fasta_file','gff_file'):
            file = params[key]
            file_path = None
            if 'path' in file and file['path'] is not None:
                # copy the local file to the input staging directory
                # (NOTE: could just move it, but then this method would have the side effect of moving your
                # file which another SDK module might have an open handle on)
                local_file_path = file['path']
                file_path = os.path.join(input_directory, os.path.basename(local_file_path))
                shutil.copy2(local_file_path, file_path)

            if 'shock_id' in file and file['shock_id'] is not None:
                # handle shock file
                print('Downloading file from SHOCK node: ' + str(self.cfg.shockURL) + ' - ' + str(file['shock_id']))
                sys.stdout.flush()
                dfUtil = DataFileUtil(self.cfg.callbackURL)
                file_name = dfUtil.shock_to_file({'file_path': input_directory,
                                                  'shock_id': file['shock_id']
                                                  })['node_file_name']
                file_path = os.path.join(input_directory, file_name)

            if 'ftp_url' in file and file['ftp_url'] is not None:
                # Note that the Transform originally had a script_utils.download_from_urls method
                # that, if the url is a folder, pulls all subfiles.  That code recently broke when
                # fetching from NCBI (not clear if it is our issue or NCBI), but for now just
                # support the most common case- an FTP to a single file.
                print('Downloading file from: ' + str(file['ftp_url']))
                sys.stdout.flush()

                url = urlparse(file['ftp_url'])
                if url.scheme != 'ftp' and url.scheme != 'http':
                    raise ValueError('Only FTP/HTTP servers are supported')
                file_name = 'assembly.fasta'
                if(key == 'gff_file'):
                    file_name = 'features.gff'
                if url.path != '':
                    file_name = url.path.split('/')[-1]

                req = urllib2.Request(file['ftp_url'])
                response = urllib2.urlopen(req)
                file_data = response.read()

                file_path = os.path.join(input_directory, file_name)
                with open(file_path, "w") as file_handle:
                    file_handle.write(file_data)

            # extract the file if it is compressed
            if file_path is not None:
                print("staged input file =" + file_path)
                sys.stdout.flush()
                dfUtil = DataFileUtil(self.cfg.callbackURL)
                dfUtil_result = dfUtil.unpack_file({ 'file_path': file_path })
                file_paths[key]=dfUtil_result['file_path']
            else:
                raise ValueError('No valid files could be extracted based on the input')

        return file_paths
