# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
import sys
import shutil
import traceback
import uuid
from pprint import pprint, pformat

from GenomeFileUtil.core.GenbankToGenome import GenbankToGenome
from GenomeFileUtil.core.GenomeToGFF import GenomeToGFF
from GenomeFileUtil.core.GenomeToGenbank import GenomeToGenbank
from GenomeFileUtil.core.FastaGFFToGenome import FastaGFFToGenome

from biokbase.workspace.client import Workspace

from DataFileUtil.DataFileUtilClient import DataFileUtil

# Used to store and pass around configuration URLs more easily
class SDKConfig:
    def __init__(self, config):
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.sharedFolder = config['scratch']
        self.raw = config

#END_HEADER


class GenomeFileUtil:
    '''
    Module Name:
    GenomeFileUtil

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.5.5"
    GIT_URL = "ssh://git@github.com/samseaver/GenomeFileUtil"
    GIT_COMMIT_HASH = "f476ab49fdee619be27cd64fa1f78f95052da7c3"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.cfg = SDKConfig(config)
        #END_CONSTRUCTOR
        pass


    def genbank_to_genome(self, ctx, params):
        """
        :param params: instance of type "GenbankToGenomeParams" (genome_name
           - becomes the name of the object workspace_name - the name of the
           workspace it gets saved to. source - Source of the file typically
           something like RefSeq or Ensembl taxon_ws_name - where the
           reference taxons are : ReferenceTaxons taxon_reference - if
           defined, will try to link the Genome to the specified taxonomy
           object insteas of performing the lookup during upload release -
           Release or version number of the data per example Ensembl has
           numbered releases of all their data: Release 31
           generate_ids_if_needed - If field used for feature id is not
           there, generate ids (default behavior is raising an exception)
           genetic_code - Genetic code of organism. Overwrites determined GC
           from taxon object type - Reference, Representative or User upload)
           -> structure: parameter "file" of type "File" -> structure:
           parameter "path" of String, parameter "shock_id" of String,
           parameter "ftp_url" of String, parameter "genome_name" of String,
           parameter "workspace_name" of String, parameter "source" of
           String, parameter "taxon_wsname" of String, parameter
           "taxon_reference" of String, parameter "release" of String,
           parameter "generate_ids_if_needed" of String, parameter
           "genetic_code" of Long, parameter "type" of String, parameter
           "metadata" of type "usermeta" -> mapping from String to String
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "genome_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genbank_to_genome
        print('genbank_to_genome -- paramaters = ')
        pprint(params)

        importer = GenbankToGenome(self.cfg)
        result = importer.import_file(ctx, params)

        print('import complete -- result = ')
        pprint(result)
        #END genbank_to_genome

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genbank_to_genome return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def fasta_gff_to_genome(self, ctx, params):
        """
        :param params: instance of type "FastaGFFToGenomeParams" (genome_name
           - becomes the name of the object workspace_name - the name of the
           workspace it gets saved to. source - Source of the file typically
           something like RefSeq or Ensembl taxon_ws_name - where the
           reference taxons are : ReferenceTaxons taxon_reference - if
           defined, will try to link the Genome to the specified taxonomy
           object insteas of performing the lookup during upload release -
           Release or version number of the data per example Ensembl has
           numbered releases of all their data: Release 31
           generate_ids_if_needed - If field used for feature id is not
           there, generate ids (default behavior is raising an exception)
           genetic_code - Genetic code of organism. Overwrites determined GC
           from taxon object type - Reference, Representative or User upload)
           -> structure: parameter "fasta_file" of type "File" -> structure:
           parameter "path" of String, parameter "shock_id" of String,
           parameter "ftp_url" of String, parameter "gff_file" of type "File"
           -> structure: parameter "path" of String, parameter "shock_id" of
           String, parameter "ftp_url" of String, parameter "genome_name" of
           String, parameter "workspace_name" of String, parameter "source"
           of String, parameter "taxon_wsname" of String, parameter
           "taxon_reference" of String, parameter "release" of String,
           parameter "genetic_code" of Long, parameter "type" of String,
           parameter "metadata" of type "usermeta" -> mapping from String to
           String
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "genome_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN fasta_gff_to_genome
        print('fasta_gff_to_genome -- paramaters = ')
        pprint(params)

        importer = FastaGFFToGenome(self.cfg)
        result = importer.import_file(ctx, params)
        pprint(result)
        #END fasta_gff_to_genome

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method fasta_gff_to_genome return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def genome_to_gff(self, ctx, params):
        """
        :param params: instance of type "GenomeToGFFParams" (is_gtf -
           optional flag switching export to GTF format (default is 0, which
           means GFF) target_dir - optional target directory to create file
           in (default is temporary folder with name 'gff_<timestamp>'
           created in scratch)) -> structure: parameter "genome_ref" of
           String, parameter "ref_path_to_genome" of list of String,
           parameter "is_gtf" of type "boolean" (A boolean - 0 for false, 1
           for true. @range (0, 1)), parameter "target_dir" of String
        :returns: instance of type "GenomeToGFFResult" (from_cache is 1 if
           the file already exists and was just returned, 0 if the file was
           generated during this call.) -> structure: parameter "gff_file" of
           type "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "ftp_url" of String, parameter
           "from_cache" of type "boolean" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_to_gff
        print('genome_to_gff -- paramaters = ')
        pprint(params)

        exporter = GenomeToGFF(self.cfg)
        result = exporter.export(ctx, params)

        print('export complete -- result = ')
        pprint(result)
        #END genome_to_gff

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genome_to_gff return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def genome_to_genbank(self, ctx, params):
        """
        :param params: instance of type "GenomeToGenbankParams" -> structure:
           parameter "genome_ref" of String, parameter "ref_path_to_genome"
           of list of String
        :returns: instance of type "GenomeToGenbankResult" (from_cache is 1
           if the file already exists and was just returned, 0 if the file
           was generated during this call.) -> structure: parameter
           "genbank_file" of type "File" -> structure: parameter "path" of
           String, parameter "shock_id" of String, parameter "ftp_url" of
           String, parameter "from_cache" of type "boolean" (A boolean - 0
           for false, 1 for true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_to_genbank
        print('genome_to_genbank -- paramaters = ')
        pprint(params)

        exporter = GenomeToGenbank(self.cfg)
        result = exporter.export(ctx, params)

        print('export complete -- result = ')
        pprint(result)
        #END genome_to_genbank

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genome_to_genbank return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def export_genome_as_genbank(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_genome_as_genbank
        print('export_genome_as_genbank -- paramaters = ')

        # validate parameters
        if 'input_ref' not in params:
            raise ValueError('Cannot run export_genome_as_genbank- no "input_ref" field defined.')

        # get WS metadata to get ws_name and obj_name
        ws = Workspace(url=self.cfg.workspaceURL)
        info = ws.get_object_info_new({'objects':[{'ref': params['input_ref'] }],'includeMetadata':0, 'ignoreErrors':0})[0]

        genome_to_genbank_params = {
          'genome_ref': params['input_ref']
        }

        # export to file
        result = self.genome_to_genbank(ctx, genome_to_genbank_params)[0]['genbank_file'];

        # create the output directory and move the file there
        export_package_dir = os.path.join(self.cfg.sharedFolder, info[1])
        os.makedirs(export_package_dir)
        shutil.move(
          result['file_path'],
          os.path.join(export_package_dir, os.path.basename(result['file_path'])))

        # package it up and be done
        dfUtil = DataFileUtil(self.cfg.callbackURL)
        package_details = dfUtil.package_for_download({
                                    'file_path': export_package_dir,
                                    'ws_refs': [ params['input_ref'] ]
                                })

        output = { 'shock_id': package_details['shock_id'] }

        print('export complete -- result = ')
        pprint(output)
        #END export_genome_as_genbank

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_genome_as_genbank return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
