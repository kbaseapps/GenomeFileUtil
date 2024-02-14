# -*- coding: utf-8 -*-
#BEGIN_HEADER

import json
import logging
import os
import shutil
from pprint import pprint

from GenomeFileUtil.core.FastaGFFToGenome import FastaGFFToGenome
from GenomeFileUtil.core.GenbankToGenome import (
    GenbankToGenome,
    MAX_THREADS_DEFAULT,
    THREADS_PER_CPU_DEFAULT,
)
from GenomeFileUtil.core.GenomeFeaturesToFasta import GenomeFeaturesToFasta
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from GenomeFileUtil.core.GenomeToGFF import GenomeToGFF
from GenomeFileUtil.core.GenomeToGenbank import GenomeToGenbank
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


# Used to store and pass around configuration URLs more easily
class SDKConfig:
    def __init__(self, config, version=None):
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.sharedFolder = config['scratch']
        self.srvWizURL = config['srv-wiz-url']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.authServiceUrl = config['auth-service-url']
        self.re_api_url = config['re-api-url']
        self.raw = config
        self.version = version


def _validate_catalog_param_type(threads_count, var_name, default_val, expected_type):
    """Helper function to validate catalog params type"""
    if threads_count is None:
        print(f"Cannot retrieve {var_name} from the catalog, set {var_name}={default_val}")
        return default_val
    print(f"Successfully retrieved {var_name} from the catalog!")
    try:
        threads_count = expected_type(threads_count)
    except ValueError as e:
        raise ValueError(f"{var_name} must be of type {expected_type.__name__}") from e
    return threads_count
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
    VERSION = "0.11.6"
    GIT_URL = "git@github.com:kbaseapps/GenomeFileUtil.git"
    GIT_COMMIT_HASH = "adfb3a4a5e9fa66309fc38f6f585ca2e8c1055a4"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.cfg = SDKConfig(config, self.VERSION)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        max_threads = os.environ.get("KBASE_SECURE_CONFIG_PARAM_MAX_THREADS")
        threads_per_cpu = os.environ.get("KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU")
        self.max_threads = _validate_catalog_param_type(
            max_threads, "MAX_THREADS", MAX_THREADS_DEFAULT, int
        )
        self.threads_per_cpu = _validate_catalog_param_type(
            threads_per_cpu, "THREADS_PER_CPU", THREADS_PER_CPU_DEFAULT, float
        )
        #END_CONSTRUCTOR
        pass


    def genbank_to_genome(self, ctx, params):
        """
        :param params: instance of type "GenbankToGenomeParams" (genome_name
           - becomes the name of the object workspace_id - the immutable,
           numeric ID of the target workspace. Always prefer providing the ID
           over the name. workspace_name - the name of the workspace it gets
           saved to. source - Source of the file typically something like
           RefSeq or Ensembl taxon_ws_name - where the reference taxons are :
           ReferenceTaxons taxon_id - if defined, will try to link the Genome
           to the specified taxonomy id in lieu of performing the lookup
           during upload release - Release or version number of the data per
           example Ensembl has numbered releases of all their data: Release
           31 generate_ids_if_needed - If field used for feature id is not
           there, generate ids (default behavior is raising an exception)
           genetic_code - Genetic code of organism. Overwrites determined GC
           from taxon object scientific_name - will be used to set the
           scientific name of the genome and link to a taxon
           generate_missing_genes - If the file has CDS or mRNA with no
           corresponding gene, generate a spoofed gene. use_existing_assembly
           - Supply an existing assembly reference) -> structure: parameter
           "file" of type "File" -> structure: parameter "path" of String,
           parameter "shock_id" of String, parameter "ftp_url" of String,
           parameter "workspace_id" of Long, parameter "genome_name" of
           String, parameter "workspace_name" of String, parameter "source"
           of String, parameter "taxon_wsname" of String, parameter
           "taxon_id" of String, parameter "release" of String, parameter
           "generate_ids_if_needed" of String, parameter "genetic_code" of
           Long, parameter "scientific_name" of String, parameter "metadata"
           of type "usermeta" -> mapping from String to String, parameter
           "generate_missing_genes" of type "boolean" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter
           "use_existing_assembly" of String
        :returns: instance of type "GenbankToGenomeSaveResult" -> structure:
           parameter "genome_ref" of String, parameter "assembly_ref" of
           String, parameter "assembly_path" of String, parameter
           "assembly_info" of type "object_info" (Information about an
           object, including user provided metadata. objid - the numerical id
           of the object. name - the name of the object. type - the type of
           the object. save_date - the save date of the object. ver - the
           version of the object. saved_by - the user that saved or copied
           the object. wsid - the id of the workspace containing the object.
           workspace - the name of the workspace containing the object. chsum
           - the md5 checksum of the object. size - the size of the object in
           bytes. meta - arbitrary user-supplied metadata about the object.)
           -> tuple of size 11: parameter "objid" of Long, parameter "name"
           of String, parameter "type" of String, parameter "save_date" of
           String, parameter "version" of Long, parameter "saved_by" of
           String, parameter "wsid" of Long, parameter "workspace" of String,
           parameter "chsum" of String, parameter "size" of Long, parameter
           "meta" of mapping from String to String, parameter "genome_info"
           of type "object_info" (Information about an object, including user
           provided metadata. obj_id objid - the numerical id of the object.
           obj_name name - the name of the object. type_string type - the
           type of the object. timestamp save_date - the save date of the
           object. obj_ver ver - the version of the object. username saved_by
           - the user that saved or copied the object. ws_id wsid - the
           workspace containing the object. ws_name workspace - the workspace
           containing the object. string chsum - the md5 checksum of the
           object. int size - the size of the object in bytes. usermeta meta
           - arbitrary user-supplied metadata about the object.) -> tuple of
           size 11: parameter "objid" of type "obj_id" (The unique, permanent
           numerical ID of an object.), parameter "name" of type "obj_name"
           (A string used as a name for an object. Any string consisting of
           alphanumeric characters and the characters |._- that is not an
           integer is acceptable.), parameter "type" of type "type_string" (A
           type string. Specifies the type and its version in a single string
           in the format [module].[typename]-[major].[minor]: module - a
           string. The module name of the typespec containing the type.
           typename - a string. The name of the type as assigned by the
           typedef statement. major - an integer. The major version of the
           type. A change in the major version implies the type has changed
           in a non-backwards compatible way. minor - an integer. The minor
           version of the type. A change in the minor version implies that
           the type has changed in a way that is backwards compatible with
           previous type definitions. In many cases, the major and minor
           versions are optional, and if not provided the most recent version
           will be used. Example: MyModule.MyType-3.1), parameter "save_date"
           of type "timestamp" (A time in the format YYYY-MM-DDThh:mm:ssZ,
           where Z is either the character Z (representing the UTC timezone)
           or the difference in time to UTC in the format +/-HHMM, eg:
           2012-12-17T23:24:06-0500 (EST time) 2013-04-03T08:56:32+0000 (UTC
           time) 2013-04-03T08:56:32Z (UTC time)), parameter "version" of
           Long, parameter "saved_by" of type "username" (Login name of a
           KBase user account.), parameter "wsid" of type "ws_id" (The
           unique, permanent numerical ID of a workspace.), parameter
           "workspace" of type "ws_name" (A string used as a name for a
           workspace. Any string consisting of alphanumeric characters and
           "_", ".", or "-" that is not an integer is acceptable. The name
           may optionally be prefixed with the workspace owner's user name
           and a colon, e.g. kbasetest:my_workspace.), parameter "chsum" of
           String, parameter "size" of Long, parameter "meta" of type
           "usermeta" (User provided metadata about an object. Arbitrary
           key-value pairs provided by the user.) -> mapping from String to
           String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genbank_to_genome
        print('genbank_to_genome -- parameters = ')
        pprint(params)

        importer = GenbankToGenome(self.cfg)
        result = importer.import_genbank(params)

        print('import complete -- result = ')
        pprint(result)
        #END genbank_to_genome

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genbank_to_genome return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def genbanks_to_genomes(self, ctx, params):
        """
        :param params: instance of type "GenbanksToGenomesParams" ->
           structure: parameter "workspace_id" of Long, parameter "inputs" of
           list of type "GenbankToGenomeInput" -> structure: parameter "file"
           of type "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "ftp_url" of String, parameter
           "genome_name" of String, parameter "source" of String, parameter
           "taxon_wsname" of String, parameter "taxon_id" of String,
           parameter "release" of String, parameter "generate_ids_if_needed"
           of String, parameter "genetic_code" of Long, parameter
           "scientific_name" of String, parameter "metadata" of type
           "usermeta" -> mapping from String to String, parameter
           "generate_missing_genes" of type "boolean" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter
           "use_existing_assembly" of String
        :returns: instance of type "GenbanksToGenomesSaveResults" (Results
           for the genbanks_to_genomes function. results - the results of the
           save operation in the same order as the input.) -> structure:
           parameter "results" of list of type "GenbankToGenomeSaveResult" ->
           structure: parameter "genome_ref" of String, parameter
           "assembly_ref" of String, parameter "assembly_path" of String,
           parameter "assembly_info" of type "object_info" (Information about
           an object, including user provided metadata. objid - the numerical
           id of the object. name - the name of the object. type - the type
           of the object. save_date - the save date of the object. ver - the
           version of the object. saved_by - the user that saved or copied
           the object. wsid - the id of the workspace containing the object.
           workspace - the name of the workspace containing the object. chsum
           - the md5 checksum of the object. size - the size of the object in
           bytes. meta - arbitrary user-supplied metadata about the object.)
           -> tuple of size 11: parameter "objid" of Long, parameter "name"
           of String, parameter "type" of String, parameter "save_date" of
           String, parameter "version" of Long, parameter "saved_by" of
           String, parameter "wsid" of Long, parameter "workspace" of String,
           parameter "chsum" of String, parameter "size" of Long, parameter
           "meta" of mapping from String to String, parameter "genome_info"
           of type "object_info" (Information about an object, including user
           provided metadata. obj_id objid - the numerical id of the object.
           obj_name name - the name of the object. type_string type - the
           type of the object. timestamp save_date - the save date of the
           object. obj_ver ver - the version of the object. username saved_by
           - the user that saved or copied the object. ws_id wsid - the
           workspace containing the object. ws_name workspace - the workspace
           containing the object. string chsum - the md5 checksum of the
           object. int size - the size of the object in bytes. usermeta meta
           - arbitrary user-supplied metadata about the object.) -> tuple of
           size 11: parameter "objid" of type "obj_id" (The unique, permanent
           numerical ID of an object.), parameter "name" of type "obj_name"
           (A string used as a name for an object. Any string consisting of
           alphanumeric characters and the characters |._- that is not an
           integer is acceptable.), parameter "type" of type "type_string" (A
           type string. Specifies the type and its version in a single string
           in the format [module].[typename]-[major].[minor]: module - a
           string. The module name of the typespec containing the type.
           typename - a string. The name of the type as assigned by the
           typedef statement. major - an integer. The major version of the
           type. A change in the major version implies the type has changed
           in a non-backwards compatible way. minor - an integer. The minor
           version of the type. A change in the minor version implies that
           the type has changed in a way that is backwards compatible with
           previous type definitions. In many cases, the major and minor
           versions are optional, and if not provided the most recent version
           will be used. Example: MyModule.MyType-3.1), parameter "save_date"
           of type "timestamp" (A time in the format YYYY-MM-DDThh:mm:ssZ,
           where Z is either the character Z (representing the UTC timezone)
           or the difference in time to UTC in the format +/-HHMM, eg:
           2012-12-17T23:24:06-0500 (EST time) 2013-04-03T08:56:32+0000 (UTC
           time) 2013-04-03T08:56:32Z (UTC time)), parameter "version" of
           Long, parameter "saved_by" of type "username" (Login name of a
           KBase user account.), parameter "wsid" of type "ws_id" (The
           unique, permanent numerical ID of a workspace.), parameter
           "workspace" of type "ws_name" (A string used as a name for a
           workspace. Any string consisting of alphanumeric characters and
           "_", ".", or "-" that is not an integer is acceptable. The name
           may optionally be prefixed with the workspace owner's user name
           and a colon, e.g. kbasetest:my_workspace.), parameter "chsum" of
           String, parameter "size" of Long, parameter "meta" of type
           "usermeta" (User provided metadata about an object. Arbitrary
           key-value pairs provided by the user.) -> mapping from String to
           String
        """
        # ctx is the context object
        # return variables are: results
        #BEGIN genbanks_to_genomes
        results = {
            'results': GenbankToGenome(self.cfg).import_genbank_mass(params)
        }
        #END genbanks_to_genomes

        # At some point might do deeper type checking...
        if not isinstance(results, dict):
            raise ValueError('Method genbanks_to_genomes return value ' +
                             'results is not type dict as required.')
        # return the results
        return [results]

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
           generated during this call.) -> structure: parameter "file_path"
           of String, parameter "from_cache" of type "boolean" (A boolean - 0
           for false, 1 for true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_to_gff
        print('genome_to_gff -- parameters = ')
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

    def metagenome_to_gff(self, ctx, params):
        """
        :param params: instance of type "MetagenomeToGFFParams" (is_gtf -
           optional flag switching export to GTF format (default is 0, which
           means GFF) target_dir - optional target directory to create file
           in (default is temporary folder with name 'gff_<timestamp>'
           created in scratch)) -> structure: parameter "metagenome_ref" of
           String, parameter "ref_path_to_genome" of list of String,
           parameter "is_gtf" of type "boolean" (A boolean - 0 for false, 1
           for true. @range (0, 1)), parameter "target_dir" of String
        :returns: instance of type "MetagenomeToGFFResult" -> structure:
           parameter "file_path" of String, parameter "from_cache" of type
           "boolean" (A boolean - 0 for false, 1 for true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN metagenome_to_gff
        if not params.get('metagenome_ref'):
            raise ValueError(f"argument 'metagenome_ref' required for function metagenome_to_gff")
        
        # quick hack to make function work with other code.
        params['genome_ref'] = params['metagenome_ref']

        exporter = GenomeToGFF(self.cfg)
        result = exporter.export(ctx, params)

        print('export complete -- result = ')
        pprint(result)
        #END metagenome_to_gff

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method metagenome_to_gff return value ' +
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
           "genbank_file" of type "GBFile" -> structure: parameter
           "file_path" of String, parameter "from_cache" of type "boolean" (A
           boolean - 0 for false, 1 for true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_to_genbank
        print('genome_to_genbank -- parameters = ')
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

    def genome_features_to_fasta(self, ctx, params):
        """
        :param params: instance of type "GenomeFeaturesToFastaParams"
           (Produce a FASTA file with the nucleotide sequences of features in
           a genome. string genome_ref: reference to a genome object
           list<string> feature_lists: Optional, which features lists
           (features, mrnas, cdss, non_coding_features) to provide sequences.
           Defaults to "features". list<string> filter_ids: Optional, if
           provided only return sequences for matching features. boolean
           include_functions: Optional, add function to header line. Defaults
           to True. boolean include_aliases: Optional, add aliases to header
           line. Defaults to True.) -> structure: parameter "genome_ref" of
           String, parameter "feature_lists" of list of String, parameter
           "filter_ids" of list of String, parameter "include_functions" of
           type "boolean" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "include_aliases" of type "boolean" (A boolean - 0
           for false, 1 for true. @range (0, 1))
        :returns: instance of type "FASTAResult" -> structure: parameter
           "file_path" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_features_to_fasta
        logging.info(f"Running genome_features_to_fasta with the following params: {params}")
        exporter = GenomeFeaturesToFasta(self.cfg)
        result = exporter.export(ctx, params, protein=False)
        #END genome_features_to_fasta

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genome_features_to_fasta return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def genome_proteins_to_fasta(self, ctx, params):
        """
        :param params: instance of type "GenomeProteinToFastaParams" (Produce
           a FASTA file with the protein sequences of CDSs in a genome.
           string genome_ref: reference to a genome object list<string>
           filter_ids: Optional, if provided only return sequences for
           matching features. boolean include_functions: Optional, add
           function to header line. Defaults to True. boolean
           include_aliases: Optional, add aliases to header line. Defaults to
           True.) -> structure: parameter "genome_ref" of String, parameter
           "filter_ids" of list of String, parameter "include_functions" of
           type "boolean" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "include_aliases" of type "boolean" (A boolean - 0
           for false, 1 for true. @range (0, 1))
        :returns: instance of type "FASTAResult" -> structure: parameter
           "file_path" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_proteins_to_fasta
        logging.info(f"Running genome_proteins_to_fasta with the following params: {params}")
        exporter = GenomeFeaturesToFasta(self.cfg)
        result = exporter.export(ctx, params, protein=True)
        #END genome_proteins_to_fasta

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genome_proteins_to_fasta return value ' +
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
        print('export_genome_as_genbank -- parameters = ')

        # validate parameters
        if 'input_ref' not in params:
            raise ValueError('Cannot run export_genome_as_genbank- no "input_ref" field defined.')

        # get WS metadata to get ws_name and obj_name
        ws = Workspace(url=self.cfg.workspaceURL)
        info = ws.get_object_info_new({'objects': [{'ref': params['input_ref']}],
                                       'includeMetadata': 0, 'ignoreErrors': 0})[0]

        genome_to_genbank_params = {
            'genome_ref': params['input_ref']
        }

        # export to file (building from KBase Genome Object)
        result = self.genome_to_genbank(ctx, genome_to_genbank_params)[0]['genbank_file']

        # create the output directory and move the file there
        export_package_dir = os.path.join(self.cfg.sharedFolder, info[1])
        os.makedirs(export_package_dir)
        shutil.move(
          result['file_path'],
          os.path.join(export_package_dir, os.path.basename(result['file_path'])))

        # export original uploaded GenBank file if it existed.
        exporter = GenomeToGenbank(self.cfg)
        original_result_full = exporter.export_original_genbank(ctx, genome_to_genbank_params)
        if original_result_full is not None:
            original_result = original_result_full['genbank_file']
            shutil.move(
              original_result['file_path'],
              os.path.join(export_package_dir, os.path.basename(original_result['file_path'])))

        # Make warning file about genes only.
        warning_filename = "README.txt"
        with open(os.path.join(export_package_dir, warning_filename), 'w') as temp_file:
            temp_file.write('This directory includes the KBase-derived GenBank file and also ' +
                            '(if you originally uploaded the genome from an annotated ' +
                            'GenBank file) the original GenBank input.')

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

    def export_genome_as_gff(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_genome_as_gff
        if 'input_ref' not in params:
            raise ValueError('Cannot run export_genome_as_gff- no "input_ref" '
                             'field defined.')

        # get WS metadata to get ws_name and obj_name
        ws = Workspace(url=self.cfg.workspaceURL)
        info = ws.get_objects2({'objects': [{
            'ref': params['input_ref'],
            'included':['/assembly_ref', '/contigset_ref', '/id', '/gff_handle_ref']}
        ]})['data'][0]['data']

        # export to file (building from KBase Genome Object)
        result = self.genome_to_gff(ctx, {
            'genome_ref': params['input_ref']
        })[0]

        # get assembly
        if 'assembly_ref' in info:
            assembly_ref = info['assembly_ref']
        else:
            assembly_ref = info['contigset_ref']
        print(('Assembly reference = ' + assembly_ref))
        print('Downloading assembly')
        au = AssemblyUtil(self.cfg.callbackURL)
        assembly_file_path = au.get_assembly_as_fasta(
            {'ref': params['input_ref'] + ";" + assembly_ref}
        )['path']

        # create the output directory and move the files there
        export_package_dir = os.path.join(self.cfg.sharedFolder, info['id'])
        os.makedirs(export_package_dir)
        shutil.move(
            result['file_path'],
            os.path.join(export_package_dir,
                         'KBase_derived_' + os.path.basename(result['file_path'])))
        shutil.move(
            assembly_file_path,
            os.path.join(export_package_dir,
                         os.path.basename(assembly_file_path)))

        # add cached genome if appropriate
        exporter = GenomeToGFF(self.cfg)
        cached = exporter.get_gff_handle(info, export_package_dir)

        # package it up
        dfUtil = DataFileUtil(self.cfg.callbackURL)
        package_details = dfUtil.package_for_download({
                                    'file_path': export_package_dir,
                                    'ws_refs': [params['input_ref']]
                                })

        output = {'shock_id': package_details['shock_id']}
        #END export_genome_as_gff

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_genome_as_gff return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def export_genome_features_protein_to_fasta(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_genome_features_protein_to_fasta
        print('export_genome_features_protein_to_fasta -- parameters = ')

        # validate parameters
        if 'input_ref' not in params:
            raise ValueError('Cannot run export_genome_features_protein_to_fasta - no "input_ref" field defined.')

        # get WS metadata to get ws_name and obj_name
        ws = Workspace(url=self.cfg.workspaceURL)
        info = ws.get_object_info_new({'objects': [{'ref': params['input_ref']}],
                                       'includeMetadata': 0, 'ignoreErrors': 0})[0]

        genome_to_protein_fasta_params = {
            'genome_ref': params['input_ref']
        }

        # export to file (building from KBase Genome Object)
        result = self.genome_to_genbank(ctx, genome_to_protein_fasta_params)[0]['genbank_file']



        #END export_genome_features_protein_to_fasta

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_genome_features_protein_to_fasta return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def export_metagenome_as_gff(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (input and output
           structure functions for standard downloaders) -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_metagenome_as_gff
        if 'input_ref' not in params:
            raise ValueError('Cannot run export_metagenome_as_gff- no "input_ref" '
                             'field defined.')

        # get WS metadata to get ws_name and obj_name
        ws = Workspace(url=self.cfg.workspaceURL)
        ws_obj = ws.get_objects2({'objects': [{
            'ref': params['input_ref'],
            'included':['/assembly_ref', '/gff_handle_ref']}
        ]})['data'][0]

        info = ws_obj['data']
        ws_info = ws_obj['info']

        # export to file (building from KBase Genome Object)
        result = self.metagenome_to_gff(ctx, {
            'metagenome_ref': params['input_ref']
        })[0]

        # get assembly
        if 'assembly_ref' in info:
            assembly_ref = info['assembly_ref']
        else:
            raise ValueError("No Assembly associated with this AnnotatedMetagenomeAssembly "
                             "object. Cannot retrieve fasta file. ")
        au = AssemblyUtil(self.cfg.callbackURL)
        assembly_file_path = au.get_assembly_as_fasta(
            {'ref': params['input_ref'] + ";" + assembly_ref}
        )['path']

        # create the output directory and move the files there
        export_package_dir = os.path.join(self.cfg.sharedFolder, '_'.join(ws_info[1].split()))
        os.makedirs(export_package_dir)
        shutil.move(
            result['file_path'],
            os.path.join(export_package_dir,
                         'KBase_derived_' + os.path.basename(result['file_path'])))
        shutil.move(
            assembly_file_path,
            os.path.join(export_package_dir,
                         os.path.basename(assembly_file_path)))

        # package it up
        dfUtil = DataFileUtil(self.cfg.callbackURL)
        package_details = dfUtil.package_for_download({
                                    'file_path': export_package_dir,
                                    'ws_refs': [params['input_ref']]
                                })

        output = {'shock_id': package_details['shock_id']}
        #END export_metagenome_as_gff

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_metagenome_as_gff return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def fasta_gff_to_genome(self, ctx, params):
        """
        :param params: instance of type "FastaGFFToGenomeParams" (genome_name
           - becomes the name of the object workspace_name - the name of the
           workspace it gets saved to. source - Source of the file typically
           something like RefSeq or Ensembl taxon_ws_name - where the
           reference taxons are : ReferenceTaxons taxon_id - if defined, will
           try to link the Genome to the specified taxonomy id in lieu of
           performing the lookup during upload release - Release or version
           number of the data per example Ensembl has numbered releases of
           all their data: Release 31 genetic_code - Genetic code of
           organism. Overwrites determined GC from taxon object
           scientific_name - will be used to set the scientific name of the
           genome and link to a taxon generate_missing_genes - If the file
           has CDS or mRNA with no corresponding gene, generate a spoofed
           gene. Off by default existing_assembly_ref - a KBase assembly upa,
           to associate the genome with. Avoids saving a new assembly when
           specified.) -> structure: parameter "fasta_file" of type "File" ->
           structure: parameter "path" of String, parameter "shock_id" of
           String, parameter "ftp_url" of String, parameter "gff_file" of
           type "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "ftp_url" of String, parameter
           "genome_name" of String, parameter "workspace_name" of String,
           parameter "source" of String, parameter "taxon_wsname" of String,
           parameter "taxon_id" of String, parameter "release" of String,
           parameter "genetic_code" of Long, parameter "scientific_name" of
           String, parameter "metadata" of type "usermeta" -> mapping from
           String to String, parameter "generate_missing_genes" of type
           "boolean" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "existing_assembly_ref" of String
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "genome_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN fasta_gff_to_genome
        print('--->\nRunning GenomeFileUtil.fasta_gff_to_genome\nparams:')
        print(json.dumps(params, indent=1))

        for key in list(params.keys()):
            if params[key] is None:
                del params[key]

        for key, value in params.items():
            if isinstance(value, str):
                params[key] = value.strip()

        importer = FastaGFFToGenome(self.cfg)
        returnVal = importer.import_file(params)
        #END fasta_gff_to_genome

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method fasta_gff_to_genome return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def fasta_gff_to_genome_json(self, ctx, params):
        """
        As above but returns the genome instead
        :param params: instance of type "FastaGFFToGenomeParams" (genome_name
           - becomes the name of the object workspace_name - the name of the
           workspace it gets saved to. source - Source of the file typically
           something like RefSeq or Ensembl taxon_ws_name - where the
           reference taxons are : ReferenceTaxons taxon_id - if defined, will
           try to link the Genome to the specified taxonomy id in lieu of
           performing the lookup during upload release - Release or version
           number of the data per example Ensembl has numbered releases of
           all their data: Release 31 genetic_code - Genetic code of
           organism. Overwrites determined GC from taxon object
           scientific_name - will be used to set the scientific name of the
           genome and link to a taxon generate_missing_genes - If the file
           has CDS or mRNA with no corresponding gene, generate a spoofed
           gene. Off by default existing_assembly_ref - a KBase assembly upa,
           to associate the genome with. Avoids saving a new assembly when
           specified.) -> structure: parameter "fasta_file" of type "File" ->
           structure: parameter "path" of String, parameter "shock_id" of
           String, parameter "ftp_url" of String, parameter "gff_file" of
           type "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "ftp_url" of String, parameter
           "genome_name" of String, parameter "workspace_name" of String,
           parameter "source" of String, parameter "taxon_wsname" of String,
           parameter "taxon_id" of String, parameter "release" of String,
           parameter "genetic_code" of Long, parameter "scientific_name" of
           String, parameter "metadata" of type "usermeta" -> mapping from
           String to String, parameter "generate_missing_genes" of type
           "boolean" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "existing_assembly_ref" of String
        :returns: instance of unspecified object
        """
        # ctx is the context object
        # return variables are: genome
        #BEGIN fasta_gff_to_genome_json
        importer = FastaGFFToGenome(self.cfg)
        genome = importer.generate_genome_json(params)
        #END fasta_gff_to_genome_json

        # At some point might do deeper type checking...
        if not isinstance(genome, object):
            raise ValueError('Method fasta_gff_to_genome_json return value ' +
                             'genome is not type object as required.')
        # return the results
        return [genome]

    def fasta_gff_to_metagenome(self, ctx, params):
        """
        :param params: instance of type "FastaGFFToMetagenomeParams"
           (genome_name - becomes the name of the object workspace_name - the
           name of the workspace it gets saved to. source - Source of the
           file typically something like RefSeq or Ensembl taxon_ws_name -
           where the reference taxons are : ReferenceTaxons taxon_id - if
           defined, will try to link the Genome to the specified taxonomy id
           in lieu of performing the lookup during upload release - Release
           or version number of the data per example Ensembl has numbered
           releases of all their data: Release 31 genetic_code - Genetic code
           of organism. Overwrites determined GC from taxon object
           scientific_name - will be used to set the scientific name of the
           genome and link to a taxon generate_missing_genes - If the file
           has CDS or mRNA with no corresponding gene, generate a spoofed
           gene. Off by default existing_assembly_ref - a KBase assembly upa,
           to associate the metagenome with. Avoids saving a new assembly
           when specified.) -> structure: parameter "fasta_file" of type
           "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "ftp_url" of String, parameter
           "gff_file" of type "File" -> structure: parameter "path" of
           String, parameter "shock_id" of String, parameter "ftp_url" of
           String, parameter "genome_name" of String, parameter
           "workspace_name" of String, parameter "source" of String,
           parameter "metadata" of type "usermeta" -> mapping from String to
           String, parameter "generate_missing_genes" of type "boolean" (A
           boolean - 0 for false, 1 for true. @range (0, 1)), parameter
           "existing_assembly_ref" of String
        :returns: instance of type "MetagenomeSaveResult" -> structure:
           parameter "metagenome_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN fasta_gff_to_metagenome

        for key in list(params.keys()):
            if params[key] is None:
                del params[key]

        for key, value in params.items():
            if isinstance(value, str):
                params[key] = value.strip()

        params['is_metagenome'] = True

        importer = FastaGFFToGenome(self.cfg)
        returnVal = importer.import_file(params)
        #END fasta_gff_to_metagenome

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method fasta_gff_to_metagenome return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def save_one_genome(self, ctx, params):
        """
        :param params: instance of type "SaveOneGenomeParams" -> structure:
           parameter "workspace" of String, parameter "name" of String,
           parameter "data" of type "Genome" (Genome type -- annotated and
           assembled genome data. Field descriptions: id - string - KBase
           legacy data ID scientific_name - string - human readable species
           name domain - string - human readable phylogenetic domain name
           (eg. "Bacteria") warnings - list of string - genome-level warnings
           generated in the annotation process genome_tiers - list of string
           - controlled vocabulary (based on app input and checked by
           GenomeFileUtil) A list of labels describing the data source for
           this genome. Allowed values - Representative, Reference,
           ExternalDB, User Tier assignments based on genome source: * All
           phytozome - Representative and ExternalDB * Phytozome flagship
           genomes - Reference, Representative and ExternalDB * Ensembl -
           Representative and ExternalDB * RefSeq Reference - Reference,
           Representative and ExternalDB * RefSeq Representative -
           Representative and ExternalDB * RefSeq Latest or All Assemblies
           folder - ExternalDB * User Data - User tagged feature_counts - map
           of string to integer - total counts of each type of feature keys
           are a controlled vocabulary of - "CDS", "gene", "misc_feature",
           "misc_recomb", "mobile_element", "ncRNA" - 72,
           "non_coding_features", "non_coding_genes",
           "protein_encoding_gene", "rRNA", "rep_origin", "repeat_region",
           "tRNA" genetic_code - int - An NCBI-assigned taxonomic category
           for the organism See here -
           https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi dna_size
           - integer - total number of nucleotides num_contigs - integer -
           total number of contigs in the genome molecule_type - string -
           controlled vocab - the type of molecule sequenced Possible values
           are "Unknown", "DNA", "RNA", "genomic DNA", "genomic RNA", "mRNA",
           "tRNA", "rRNA", "other RNA", "other DNA", "transcribed RNA",
           "viral cRNA", "unassigned DNA", "unassigned RNA" contig_lengths -
           list of int - nucleotide length of each contig in the genome
           Indexes in this list correspond to indexes in the `contig_ids`
           list. contig_ids - list of str - external database identifiers for
           each contig (eg. "NC_000913.3") source - str - controlled vocab -
           descriptor of where this data came from (eg. "RefSeq") Allowed
           entries RefSeq, Ensembl, Phytozome, RAST, Prokka, User_upload
           source_id - string - identifier of this genome from the source
           database (eg. the RefSeq ID such as "NC_000913") md5 - string -
           checksum of the underlying assembly sequence taxonomy - string -
           semicolon-delimited taxonomy lineage, in order of parent to child
           taxon_assignments - mapping of taxonomy namespace to taxon ID.
           example - {"ncbi": "286", "gtdb": "s__staphylococcus_devriesei"}
           gc_content - float - ratio of GC count to AT in the genome
           publications - tuple of (pubmedid, source, title, web_addr, year,
           authors, journal). See typedef above. ontology_events - A record
           of the service and method used for a set of ontology assignments
           on the genome. ontologies_present - a mapping of ontology source
           id (eg. "GO") to a mapping of term IDs (eg "GO:16209") to term
           names (eg. "histidine biosynthetic process"). features - array of
           Feature - protein coding genes (see the separate Feature spec)
           cdss - array of protein-coding sequences mrnas - array of
           transcribed messenger RNA sequences (equal to cdss plus 5' and 3'
           UTRs) non_coding_features - array of features that does not
           include mRNA, CDS, and protein-encoding genes assembly_ref -
           workspace reference to an assembly object from which this
           annotated genome was derived. taxon_ref - workspace reference to a
           taxon object that classifies the species or strain of this genome.
           genbank_handle_ref - file server handle reference to the source
           genbank file for this genome. gff_handle_ref - file server handle
           reference to the source GFF file for this genome.
           external_source_origination_date - TODO look at GFU for this
           release - string - User-supplied release or version of the source
           data. This most likely will come from an input field in the import
           app. original_source_file_name - filename from which this genome
           was derived (eg. genbank or gff filename). notes - TODO
           quality_scores - TODO suspect - bool - flag of whether this
           annotation is problematic due to some warning genome_type - string
           - controlled vocab - One of "draft isolate", "finished isolate",
           "mag", "sag", "virus", "plasmid", "construct" Features vs. coding
           sequences: a feature is a sequence in the DNA that codes for a
           protein, including non-transcribed introns. A coding sequence
           (stored as `cdss`) includes **only** the sections of the feature
           that codes for a protein, minus introns and UTRs. @optional
           warnings contig_lengths contig_ids source_id taxonomy publications
           @optional ontology_events ontologies_present non_coding_features
           mrnas genome_type @optional genbank_handle_ref gff_handle_ref
           external_source_origination_date @optional release
           original_source_file_name notes quality_scores suspect
           assembly_ref @optional taxon_ref taxon_assignments @metadata ws
           gc_content as GC content @metadata ws taxonomy as Taxonomy
           @metadata ws md5 as MD5 @metadata ws dna_size as Size @metadata ws
           genetic_code as Genetic code @metadata ws domain as Domain
           @metadata ws source_id as Source ID @metadata ws source as Source
           @metadata ws scientific_name as Name @metadata ws genome_type as
           Type @metadata ws length(features) as Number of Protein Encoding
           Genes @metadata ws length(cdss) as Number of CDS @metadata ws
           assembly_ref as Assembly Object @metadata ws num_contigs as Number
           contigs @metadata ws length(warnings) as Number of Genome Level
           Warnings @metadata ws suspect as Suspect Genome) -> structure:
           parameter "id" of type "Genome_id" (KBase legacy data ID @id kb),
           parameter "scientific_name" of String, parameter "domain" of
           String, parameter "warnings" of list of String, parameter
           "genome_tiers" of list of String, parameter "feature_counts" of
           mapping from String to Long, parameter "genetic_code" of Long,
           parameter "dna_size" of Long, parameter "num_contigs" of Long,
           parameter "molecule_type" of String, parameter "contig_lengths" of
           list of Long, parameter "contig_ids" of list of String, parameter
           "source" of String, parameter "source_id" of type "source_id"
           (Reference to a source_id @id external), parameter "md5" of
           String, parameter "taxonomy" of String, parameter
           "taxon_assignments" of mapping from String to String, parameter
           "gc_content" of Double, parameter "publications" of list of type
           "publication" (Structure for a publication Elements: (0) pubmedid
           - float (1) source - string - (ex. Pubmed) (2) title - string (3)
           string web address - string (4) publication year - string (5)
           authors - string (6) journal - string) -> tuple of size 7:
           parameter "pubmedid" of Double, parameter "source" of String,
           parameter "title" of String, parameter "url" of String, parameter
           "year" of String, parameter "authors" of String, parameter
           "journal" of String, parameter "ontology_events" of list of type
           "Ontology_event" (@optional ontology_ref method_version eco) ->
           structure: parameter "id" of String, parameter "ontology_ref" of
           type "Ontology_ref" (Reference to a ontology object @id ws
           KBaseOntology.OntologyDictionary), parameter "method" of String,
           parameter "method_version" of String, parameter "timestamp" of
           String, parameter "eco" of String, parameter "ontologies_present"
           of mapping from String to mapping from String to String, parameter
           "features" of list of type "Feature" (Structure for a single CDS
           encoding “gene” of a genome ONLY PUT GENES THAT HAVE A
           CORRESPONDING CDS IN THIS ARRAY NOTE: Sequence is optional.
           Ideally we can keep it in here, but Recognize due to space
           constraints another solution may be needed. We may want to add
           additional fields for other CDM functions (e.g., atomic regulons,
           coexpressed fids, co_occurring fids,...)
           protein_translation_length and protein_translation are for longest
           coded protein (representative protein for splice variants) NOTE:
           New Aliases field definitely breaks compatibility. As Does
           Function. flags are flag fields in GenBank format. This will be a
           controlled vocabulary. Initially Acceptable values are pseudo,
           ribosomal_slippage, and trans_splicing Md5 is the md5 of
           dna_sequence. @optional functions ontology_terms note
           protein_translation mrnas flags warnings @optional inference_data
           dna_sequence aliases db_xrefs children functional_descriptions) ->
           structure: parameter "id" of type "Feature_id" (KBase Feature ID
           @id external), parameter "location" of list of tuple of size 4:
           type "Contig_id" (ContigSet contig ID @id external), Long, String,
           Long, parameter "functions" of list of String, parameter
           "functional_descriptions" of list of String, parameter
           "ontology_terms" of mapping from String to mapping from String to
           list of Long, parameter "note" of String, parameter "md5" of
           String, parameter "protein_translation" of String, parameter
           "protein_translation_length" of Long, parameter "cdss" of list of
           String, parameter "mrnas" of list of String, parameter "children"
           of list of String, parameter "flags" of list of String, parameter
           "warnings" of list of String, parameter "inference_data" of list
           of type "InferenceInfo" (Type spec for the "InferenceInfo" object.
           TODO docs Found in the `inference_data` fields in mRNAs and CDSs
           Fields: category - string - TODO type - string - TODO evidence -
           string - TODO) -> structure: parameter "category" of String,
           parameter "type" of String, parameter "evidence" of String,
           parameter "dna_sequence" of String, parameter
           "dna_sequence_length" of Long, parameter "aliases" of list of
           tuple of size 2: parameter "fieldname" of String, parameter
           "alias" of String, parameter "db_xrefs" of list of tuple of size
           2: parameter "db_source" of String, parameter "db_identifier" of
           String, parameter "non_coding_features" of list of type
           "NonCodingFeature" (Structure for a single feature that is NOT one
           of the following: - Protein encoding gene (gene that has a
           corresponding CDS) - mRNA - CDS Note pseudo-genes and Non protein
           encoding genes are put into this flags are flag fields in GenBank
           format. This will be a controlled vocabulary. Initially Acceptable
           values are pseudo, ribosomal_slippage, and trans_splicing Md5 is
           the md5 of dna_sequence. @optional functions ontology_terms note
           flags warnings functional_descriptions @optional inference_data
           dna_sequence aliases db_xrefs children parent_gene) -> structure:
           parameter "id" of type "Feature_id" (KBase Feature ID @id
           external), parameter "location" of list of tuple of size 4: type
           "Contig_id" (ContigSet contig ID @id external), Long, String,
           Long, parameter "type" of String, parameter "functions" of list of
           String, parameter "functional_descriptions" of list of String,
           parameter "ontology_terms" of mapping from String to mapping from
           String to list of Long, parameter "note" of String, parameter
           "md5" of String, parameter "parent_gene" of String, parameter
           "children" of list of String, parameter "flags" of list of String,
           parameter "warnings" of list of String, parameter "inference_data"
           of list of type "InferenceInfo" (Type spec for the "InferenceInfo"
           object. TODO docs Found in the `inference_data` fields in mRNAs
           and CDSs Fields: category - string - TODO type - string - TODO
           evidence - string - TODO) -> structure: parameter "category" of
           String, parameter "type" of String, parameter "evidence" of
           String, parameter "dna_sequence" of String, parameter
           "dna_sequence_length" of Long, parameter "aliases" of list of
           tuple of size 2: parameter "fieldname" of String, parameter
           "alias" of String, parameter "db_xrefs" of list of tuple of size
           2: parameter "db_source" of String, parameter "db_identifier" of
           String, parameter "cdss" of list of type "CDS" (Structure for a
           single coding sequence. Coding sequences are the sections of a
           feature's sequence that are translated to a protein (minus introns
           and UTRs). Fields: id - string - identifier of the coding
           sequence, such as "b0001_CDS_1" location - list<tuple<string, int,
           string, int>> - list of locations from where this sequence
           originates in the original assembly. Each sub-sequence in the list
           constitutes a section of the resulting CDS. The first element in
           the tuple corresponds to the "contig_id", such as "NC_000913.3".
           The second element in the tuple is an index in the contig of where
           the sequence starts. The third element is either a plus or minus
           sign indicating whether it is on the 5' to 3' leading strand ("+")
           or on the 3' to 5' lagging strand ("-"). The last element is the
           length of the sub-sequence. For a location on the leading strand
           (denoted by "+"), the index is of the leftmost base, and the
           sequence extends to the right. For a location on the lagging
           strand (denoted by "-"), the index is of the rightmost base, and
           the sequence extends to the left. NOTE: the last element in each
           tuple is the *length* of each sub-sequence. If you have a location
           such as ("xyz", 100, "+", 50), then your sequence will go from
           index 100 to index 149 (this has a length of 50). It *does not* go
           from index 100 to index 150, as that would have a length of 51.
           Likewise, if you have the location ("xyz", 100, "-", 50), then the
           sequence extends from 100 down to 51, which has a length of 50
           bases. It does not go from index 100 to 50, as that would have a
           length of 51. md5 - string - md5 of the dna sequence - TODO
           clarification protein_md5 - string - hash of the protein sequence
           that this CDS encodes parent_gene - string - gene (feature) from
           which this CDS comes from, including introns and UTRs that have
           been removed to create this CDS. parent_mrna - string - mRNA
           sequence from which this sequence is derived, including UTRs but
           not introns. note - string - TODO functions - list<string> - list
           of protein products or chemical processes that this sequence
           creates, facilitates, or influences. functional_descriptions -
           list<string> - TODO list of protein products or chemical processes
           that sequence creates, facilitates, or influences. ontology_terms
           - mapping<string, mapping<string, list<int>>> - a mapping of
           ontology source id (eg. "GO") to a mapping of term IDs (eg
           "GO:16209") to a list of indexes into the ontology_events data
           (found in the top level of the genome object). The index into an
           ontology event indicates what service and method created this term
           assignment. flags - list<string>  - (controlled vocab) fields from
           the genbank source. A common example is "pseudo" for pseudo-genes
           that do not encode proteins, which shows up as "/pseudo" in the
           genbank. Values can be: "pseudo", "ribosomal_slippage",
           "trans_splicing" warnings - list<string> - TODO inference_data -
           list<InferenceInfo> - TODO protein_translation - string - amino
           acid sequence that this CDS gets translated into.
           protein_translation_length - int - length of the above aliases -
           list<(string, string)> - alternative list of names or identifiers
           eg: [["gene", "thrA"], ["locus_tag", "b0002"]] db_xrefs -
           list<(string, string)> - Identifiers from other databases
           (database cross-references) The first string is the database name,
           the second is the database identifier. eg: [["ASAP",
           "ABE-0000006"], ["EcoGene", "EG11277"]] dna_sequence - string -
           sequence of exons from the genome that constitute this protein
           encoding sequence. dna_sequence_length - int - length of the above
           @optional parent_gene parent_mrna functions ontology_terms note
           flags warnings @optional inference_data dna_sequence aliases
           db_xrefs functional_descriptions) -> structure: parameter "id" of
           type "cds_id" (KBase CDS ID @id external), parameter "location" of
           list of tuple of size 4: type "Contig_id" (ContigSet contig ID @id
           external), Long, String, Long, parameter "md5" of String,
           parameter "protein_md5" of String, parameter "parent_gene" of type
           "Feature_id" (KBase Feature ID @id external), parameter
           "parent_mrna" of type "mrna_id" (KBase mRNA ID @id external),
           parameter "note" of String, parameter "functions" of list of
           String, parameter "functional_descriptions" of list of String,
           parameter "ontology_terms" of mapping from String to mapping from
           String to list of Long, parameter "flags" of list of String,
           parameter "warnings" of list of String, parameter "inference_data"
           of list of type "InferenceInfo" (Type spec for the "InferenceInfo"
           object. TODO docs Found in the `inference_data` fields in mRNAs
           and CDSs Fields: category - string - TODO type - string - TODO
           evidence - string - TODO) -> structure: parameter "category" of
           String, parameter "type" of String, parameter "evidence" of
           String, parameter "protein_translation" of String, parameter
           "protein_translation_length" of Long, parameter "aliases" of list
           of tuple of size 2: parameter "fieldname" of String, parameter
           "alias" of String, parameter "db_xrefs" of list of tuple of size
           2: parameter "db_source" of String, parameter "db_identifier" of
           String, parameter "dna_sequence" of String, parameter
           "dna_sequence_length" of Long, parameter "mrnas" of list of type
           "mRNA" (The mRNA is the transcribed sequence from the original
           feature, minus the introns, but including the UTRs. Fields: id -
           string - identifying string for the mRNA location -
           list<tuple<string, int, string, int>> - list of locations from
           where this sequence originates in the original assembly. Each
           sub-sequence in the list constitutes a section of the resulting
           CDS. The first element in the tuple corresponds to the
           "contig_id", such as "NC_000913.3". The second element in the
           tuple is an index in the contig of where the sequence starts. The
           third element is either a plus or minus sign indicating whether it
           is on the 5' to 3' leading strand ("+") or on the 3' to 5' lagging
           strand ("-"). The last element is the length of the sub-sequence.
           For a location on the leading strand (denoted by "+"), the index
           is of the leftmost base, and the sequence extends to the right.
           For a location on the lagging strand (denoted by "-"), the index
           is of the rightmost base, and the sequence extends to the left.
           NOTE: the last element in each tuple is the *length* of each
           sub-sequence. If you have a location such as ("xyz", 100, "+",
           50), then your sequence will go from index 100 to index 149 (this
           has a length of 50). It *does not* go from index 100 to index 150,
           as that would have a length of 51. Likewise, if you have the
           location ("xyz", 100, "-", 50), then the sequence extends from 100
           down to 51, which has a length of 50 bases. It does not go from
           index 100 to 50, as that would have a length of 51. md5 - string -
           md5 of the dna sequence - TODO clarification parent_gene -
           Feature_id - corresponding feature for this sequence, including
           introns and UTRs cds - string - corresponding coding sequence for
           this mRNA (the sequence minus UTRs) dna_sequence - string -
           sequence of UTRs and exons from the genome that constitute this
           mRNA dna_sequence_length - int - length of the above note - string
           - TODO functions - list<string> - TODO list of protein products or
           chemical processes that sequence creates, facilitates, or
           influences. functional_descriptions - list<string> - TODO list of
           protein products or chemical processes that sequence creates,
           facilitates, or influences. ontology_terms - mapping<string,
           mapping<string, list<int>>> - a mapping of ontology source id (eg.
           "GO") to a mapping of term IDs (eg "GO:16209") to a list of
           indexes into the ontology_events data (found in the top level of
           the genome object). The index into an ontology event indicates
           what service and method created this term assignment. flags -
           list<string> - controlled vocab - fields from the genbank source.
           A common example is "pseudo" for pseudo-genes that do not encode
           proteins, which shows up as "/pseudo" in the genbank. Values can
           be: "pseudo", "ribosomal_slippage", "trans_splicing" warnings -
           list<string> - TODO inference_data - list<InferenceInfo> - TODO
           aliases - list<(string, string)> - alternative list of names or
           identifiers eg: [["gene", "thrA"], ["locus_tag", "b0002"]]
           db_xrefs - list<(string, string)> - Identifiers from other
           databases (database cross-references). The first string is the
           database name, the second is the database identifier. eg:
           [["ASAP", "ABE-0000006"], ["EcoGene", "EG11277"]] @optional
           parent_gene cds functions ontology_terms note flags warnings
           @optional inference_data dna_sequence aliases db_xrefs
           functional_descriptions) -> structure: parameter "id" of type
           "mrna_id" (KBase mRNA ID @id external), parameter "location" of
           list of tuple of size 4: type "Contig_id" (ContigSet contig ID @id
           external), Long, String, Long, parameter "md5" of String,
           parameter "parent_gene" of type "Feature_id" (KBase Feature ID @id
           external), parameter "cds" of type "cds_id" (KBase CDS ID @id
           external), parameter "dna_sequence" of String, parameter
           "dna_sequence_length" of Long, parameter "note" of String,
           parameter "functions" of list of String, parameter
           "functional_descriptions" of list of String, parameter
           "ontology_terms" of mapping from String to mapping from String to
           list of Long, parameter "flags" of list of String, parameter
           "warnings" of list of String, parameter "inference_data" of list
           of type "InferenceInfo" (Type spec for the "InferenceInfo" object.
           TODO docs Found in the `inference_data` fields in mRNAs and CDSs
           Fields: category - string - TODO type - string - TODO evidence -
           string - TODO) -> structure: parameter "category" of String,
           parameter "type" of String, parameter "evidence" of String,
           parameter "aliases" of list of tuple of size 2: parameter
           "fieldname" of String, parameter "alias" of String, parameter
           "db_xrefs" of list of tuple of size 2: parameter "db_source" of
           String, parameter "db_identifier" of String, parameter
           "assembly_ref" of type "Assembly_ref" (Reference to an Assembly
           object in the workspace @id ws KBaseGenomeAnnotations.Assembly),
           parameter "taxon_ref" of type "Taxon_ref" (Reference to a taxon
           object @id ws KBaseGenomeAnnotations.Taxon), parameter
           "genbank_handle_ref" of type "genbank_handle_ref" (Reference to a
           handle to the Genbank file on shock @id handle), parameter
           "gff_handle_ref" of type "gff_handle_ref" (Reference to a handle
           to the GFF file on shock @id handle), parameter
           "external_source_origination_date" of String, parameter "release"
           of String, parameter "original_source_file_name" of String,
           parameter "notes" of String, parameter "quality_scores" of list of
           type "GenomeQualityScore" (Genome quality score Fields: method -
           string - TODO method_report_ref - string - TODO method_version -
           string - TODO score: string - TODO score_interpretation - string -
           TODO timestamp - string - TODO Score_interpretation -
           fraction_complete - controlled vocabulary managed by API @optional
           method_report_ref method_version) -> structure: parameter "method"
           of String, parameter "method_report_ref" of type
           "Method_report_ref" (Reference to a report object @id ws
           KBaseReport.Report), parameter "method_version" of String,
           parameter "score" of String, parameter "score_interpretation" of
           String, parameter "timestamp" of String, parameter "suspect" of
           type "Bool", parameter "genome_type" of String, parameter "hidden"
           of type "boolean" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "upgrade" of type "boolean" (A boolean - 0 for
           false, 1 for true. @range (0, 1))
        :returns: instance of type "SaveGenomeResult" -> structure: parameter
           "info" of type "object_info" (Information about an object,
           including user provided metadata. obj_id objid - the numerical id
           of the object. obj_name name - the name of the object. type_string
           type - the type of the object. timestamp save_date - the save date
           of the object. obj_ver ver - the version of the object. username
           saved_by - the user that saved or copied the object. ws_id wsid -
           the workspace containing the object. ws_name workspace - the
           workspace containing the object. string chsum - the md5 checksum
           of the object. int size - the size of the object in bytes.
           usermeta meta - arbitrary user-supplied metadata about the
           object.) -> tuple of size 11: parameter "objid" of type "obj_id"
           (The unique, permanent numerical ID of an object.), parameter
           "name" of type "obj_name" (A string used as a name for an object.
           Any string consisting of alphanumeric characters and the
           characters |._- that is not an integer is acceptable.), parameter
           "type" of type "type_string" (A type string. Specifies the type
           and its version in a single string in the format
           [module].[typename]-[major].[minor]: module - a string. The module
           name of the typespec containing the type. typename - a string. The
           name of the type as assigned by the typedef statement. major - an
           integer. The major version of the type. A change in the major
           version implies the type has changed in a non-backwards compatible
           way. minor - an integer. The minor version of the type. A change
           in the minor version implies that the type has changed in a way
           that is backwards compatible with previous type definitions. In
           many cases, the major and minor versions are optional, and if not
           provided the most recent version will be used. Example:
           MyModule.MyType-3.1), parameter "save_date" of type "timestamp" (A
           time in the format YYYY-MM-DDThh:mm:ssZ, where Z is either the
           character Z (representing the UTC timezone) or the difference in
           time to UTC in the format +/-HHMM, eg: 2012-12-17T23:24:06-0500
           (EST time) 2013-04-03T08:56:32+0000 (UTC time)
           2013-04-03T08:56:32Z (UTC time)), parameter "version" of Long,
           parameter "saved_by" of type "username" (Login name of a KBase
           user account.), parameter "wsid" of type "ws_id" (The unique,
           permanent numerical ID of a workspace.), parameter "workspace" of
           type "ws_name" (A string used as a name for a workspace. Any
           string consisting of alphanumeric characters and "_", ".", or "-"
           that is not an integer is acceptable. The name may optionally be
           prefixed with the workspace owner's user name and a colon, e.g.
           kbasetest:my_workspace.), parameter "chsum" of String, parameter
           "size" of Long, parameter "meta" of type "usermeta" (User provided
           metadata about an object. Arbitrary key-value pairs provided by
           the user.) -> mapping from String to String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN save_one_genome

        genome_interface = GenomeInterface(self.cfg)
        returnVal = genome_interface.save_one_genome(params)
        #END save_one_genome

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method save_one_genome return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def ws_obj_gff_to_genome(self, ctx, params):
        """
        This function takes in a workspace object of type KBaseGenomes.Genome or KBaseGenomeAnnotations.Assembly and a gff file and produces a KBaseGenomes.Genome reanotated according to the the input gff file.
        :param params: instance of type "WsObjGFFToGenomeParams" (gff_file -
           object containing path to gff_file ws_ref - input Assembly or
           Genome reference genome_name - becomes the name of the object
           workspace_name - the name of the workspace it gets saved to.
           source - Source of the file typically something like RefSeq or
           Ensembl taxon_ws_name - where the reference taxons are :
           ReferenceTaxons taxon_id - if defined, will try to link the Genome
           to the specified taxonomy id in lieu of performing the lookup
           during upload release - Release or version number of the data per
           example Ensembl has numbered releases of all their data: Release
           31 genetic_code - Genetic code of organism. Overwrites determined
           GC from taxon object scientific_name - will be used to set the
           scientific name of the genome and link to a taxon metadata - any
           user input metadata generate_missing_genes - If the file has CDS
           or mRNA with no corresponding gene, generate a spoofed gene. Off
           by default) -> structure: parameter "ws_ref" of String, parameter
           "gff_file" of type "File" -> structure: parameter "path" of
           String, parameter "shock_id" of String, parameter "ftp_url" of
           String, parameter "genome_name" of String, parameter
           "workspace_name" of String, parameter "source" of String,
           parameter "taxon_wsname" of String, parameter "taxon_id" of
           String, parameter "release" of String, parameter "genetic_code" of
           Long, parameter "scientific_name" of String, parameter "metadata"
           of type "usermeta" -> mapping from String to String, parameter
           "generate_missing_genes" of type "boolean" (A boolean - 0 for
           false, 1 for true. @range (0, 1))
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "genome_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN ws_obj_gff_to_genome
        au = AssemblyUtil(self.cfg.callbackURL)
        ws = Workspace(url=self.cfg.workspaceURL)

        if params.get('ws_ref'):
            input_ref = params['ws_ref']
        else:
            raise ValueError(f"ws_ref argument required for ws_obj_gff_to_metagenome method") 
        # input should only be of workspace object types:
        #   - 'KBaseGenomes.Genome'
        #   - 'KBaseGenomeAnnotations.Assembly'
        #   - 'KBaseGenomes.ContigSet'
        obj_type = ws.get_object_info3({'objects': [{'ref': input_ref}]})['infos'][0][2]
        obj_type = obj_type.split('-')[0]
        valid_types = ["KBaseGenomes.Genome", "KBaseGenomes.ContigSet", "KBaseGenomeAnnotations.Assembly"]
        if obj_type not in valid_types:
            raise ValueError(f"input reference must be of type KBaseGenomes.Genome or "
                             f"KBaseGenomes.ContigSet or KBaseGenomeAnnotations.Assembly "
                             f"not {obj_type}")

        fasta_file = [val['paths'][0] for key, val in au.get_fastas({'ref_lst': [input_ref]}).items()][0]

        params['fasta_file'] = {
            'path': fasta_file
        }

        if obj_type == "KBaseGenomes.Genome":
            assembly_ref = ws.get_objects2(
                {"objects":
                    [{"ref": input_ref, "included": ["assembly_ref"]}]
                }
            )['data'][0]['data']['assembly_ref']
            params['existing_assembly_ref'] = assembly_ref
        else:
            params['existing_assembly_ref'] = input_ref

        importer = FastaGFFToGenome(self.cfg)
        returnVal = importer.import_file(params)

        #END ws_obj_gff_to_genome

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method ws_obj_gff_to_genome return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def ws_obj_gff_to_metagenome(self, ctx, params):
        """
        This function takes in a workspace object of type KBaseMetagenomes.AnnotatedMetagenomeAssembly or KBaseGenomeAnnotations.Assembly and a gff file and produces a KBaseMetagenomes.AnnotatedMetagenomeAssembly reanotated according to the the input gff file.
        :param params: instance of type "WsObjGFFToMetagenomeParams"
           (gff_file - object containing path to gff_file ws_ref - input
           Assembly or AnnotatedMetagenomeAssembly reference genome_name -
           becomes the name of the object workspace_name - the name of the
           workspace it gets saved to. source - Source of the file typically
           something like RefSeq or Ensembl genetic_code - Genetic code of
           organism. Overwrites determined GC from taxon object metadata -
           any user input metadata generate_missing_genes - If the file has
           CDS or mRNA with no corresponding gene, generate a spoofed gene.
           Off by default) -> structure: parameter "ws_ref" of String,
           parameter "gff_file" of type "File" -> structure: parameter "path"
           of String, parameter "shock_id" of String, parameter "ftp_url" of
           String, parameter "genome_name" of String, parameter
           "workspace_name" of String, parameter "source" of String,
           parameter "metadata" of type "usermeta" -> mapping from String to
           String, parameter "generate_missing_genes" of type "boolean" (A
           boolean - 0 for false, 1 for true. @range (0, 1))
        :returns: instance of type "MetagenomeSaveResult" -> structure:
           parameter "metagenome_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN ws_obj_gff_to_metagenome
        au = AssemblyUtil(self.cfg.callbackURL)
        ws = Workspace(url=self.cfg.workspaceURL)

        if params.get('ws_ref'):
            input_ref = params['ws_ref']
        else:
            raise ValueError(f"ws_ref argument required for ws_obj_gff_to_metagenome method")
        # input should only be of workspace object types:
        #   - 'KBaseMetagenomes.AnnotatedMetagenomeAssembly'
        #   - 'KBaseGenomeAnnotations.Assembly'
        #   - 'KBaseGenomes.ContigSet'
        obj_type = ws.get_object_info3({'objects': [{'ref': input_ref}]})['infos'][0][2]
        obj_type = obj_type.split('-')[0]
        valid_types = ["KBaseMetagenomes.AnnotatedMetagenomeAssembly", "KBaseGenomes.ContigSet", "KBaseGenomeAnnotations.Assembly"]
        if obj_type not in valid_types:
            raise ValueError(f"input reference must be of type KBaseMetagenomes.AnnotatedMetagenomeAssembly or "
                             f"KBaseGenomes.ContigSet or KBaseGenomeAnnotations.Assembly "
                             f"not {obj_type}")

        fasta_file = [val['paths'][0] for key, val in au.get_fastas({'ref_lst': [input_ref]}).items()][0]
        params['fasta_file'] = {
          'path': fasta_file
        }
        if obj_type == "KBaseMetagenomes.AnnotatedMetagenomeAssembly":
            assembly_ref = ws.get_objects2(
                {"objects":
                    [{"ref": input_ref, "included": ["assembly_ref"]}]
                }
            )['data'][0]['data']['assembly_ref']
            params['existing_assembly_ref'] = assembly_ref
        else:
            params['existing_assembly_ref'] = input_ref
        params['is_metagenome'] = True

        importer = FastaGFFToGenome(self.cfg)
        returnVal = importer.import_file(params)
        #END ws_obj_gff_to_metagenome

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method ws_obj_gff_to_metagenome return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def update_taxon_assignments(self, ctx, params):
        """
        Add, replace, or remove taxon assignments for a Genome object.
        :param params: instance of type "UpdateTaxonAssignmentsParams"
           (Parameters for the update_taxon_assignments function. Fields:
           workspace_id: a workspace UPA of a Genome object
           taxon_assignments: an optional mapping of assignments to add or
           replace. This will perform a merge on the existing assignments.
           Any new assignments are added, while any existing assignments are
           replaced. remove_assignments: an optional list of assignment names
           to remove. @optional taxon_assignments remove_assignments) ->
           structure: parameter "workspace_id" of Long, parameter "object_id"
           of Long, parameter "taxon_assignments" of mapping from String to
           String, parameter "remove_assignments" of list of String
        :returns: instance of type "UpdateTaxonAssignmentsResult" (Result of
           the update_taxon_assignments function. Fields: ws_obj_ref: a
           workspace UPA of a Genome object) -> structure: parameter
           "ws_obj_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN update_taxon_assignments
        for key in ('workspace_id', 'object_id'):
            if not params.get(key):
                raise TypeError(f"`{key}` is required but missing.")
            if not isinstance(params[key], int):
                raise TypeError(f"`{key}` must be an integer.")
        dfu = DataFileUtil(self.cfg.callbackURL)
        obj_ref = f"{params['workspace_id']}/{params['object_id']}"
        result = dfu.get_objects({'object_refs': [obj_ref]})['data'][0]
        obj_data = result['data']
        obj_info = result['info']
        if 'taxon_assignments' not in obj_data:
            obj_data['taxon_assignments'] = {}
        # Merge in params['taxon_assignments']
        for key, val in params.get('taxon_assignments', {}).items():
            obj_data['taxon_assignments'][key] = val
        # Remove all keys from `remove_assignments`
        for key in params.get('remove_assignments', []):
            if key in obj_data['taxon_assignments']:
                del obj_data['taxon_assignments'][key]
        new_obj = {
            'type': obj_info[2],
            'data': obj_data,
            'objid': obj_info[0]
        }
        infos = dfu.save_objects({
            'id': params['workspace_id'],
            'objects': [new_obj]
        })
        obj_ref = f"{infos[0][6]}/{infos[0][0]}/{infos[0][4]}"
        returnVal = {'ws_obj_ref': obj_ref}
        #END update_taxon_assignments

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method update_taxon_assignments return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION,
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
