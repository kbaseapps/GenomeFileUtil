# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
import sys
import shutil
import traceback
import uuid
from pprint import pprint, pformat
import re

from GenomeFileUtil.core.GenbankToGenome import GenbankToGenome
from GenomeFileUtil.core.GenomeToGFF import GenomeToGFF

from biokbase.workspace.client import Workspace as workspaceService
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI as GenomeAnnotationAPI

# silence whining
import requests
requests.packages.urllib3.disable_warnings()


# Used to store and pass around configuration URLs more easily
class SDKConfig:
    def __init__(self, config):
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.sharedFolder = config['scratch']

#END_HEADER


class GenomeFileUtil:
    '''
    Module Name:
    GenomeFileUtil

    Module Description:
    ** A KBase module: GenomeFileUtil
**
** This module contains methods for converting KBase Data Objects to common bioinformatics file formats
**
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.1.0"
    GIT_URL = "https://github.com/dcchivian/GenomeFileUtil"
    GIT_COMMIT_HASH = "aff5e42bff205da10c3a2a2711f6f60d82d12d9b"
    
    #BEGIN_CLASS_HEADER

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

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
        :param params: instance of type "GenbankToGenomeParams" -> structure:
           parameter "file" of type "File" -> structure: parameter "path" of
           String, parameter "shock_id" of String, parameter "ftp_url" of
           String, parameter "genome_name" of String, parameter
           "workspace_name" of String, parameter "source" of String,
           parameter "taxon_wsname" of String
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

    def genome_to_gff(self, ctx, params):
        """
        :param params: instance of type "GenomeToGFFParams" -> structure:
           parameter "genome_ref" of String, parameter "ref_path_to_genome"
           of list of String
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

    def TranslateNucToProtSeq(self, ctx, params):
        """
        :param params: instance of type "TranslateNucToProtSeq_Params"
           (TranslateNucToProtSeq()) -> structure: parameter "nuc_seq" of
           String, parameter "genetic_code" of String
        :returns: instance of type "TranslateNucToProtSeq_Output" ->
           structure: parameter "prot_seq" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN TranslateNucToProtSeq
        if 'nuc_seq' != params or params['nuc_seq'] == None:
            raise ValueError('Method TranslateNucToProtSeq() requires nuc_seq parameter')
        if 'genetic_code' != params or params['genetic_code'] == None:
            params['genetic_code'] = '11'

        if params['genetic_cde'] != '11':
            raise ValueError('Method TranslateNucToProtSeq() only knows genetic code 11')
        
        nuc_seq = nuc_seq.upper()
        prot_seq = ''

        genetic_code = dict()
        genetic_code['11'] = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
            }
        prot_seq = ''.join([genetic_code[table].get(nuc_seq[3*i:3*i+3],'X') for i in range(len(nuc_seq)//3)])

        returnVal = dict()
        returnVal['prot_seq'] = prot_seq
        #END TranslateNucToProtSeq

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method TranslateNucToProtSeq return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def GenomeToFASTA(self, ctx, params):
        """
        :param params: instance of type "GenomeToFASTA_Params"
           (GenomeToFASTA()) -> structure: parameter "genome_ref" of type
           "data_obj_ref", parameter "file" of type "path_type", parameter
           "dir" of type "path_type", parameter "console" of list of type
           "log_msg", parameter "invalid_msgs" of list of type "log_msg",
           parameter "residue_type" of String, parameter "feature_type" of
           String, parameter "record_id_pattern" of type "pattern_type",
           parameter "record_desc_pattern" of type "pattern_type", parameter
           "case" of String, parameter "linewrap" of Long
        :returns: instance of type "GenomeToFASTA_Output" -> structure:
           parameter "fasta_file_path" of type "path_type", parameter
           "feature_ids" of list of type "feature_id"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeToFASTA

        # init and clean up params
        genome_ref = params['genome_ref']
        file = params['file']
        dir = params['dir']
        console = params['console']
        invalid_msgs = params['invalid_msgs']
        residue_type = params['residue_type']
        feature_type = params['feature_type']
        record_id_pattern = params['record_id_pattern']
        record_desc_pattern = params['record_desc_pattern']
        case = params['case']
        linewrap = params['linewrap']

        # Defaults
        if console == None:
            console = []
        if invalid_msgs == None:
            invalid_msgs = []
        if residue_type == None:
            residue_type = 'nuc'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = '%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_ref%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0

        # init and simplify
        feature_ids = []
        feature_sequence_found = False
        residue_type = residue_type[0:3].lower()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        
        def record_header_sub(str, feature_id, genome_id, genome_ref):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            str = str.replace('%%genome_ref%%', genome_ref)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch
        fasta_file_path = os.path.join(dir, file)
        self.log(console, 'KB SDK data2file GenomeAnnotationToFASTA(): writing fasta file: '+fasta_file_path)

        # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Perhaps compromise with limited buffer size that flushes as we go.  Check later.

        #records = []
        self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG

        with open(fasta_file_path, 'w', 0) as fasta_file_handle:

            GA = GenomeAnnotationAPI ({"workspace_service_url": self.workspaceURL,
                                       "shock_service_url": self.shockURL
                                       },
                                      token=os.environ["KB_AUTH_TOKEN"],
                                      ref=genome_ref
                                     )
            features = GA.get_features()
            if residue_type == 'pro' or residue_type == 'pep':
                proteins = GA.get_proteins()

            for fid in features.keys():
                feature = features[fid]
                
                if feature_type == 'ALL' or feature_type == feature['feature_type']:

                    # protein recs
                    if residue_type == 'pro' or residue_type == 'pep':
                        if feature['feature_type'] != 'CDS':
                            continue
                        elif fid not in proteins or 'protein_amino_acid_sequence' not in proteins[fid] or proteins[fid]['protein_amino_acid_sequence'] == None:
                            self.log(invalid_msgs, "bad CDS feature "+fid+": No protein_translation field.")
                        else:
                            feature_sequence_found = True
                            # can't get genome_id from GenomeAnnotationAPI
                            #rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            #rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
                            rec_id = record_header_sub(record_id_pattern, fid, genome_ref, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_ref, genome_ref)
                            seq = proteins[fid]['protein_amino_acid_sequence']
                            seq = seq.upper() if case == 'U' else seq.lower()

                            rec_rows = []
                            rec_rows.append('>'+rec_id+' '+rec_desc)
                            if linewrap == None or linewrap == 0:
                                rec_rows.append(seq)
                            else:
                                seq_len = len(seq)
                                base_rows_cnt = seq_len//linewrap
                                for i in range(base_rows_cnt):
                                    rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                rec_rows.append(seq[base_rows_cnt*linewrap:])
                            rec = "\n".join(rec_rows)+"\n"

                            feature_ids.append(fid)
                            fasta_file_handle.write(rec)

                    # nuc recs
                    else:
                        if 'dna_sequence' not in feature or feature['dna_sequence'] == None:
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        else:
                            feature_sequence_found = True
                            # can't get genome_id from GenomeAnnotationAPI
                            #rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            #rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
                            rec_id = record_header_sub(record_id_pattern, fid, genome_ref, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_ref, genome_ref)
                            seq = feature['dna_sequence']
                            seq = seq.upper() if case == 'U' else seq.lower()

                            rec_rows = []
                            rec_rows.append('>'+rec_id+' '+rec_desc)
                            if linewrap == None or linewrap == 0:
                                rec_rows.append(seq)
                            else:
                                seq_len = len(seq)
                                base_rows_cnt = seq_len//linewrap
                                for i in range(base_rows_cnt):
                                    rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                rec_rows.append(seq[base_rows_cnt*linewrap:])
                            rec = "\n".join(rec_rows)+"\n"

                            feature_ids.append(fid)
                            fasta_file_handle.write(rec)

        # report if no features found
        if not feature_sequence_found:
            self.log(invalid_msgs, "No sequence records found in Genome "+genome_object['id']+" of residue_type: "+residue_type+", feature_type: "+feature_type)
        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path'] = fasta_file_path
        returnVal['feature_ids'] = feature_ids
        #END GenomeToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GenomeToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def GenomeSetToFASTA(self, ctx, params):
        """
        :param params: instance of type "GenomeSetToFASTA_Params"
           (GenomeSetToFASTA()) -> structure: parameter "genomeSet_ref" of
           type "data_obj_ref", parameter "file" of type "path_type",
           parameter "dir" of type "path_type", parameter "console" of list
           of type "log_msg", parameter "invalid_msgs" of list of type
           "log_msg", parameter "residue_type" of String, parameter
           "feature_type" of String, parameter "record_id_pattern" of type
           "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long, parameter "merge_fasta_files" of type "true_false"
        :returns: instance of type "GenomeSetToFASTA_Output" -> structure:
           parameter "fasta_file_path_list" of list of type "path_type",
           parameter "feature_ids_by_genome_id" of mapping from type
           "genome_id" to list of type "feature_id"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeSetToFASTA

        # init and clean up params
        genomeSet_ref = params['genomeSet_ref']
        file = params['file']
        dir = params['dir']
        console = params['console']
        invalid_msgs = params['invalid_msgs']
        residue_type = params['residue_type']
        feature_type = params['feature_type']
        record_id_pattern = params['record_id_pattern']
        record_desc_pattern = params['record_desc_pattern']
        case = params['case']
        linewrap = params['linewrap']
        merge_fasta_files = params['merge_fasta_files']

        # Defaults
        if console == None:
            console = []
        if invalid_msgs == None:
            invalid_msgs = []
        if residue_type == None:
            residue_type = 'nuc'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = 'g:%%genome_ref%%.f:%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_ref%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0
        if merge_fasta_files == None or merge_fasta_files[0:1].upper() == 'T':
            merge_fasta_files = True
        else:
            merge_fasta_files = False

        # init and simplify
        fasta_file_path_list = []
        feature_ids_by_genome_id = dict()
        feature_sequence_found = False
        residue_type = residue_type[0:3].lower()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        
        def record_header_sub(str, feature_id, genome_id, genome_ref):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            str = str.replace('%%genome_ref%%', genome_ref)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch

        # get genomeSet object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            genomeSet_object = ws.get_objects2({'objects':[{'ref':genomeSet_ref}]})['data'][0]
        except Exception as e:
            raise ValueError('Unable to fetch GenomeSet object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # iterate through genomeSet members
        genome_ids = genomeSet_object['elements'].keys()
        for genome_i in range(len(genome_ids)):
            genome_id = genome_ids[genome_i]
            feature_ids_by_genome_id[genome_id] = []

            if 'ref' not in genomeSet_object['elements'][genome_id] or \
                    genomeSet_object['elements'][genome_id]['ref'] == None:
                raise ValueError('GenomeSetToFASTA() cannot handle GenomeSet objects with embedded genome.  Must be a set of genome references')
                #to get the full stack trace: traceback.format_exc()       
            else:
                genome_ref = genomeSet_object['elements'][genome_id]['ref']

            # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Perhaps compromise with limited buffer size that flushes as we go.  Check later.

            #records = []
            this_file = file
            if len(genome_ids) > 1 and not merge_fasta_files:
                this_file = this_file+'.'+genome_id
                this_file.replace('/', '-')
            fasta_file_path = os.path.join(dir, this_file)
            #self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG
            
            if genome_i == 0 or not merge_fasta_files:
                fasta_file_handle = open(fasta_file_path, 'w', 0)
                self.log(console, 'KB SDK data2file GenomeSet2FASTA(): writing fasta file: '+fasta_file_path)
            # DEBUG
            #self.log(console, "ADDING GENOME: "+str(genome_i+1)+" of "+str(len(genome_ids))+" "+genome_ids[genome_i])

            GA = GenomeAnnotationAPI ({"workspace_service_url": self.workspaceURL,
                                       "shock_service_url": self.shockURL
                                       },
                                      token=os.environ["KB_AUTH_TOKEN"],
                                      ref=genome_ref
                                      )
            features = GA.get_features()
            if residue_type == 'pro' or residue_type == 'pep':
                proteins = GA.get_proteins()

            for fid in features.keys():
                feature = features[fid]
                
                if feature_type == 'ALL' or feature_type == feature['feature_type']:
                    
                    # protein recs
                    if residue_type == 'pro' or residue_type == 'pep':
                        if feature['feature_type'] != 'CDS':
                            continue
                        elif fid not in proteins or 'protein_amino_acid_sequence' not in proteins[fid] or proteins[fid]['protein_amino_acid_sequence'] == None:
                            self.log(invalid_msgs, "bad CDS feature "+fid+": No protein_translation field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
                            seq = proteins[fid]['protein_amino_acid_sequence']
                            seq = seq.upper() if case == 'U' else seq.lower()

                            rec_rows = []
                            rec_rows.append('>'+rec_id+' '+rec_desc)
                            if linewrap == None or linewrap == 0:
                                rec_rows.append(seq)
                            else:
                                seq_len = len(seq)
                                base_rows_cnt = seq_len//linewrap
                                for i in range(base_rows_cnt):
                                    rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                rec_rows.append(seq[base_rows_cnt*linewrap:])
                            rec = "\n".join(rec_rows)+"\n"

                            feature_ids_by_genome_id[genome_id].append(fid)
                            fasta_file_handle.write(rec)

                    # nuc recs
                    else:
                        if 'dna_sequence' not in feature or feature['dna_sequence'] == None:
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
                            seq = feature['dna_sequence']
                            seq = seq.upper() if case == 'U' else seq.lower()
                            
                            rec_rows = []
                            rec_rows.append('>'+rec_id+' '+rec_desc)
                            if linewrap == None or linewrap == 0:
                                rec_rows.append(seq)
                            else:
                                seq_len = len(seq)
                                base_rows_cnt = seq_len//linewrap
                                for i in range(base_rows_cnt):
                                    rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                rec_rows.append(seq[base_rows_cnt*linewrap:])
                            rec = "\n".join(rec_rows)+"\n"

                            feature_ids_by_genome_id[genome_id].append(fid)
                            fasta_file_handle.write(rec)

            if genome_i == (len(genome_ids)-1) or not merge_fasta_files:
                #self.log(console,"CLOSING FILE: '"+fasta_file_path+"'")  # DEBUG
                fasta_file_handle.close()
                fasta_file_path_list.append(fasta_file_path)


        # report if no features found
        if not feature_sequence_found:
            self.log(invalid_msgs, "No sequence records found in Genome "+genome_object['id']+" of residue_type: "+residue_type+", feature_type: "+feature_type)

        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path_list'] = fasta_file_path_list
        returnVal['feature_ids_by_genome_id'] = feature_ids_by_genome_id
        #END GenomeSetToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GenomeSetToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def FeatureSetToFASTA(self, ctx, params):
        """
        :param params: instance of type "FeatureSetToFASTA_Params"
           (FeatureSetToFASTA()) -> structure: parameter "featureSet_ref" of
           type "data_obj_ref", parameter "file" of type "path_type",
           parameter "dir" of type "path_type", parameter "console" of list
           of type "log_msg", parameter "invalid_msgs" of list of type
           "log_msg", parameter "residue_type" of String, parameter
           "feature_type" of String, parameter "record_id_pattern" of type
           "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long
        :returns: instance of type "FeatureSetToFASTA_Output" -> structure:
           parameter "fasta_file_path" of type "path_type", parameter
           "feature_ids_by_genome_ref" of mapping from type "data_obj_ref" to
           list of type "feature_id"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN FeatureSetToFASTA

        # init and clean up params
        featureSet_ref = params['featureSet_ref']
        file = params['file']
        dir = params['dir']
        console = params['console']
        invalid_msgs = params['invalid_msgs']
        residue_type = params['residue_type']
        feature_type = params['feature_type']
        record_id_pattern = params['record_id_pattern']
        record_desc_pattern = params['record_desc_pattern']
        case = params['case']
        linewrap = params['linewrap']

        # Defaults
        if console == None:
            console = []
        if invalid_msgs == None:
            invalid_msgs = []
        if residue_type == None:
            residue_type = 'nuc'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = 'g:%%genome_ref%%.f:%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_ref%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0

        # init and simplify
        feature_ids_by_genome_ref = dict()
        feature_sequence_found = False
        residue_type = residue_type[0:3].lower()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        
        def record_header_sub(str, feature_id, genome_id, genome_ref):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            str = str.replace('%%genome_ref%%', genome_ref)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch
        fasta_file_path = os.path.join(dir, file)
        self.log(console, 'KB SDK data2file FeatureSetToFASTA(): writing fasta file: '+fasta_file_path)


        # get featureSet object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            featureSet_object = ws.get_objects2({'objects':[{'ref':featureSet_ref}]})['data'][0]
        except Exception as e:
            raise ValueError('Unable to fetch featureSet object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()
        featureSet_features = featureSet_object['elements']
        genome2features = {}
        for fId in featureSet_features.keys():
            for genome_ref in featureSet_features[fId]:
                genome_ref = featureSet_features[fId][0]
                if genome_ref not in genome2features.keys():
                    genome2features[genome_ref] = []
                genome2features[genome_ref].append(fId)

        # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Perhaps compromise with limited buffer size that flushes as we go.  Check later.

        #records = []
        #self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG

        with open(fasta_file_path, 'w', 0) as fasta_file_handle:

            for genome_ref in genome2features.keys():
                feature_ids_by_genome_ref[genome_ref] = []

                GA = GenomeAnnotationAPI ({"workspace_service_url": self.workspaceURL,
                                           "shock_service_url": self.shockURL
                                           },
                                          token=os.environ["KB_AUTH_TOKEN"],
                                          ref=genome_ref
                                          )
                features = GA.get_features(feature_id_list=genome2features[genome_ref])
                if residue_type == 'pro' or residue_type == 'pep':
                    # feature_id_list not working for get_proteins()
                    #proteins = GA.get_proteins(feature_id_list=genome2features[genome_ref])
                    proteins = GA.get_proteins()

#                for feature in genome_object['features']:
                for fid in genome2features[genome_ref]:
                    feature = features[fid]
                
                    if feature_type == 'ALL' or feature_type == feature['feature_type']:

                        # protein recs
                        if residue_type == 'pro' or residue_type == 'pep':
                            if feature['feature_type'] != 'CDS':
                                continue
                            elif fid not in proteins or 'protein_amino_acid_sequence' not in proteins[fid] or proteins[fid]['protein_amino_acid_sequence'] == None:
                                self.log(invalid_msgs, "bad CDS feature "+fid+": No protein_translation field.")
                            else:
                                feature_sequence_found = True
                                #rec_id = record_header_sub(record_id_pattern, fid, genome_ref, genome_ref)
                                #rec_desc = record_header_sub(record_desc_pattern, fid, genome_ref, genome_ref)
                                rec_id = record_header_sub(record_id_pattern, fid, genome_ref, genome_ref)
                                rec_desc = record_header_sub(record_desc_pattern, fid, genome_ref, genome_ref)
                                seq = proteins[fid]['protein_amino_acid_sequence']
                                seq = seq.upper() if case == 'U' else seq.lower()

                                rec_rows = []
                                rec_rows.append('>'+rec_id+' '+rec_desc)
                                if linewrap == None or linewrap == 0:
                                    rec_rows.append(seq)
                                else:
                                    seq_len = len(seq)
                                    base_rows_cnt = seq_len//linewrap
                                    for i in range(base_rows_cnt):
                                        rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                    rec_rows.append(seq[base_rows_cnt*linewrap:])
                                rec = "\n".join(rec_rows)+"\n"

                                feature_ids_by_genome_ref[genome_ref].append(fid)
                                fasta_file_handle.write(rec)

                        # nuc recs
                        else:
                            if 'dna_sequence' not in feature or feature['dna_sequence'] == None:
                                self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                            else:
                                feature_sequence_found = True
                                rec_id = record_id_pattern
                                rec_desc = record_desc_pattern
                                rec_id = record_header_sub(rec_id, fid, genome_ref)
                                rec_desc = record_header_sub(rec_desc, fid, genome_ref)
                                seq = feature['dna_sequence']
                                seq = seq.upper() if case == 'U' else seq.lower()

                                rec_rows = []
                                rec_rows.append('>'+rec_id+' '+rec_desc)
                                if linewrap == None or linewrap == 0:
                                    rec_rows.append(seq)
                                else:
                                    seq_len = len(seq)
                                    base_rows_cnt = seq_len//linewrap
                                    for i in range(base_rows_cnt):
                                        rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                    rec_rows.append(seq[base_rows_cnt*linewrap:])
                                rec = "\n".join(rec_rows)+"\n"

                                feature_ids_by_genome_ref[genome_ref].append(fid)
                                fasta_file_handle.write(rec)

        # report if no features found
        if not feature_sequence_found:
            self.log(invalid_msgs, "No sequence records found in FeatureSet "+featureSet_ref+" of residue_type: "+residue_type+", feature_type: "+feature_type)

        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path'] = fasta_file_path
        returnVal['feature_ids_by_genome_ref'] = feature_ids_by_genome_ref
        #END FeatureSetToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method FeatureSetToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
