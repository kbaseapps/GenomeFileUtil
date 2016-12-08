#!/usr/bin/env python

# standard library imports
import os
import sys
import logging
import re
import hashlib
import time 
import traceback 
import os.path 
import datetime
import shutil
#import sqlite3 
#try: 
#    import cPickle as cPickle 
#except: 
#    import pickle as cPickle 
from string import digits
from string import maketrans
from collections import OrderedDict

try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

# 3rd party imports
import simplejson
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

# KBase imports
import biokbase.Transform.script_utils as script_utils
import biokbase.Transform.TextFileDecoder as TextFileDecoder
import biokbase.workspace.client 
import trns_transform_FASTA_DNA_Assembly_to_KBaseGenomeAnnotations_Assembly as assembly
from KBaseReport.KBaseReportClient import KBaseReport
#from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI, GenomeAnnotationClientAPI

def insert_newlines(s, every): 
    lines = [] 
    for i in xrange(0, len(s), every): 
        lines.append(s[i:i+every]) 
    return "\n".join(lines)+"\n" 

def represents_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


# transformation method that can be called if this module is imported
# Note the logger has different levels it could be run.  
# See: https://docs.python.org/2/library/logging.html#logging-levels
#
# The default level is set to INFO which includes everything except DEBUG
#@profile
def upload_genome(shock_service_url=None,
                  handle_service_url=None,
                  workspace_service_url=None,
                  callback_url=None,

                  input_directory=None, 

                  shock_id=None, handle_id=None,

                  workspace_name=None,
                  core_genome_name=None,

                  taxon_wsname=None,
                  taxon_lookup_obj_name=None,
                  taxon_reference = None,

                  exclude_ontologies=None,
                  ontology_wsname=None,
                  ontology_GO_obj_name = None,
                  ontology_PO_obj_name = None,

                  release= None,
                  source=None,
                  type=None,
                  genetic_code=None,
                  generate_ids_if_needed=None,

                  provenance=None,
                  usermeta=None,

                  level=logging.INFO, logger=None):
    """
    Uploads CondensedGenomeAssembly
    Args:
        shock_service_url: A url for the KBase SHOCK service.
        input_fasta_directory: The directory where files will be read from.
        level: Logging level, defaults to logging.INFO.
        
    Returns:
        JSON file on disk that can be saved as a KBase workspace object.
    Authors:
        Jason Baumohl, Matt Henderson
    """

    if logger is None:
        logger = script_utils.stderrlogger(__file__)
    token = os.environ.get('KB_AUTH_TOKEN') 

    logger.info("WS URL : " + str(workspace_service_url))
    logger.info("HANDLE URL : " + str(handle_service_url))
    logger.info("SHOCK URL : " + str(shock_service_url))
    logger.info("CALLBACK URL : " + str(callback_url))

    ws_client = biokbase.workspace.client.Workspace(workspace_service_url)
 
    workspace_object = ws_client.get_workspace_info({'workspace':workspace_name}) 

    workspace_id = workspace_object[0] 
    workspace_name = workspace_object[1] 
 
    taxon_workspace_object = ws_client.get_workspace_info({'workspace':taxon_wsname}) 

    taxon_workspace_id = taxon_workspace_object[0] 
    taxon_workspace_name = taxon_workspace_object[1] 

    report = StringIO() #variable to put warnings report into.  For UI widget reports output tab

    if exclude_ontologies is not None:
        if exclude_ontologies != 1:
            exclude_ontologies == 0
    #Get GO OntologyDictionary
    ontology_sources = dict()

    if exclude_ontologies == 0:
        #    ontologies = ws_client.get_objects2({'objects': [{'workspace': 'KBaseOntology', 'name':'gene_ontology'}]}) 
        #    go_ontology = ontologies['data'][0]['data'] 
        logger.info("Retrieving Ontology databases.") 
        ontologies = ws_client.get_objects( [{'workspace':ontology_wsname,
                                              'name':ontology_GO_obj_name},
                                             {'workspace':ontology_wsname,
                                              'name':ontology_PO_obj_name}])
        logger.info("Ontology databases retrieved.") 
        
        ontology_sources["GO"] = ontologies[0]['data']['term_hash']
        ontology_sources["PO"] = ontologies[1]['data']['term_hash']
        del ontologies
#    go_ontologies = ws_client.get_objects( [{'workspace':'KBaseOntology',
#                                             'name':'gene_ontology'}])
#    logger.info("Retrieved GO Ontology database, starting PO") 
#    po_ontologies = ws_client.get_objects( [{'workspace':'KBaseOntology',
#                                          'name':'plant_ontology'}])
#    logger.info("Retrieved PO Ontology database") 
#    ontology_sources["GO"] = go_ontologies[0]['data']['term_hash']
#    ontology_sources["PO"] = po_ontologies[0]['data']['term_hash']
#    del go_ontologies
#    del po_ontologies

    logger.info("Scanning for Genbank Format files.") 
    logger.info("GENETIC_CODE ENTERED : {}".format(str(genetic_code)))
    valid_extensions = [".gbff",".gbk",".gb",".genbank",".dat"] 
 
    files = os.listdir(os.path.abspath(input_directory)) 
    print "FILES : " + str(files)
    genbank_files = [x for x in files if os.path.splitext(x)[-1] in valid_extensions] 

    genetic_code_supplied = False
    if genetic_code is not None:
        genetic_code_supplied = True
        valid_genetic_codes = [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26]
        if genetic_code not in valid_genetic_codes:
            raise Exception("The entered genetic code of {} is not a valid genetic code, please see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi".format(str(genetic_code)))
    else:
        genetic_code = 1


    if generate_ids_if_needed is not None:
        if generate_ids_if_needed != 1:
            generate_ids_if_needed = 0
    else:
        generate_ids_if_needed = 0

    if (len(genbank_files) == 0): 
        raise Exception("The input directory does not have one of the following extensions %s." % (",".join(valid_extensions))) 
  
    logger.info("Found {0}".format(str(genbank_files))) 
 
    source_file_name = genbank_files[0]
    input_file_name = os.path.join(input_directory,genbank_files[0]) 
 
    if len(genbank_files) > 1: 
        # TODO if multiple files - CONCATENATE FILES HERE (sort by name)? OR Change how the byte coordinates work.
        report.write("Not sure how to handle multiple Genbank files in this context. Using {0}.\n\n".format(input_file_name))
        logger.warning("Not sure how to handle multiple Genbank files in this context. Using {0}".format(input_file_name))

    print "INPUT FILE NAME :" + input_file_name + ":"

    genbank_file_boundaries = list()  
    #list of tuples: (first value record start byte position, second value record stop byte position)

    if os.path.isfile(input_file_name):
        print "Found Genbank_File" 
#        make_sql_in_memory = True
        dir_name = os.path.dirname(input_file_name)

        #take in Genbank file and remove all empty lines from it.
        os.rename(input_file_name,"%s/temp_file_name" % (dir_name))
        temp_file = "%s/temp_file_name" % (dir_name)

        with open(temp_file,'r') as f_in:
            with open(input_file_name,'w', buffering=2**20 ) as f_out:
                for line in f_in:
                    if line.strip():
                        f_out.write(line)
        os.remove(temp_file)

        #If file is over a 1GB need to do SQLLite on disc
#        if os.stat(input_file_name) > 1073741824 :
#            make_sql_in_memory = True
            
        genbank_file_handle = TextFileDecoder.open_textdecoder(input_file_name, 'ISO-8859-1') 
        start_position = 0
        current_line = genbank_file_handle.readline()
        last_line = None
        while (current_line != ''):
            last_line = current_line 
            if current_line.startswith("//"):
                end_position =  genbank_file_handle.tell() - len(current_line)
                genbank_file_boundaries.append([start_position,end_position])
#                last_start_position = start_position
                start_position = genbank_file_handle.tell()
            current_line = genbank_file_handle.readline()

        if not last_line.startswith("//"):
            end_position = genbank_file_handle.tell()
            genbank_file_boundaries.append([start_position,end_position])
    else:
        raise ValueError("NO GENBANK FILE")

    print "Number of contigs : " + str(len(genbank_file_boundaries))
   
    organism_dict = dict() 
    organism = None
    if len(genbank_file_boundaries) < 1 :
        raise ValueError("Error no genbank record found in the input file")
    else:
        byte_coordinates = genbank_file_boundaries[0]
        genbank_file_handle.seek(byte_coordinates[0]) 
        temp_record = genbank_file_handle.read(byte_coordinates[1] - byte_coordinates[0]) 

        record_lines = temp_record.split("\n")
        for record_line in record_lines:
            if record_line.startswith("  ORGANISM  "):
                organism = record_line[12:]
                print "Organism Line :" + record_line + ":"
                print "Organism :" + organism + ":"
                organism_dict[organism] = 1
                break

    tax_id = 0;
    tax_lineage = None;

    genome = {
        'notes':''
    }

    display_sc_name = None

    genomes_without_taxon_refs = list()


    logger.info("Looking up taxonomy")
    if taxon_reference is None:
        #Get the taxon_lookup_object
        taxon_lookup = ws_client.get_objects( [{'workspace':taxon_wsname,
                                                'name':taxon_lookup_obj_name}])
        if ((organism is not None) and (organism[0:3] in taxon_lookup[0]['data']['taxon_lookup'])):
            if organism in taxon_lookup[0]['data']['taxon_lookup'][organism[0:3]]:
                tax_id = taxon_lookup[0]['data']['taxon_lookup'][organism[0:3]][organism] 
                taxon_object_name = "%s_taxon" % (str(tax_id))
            else:
                genomes_without_taxon_refs.append(organism)
                taxon_object_name = "unknown_taxon"
                genome['notes'] = "Unable to find taxon for this organism : {}.".format(organism )
                genome['scientific_name'] = "Unconfirmed Organism: {}".format(organism )
                report.write("Unable to find taxon for this organism : {}.\n\n".format(organism ))
        else: 
            genomes_without_taxon_refs.append(organism)
            taxon_object_name = "unknown_taxon"
            genome['notes'] = "Unable to find taxon for this organism : {}.".format(organism )
            genome['scientific_name'] = "Unconfirmed Organism: {}".format(organism )
            report.write("Unable to find taxon for this organism : {}.\n\n".format(organism ))
        del taxon_lookup

        try: 
            logger.info("attempting to link to " + taxon_wsname + '/' + taxon_object_name)
            taxon_info = ws_client.get_objects([{"workspace": taxon_wsname, 
                                                 "name": taxon_object_name}])
        except Exception, e: 
            raise Exception("The taxon " + taxon_object_name + " from workspace " + str(taxon_workspace_id) + " does not exist. Error was:" + str(e))
    
        taxon_id = "%s/%s/%s" % (taxon_info[0]["info"][6], taxon_info[0]["info"][0], taxon_info[0]["info"][4]) 
        if not genetic_code_supplied:
            genetic_code = taxon_info[0]["data"]["genetic_code"]
        elif genetic_code != taxon_info[0]["data"]["genetic_code"]:
            #Supplied genetic code differs from taxon genetic code.  Add warning to genome notes
            temp_notes = ""
            if "notes" in genome:
                temp_notes = "{} ".format(genome["notes"])
            genome["notes"] = "{}The supplied genetic code of {} differs from the taxon genetic code of {}. The supplied genetic code is being used.".format(temp_notes,genetic_code, taxon_info[0]["data"]["genetic_code"])
            report.write("The supplied genetic code of {} differs from the taxon genetic code of {}. The supplied genetic code is being used.\n\n".format(genetic_code, taxon_info[0]["data"]["genetic_code"]))
        else:
            temp_notes = ""
            if "notes" in genome:
                temp_notes = "{} ".format(genome["notes"])
            genome["notes"] = "{}The genetic code of {} was supplied by the user.".format(temp_notes,genetic_code, taxon_info[0]["data"]["genetic_code"])
            report.write("The genetic code of {} was supplied by the user.\n\n".format(genetic_code, taxon_info[0]["data"]["genetic_code"]))

        genome['genetic_code'] = genetic_code
#            print "Found name : " + taxon_object_name + " id: " + taxon_id
#            print "TAXON OBJECT TYPE : " + taxon_info[0]["info"][2]
        if not taxon_info[0]["info"][2].startswith("KBaseGenomeAnnotations.Taxon"):
            raise Exception("The object retrieved for the taxon object is not actually a taxon object.  It is {}".format(taxon_info[0]["info"][2]))
        if 'scientific_name' not in genome:
            genome['scientific_name'] = taxon_info[0]['data']['scientific_name']
        genome['domain'] = taxon_info[0]['data']['domain']

    else:

        logger.info("Exact reference provided, using:" + str(taxon_reference))
        try: 
            taxon_info = ws_client.get_objects([{"ref": taxon_reference}])
            taxon_id = "%s/%s/%s" % (taxon_info[0]["info"][6], taxon_info[0]["info"][0], taxon_info[0]["info"][4]) 
        except Exception, e:
            raise Exception("The taxon reference " + str(taxon_reference) + " does not correspond to a workspace object. Error was:" + str(e))
    
        print "TAXON OBJECT TYPE : " + taxon_info[0]["info"][2] 
        if not taxon_info[0]["info"][2].startswith("KBaseGenomeAnnotations.Taxon"):
            raise Exception("The object retrieved for the taxon object is not actually a taxon object.  It is " + taxon_info[0]["info"][2])
        if not genetic_code_supplied:
            genetic_code = taxon_info[0]["data"]["genetic_code"]
        elif genetic_code != taxon_info[0]["data"]["genetic_code"]:
            #Supplied genetic code differs from taxon genetic code.  Add warning to genome notes
            temp_notes = ""
            if "notes" in genome:
                temp_notes = "{} ".format(genome["notes"])
            genome['notes'] ="{}  The supplied genetic code of {} differs from the taxon genetic code of {}. The supplied genetic code is being used.".format(temp_notes,genetic_code, 
                                                                                                                                                           taxon_info[0]["data"]["genetic_code"])
        report.write("The supplied genetic code of {} differs from the taxon genetic code of {}. The supplied genetic code is being used.\n\n".format(genetic_code, 
                                                                                                                                                              taxon_info[0]["data"]["genetic_code"]))
        genome['genetic_code'] = genetic_code
        genome['scientific_name'] = taxon_info[0]['data']['scientific_name']
        genome['domain'] = taxon_info[0]['data']['domain']

    genome['taxonomy'] = taxon_info[0]["data"]["scientific_lineage"]

#EARLY BAILOUT FOR TESTING
#    sys.exit(1)


#import re  
#s = "$abcdef::1234'ghi^_+jkl_mnop" 
#s = re.sub(r'[\W_]+', '_', s) 
 
#    core_scientific_name = re.sub(r'[\W_]+', '_', display_sc_name)
    core_scientific_name = re.sub(r'[\W_]+', '_', genome['scientific_name'])

    #CORE OBJECT NAME WILL BE EITHER PASSED IN GENERATED (TaxID_Source)
    #Fasta file name format is taxID_source_timestamp
    time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
    if core_genome_name is None:
        if source is None:
            source_name = "unknown_source"
        else:
            source_name = source
        if tax_id == 0:
            core_genome_name = "%s_%s" % (source_name,time_string) 
            fasta_file_name = "unknown_%s_%s.fa" % (source_name,time_string) 
        else:
            core_genome_name = "%s_%s" % (core_scientific_name,source_name) 
            fasta_file_name = "%s_%s.fa" % (core_scientific_name,time_string) 
    else:
        fasta_file_name = "%s_%s.fa" % (core_genome_name,time_string) 
        if source is None:
            source_name = "unknown_source"
        else:
            source_name = source

    print "Core Genome Name :"+ core_genome_name + ":"
    print "FASTA FILE Name :"+ fasta_file_name + ":"

    now_date = datetime.datetime.now()
        
    #Parse LOCUS line from each file and grab that meta data (also establish order of the contigs)
    locus_name_order = list() #for knowing order of the genbank files/contigs
    genbank_metadata_objects = dict() #the data structure for holding the top level metadata information of each genbank file
    contig_information_dict = dict() #the data structure for holding the top level metadata information of each genbank file for the stuff needed for making the assembly.

    #HAD TO ADD "CON" as a possible division set.  Even though that does not exist according to this documentation:
    #http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#GenBankDivisionB
    #Found this http://www.ncbi.nlm.nih.gov/Web/Newsltr/Fall99/contig.html , oddly suggests no sequence should be associated with this.
    genbank_division_set = {'PRI','ROD','MAM','VRT','INV','PLN','BCT','VRL','PHG','SYN','UNA','EST','PAT','STS','GSS','HTG','HTC','ENV','CON'}

    #Make the Fasta file for the sequences to be written to
    if not os.path.exists("temp_fasta_file_dir"):
        os.makedirs("temp_fasta_file_dir")
    fasta_file_name = "temp_fasta_file_dir/" +fasta_file_name
    fasta_file_handle = open(fasta_file_name, 'w')
    
    min_date = None
    max_date = None
    genbank_time_string = None
    genome_publication_dict = dict()

#    #Create a SQLLite database and connection.
#    if make_sql_in_memory:
#        sql_conn = sqlite3.connect(':memory:') 
#    else:
#        db_name = "GenomeAnnotation_{}.db".format(time_string) 
#        sql_conn = sqlite3.connect(db_name) 

#    sql_cursor = sql_conn.cursor() 

    #Create a protein and feature table.
#    sql_cursor.execute('''CREATE TABLE features (feature_id text, feature_type text, sequence_length integer, feature_data blob)''')
#    sql_cursor.execute('''CREATE INDEX feature_id_idx ON features (feature_id)''')
#    sql_cursor.execute('''CREATE INDEX feature_type_idx ON features (feature_type)''')
#    sql_cursor.execute('''CREATE INDEX seq_len_idx ON features (sequence_length)''')
#    sql_cursor.execute('''PRAGMA synchronous=OFF''') 

    #Feature Data structures
    list_of_features = list()
    locus_tag_index_dict = dict()

#    #Key is the gene tag (ex: gene="NAC001"), the value is a dict with feature type as the key. The value is a list of maps (one for each feature with that gene value).  
#    #The internal map stores all the key value pairs of that feature.
#    features_grouping_dict = dict() 
    
    #Feature_type_id_counter_dict  keeps track of the count of each time a specific id needed to be generated by the uploader and not the file.
    feature_type_id_counter_dict = dict()

    #Key is feature type, value is the number of occurances of this type. Lets me know the feature containers that will need 
    #to be made and a check to insure the counts are accurate.
    feature_type_counts = dict() 

    #key feature id to be used, value 1
    feature_ids = dict()

    #integers used for stripping text 
    complement_len = len("complement(")
    join_len = len("join(")
    order_len = len("order(")

    genome["num_contigs"] = len(genbank_file_boundaries)

    print "NUMBER OF GENBANK RECORDS: " + str(len(genbank_file_boundaries))
    for byte_coordinates in genbank_file_boundaries: 
        genbank_file_handle.seek(byte_coordinates[0]) 
        genbank_record = genbank_file_handle.read(byte_coordinates[1] - byte_coordinates[0]) 
        try:
            annotation_part, sequence_part = genbank_record.rsplit("ORIGIN",1)
        except Exception, e:
            #sequence does not exist.
            fasta_file_handle.close() 
            raise Exception("This Genbank file has at least one record without a sequence.")

        #done with need for variable genbank_record. Freeing up memory
        genbank_record = None
        metadata_part, features_part = annotation_part.rsplit("FEATURES             Location/Qualifiers",1) 

        metadata_lines = metadata_part.split("\n")

        ##########################################
        #METADATA PARSING PORTION
        ##########################################
        for metadata_line in metadata_lines: 
            if metadata_line.startswith("ACCESSION   "): 
                temp = metadata_line[12:]
                accession = temp.split(' ', 1)[0]
                break

        #LOCUS line parsing
        locus_line_info = metadata_lines[0].split()
        genbank_metadata_objects[accession] = dict()
        contig_information_dict[accession] = dict()
        locus_name_order.append(accession)
        genbank_metadata_objects[accession]["number_of_basepairs"] = locus_line_info[2]
        date_text = None
        if ((len(locus_line_info)!= 7) and (len(locus_line_info)!= 8)): 
            fasta_file_handle.close()
            raise Exception("Error the record with the Locus Name of %s does not have a valid Locus line.  It has %s space separated elements when 6 to 8 are expected (typically 8)." % (locus_info_line[1],str(len(locus_line_info))))
        if locus_line_info[4].upper() != 'DNA':
            if (locus_line_info[4].upper() == 'RNA') or (locus_line_info[4].upper() == 'SS-RNA') or (locus_line_info[4].upper() == 'SS-DNA'):
                if not tax_lineage.lower().startswith("viruses") and not tax_lineage.lower().startswith("viroids"):
                    fasta_file_handle.close()
                    raise Exception("Error the record with the Locus Name of %s is RNA, but the organism does not belong to Viruses or Viroids." % (locus_line_info[1]))
            else:
                fasta_file_handle.close()
                raise Exception("Error the record with the Locus Name of %s is not valid as the molecule type of '%s' , is not 'DNA' or 'RNA'.  If it is RNA it must be a virus or a viroid." % (locus_line_info[1],locus_line_info[4]))
        if ((locus_line_info[5] in genbank_division_set) and (len(locus_line_info) == 7)) :
            genbank_metadata_objects[accession]["is_circular"] = "Unknown"
            contig_information_dict[accession]["is_circular"] = "Unknown"
            date_text = locus_line_info[6]
        elif (locus_line_info[6] in genbank_division_set  and (len(locus_line_info) == 8)) :
            date_text = locus_line_info[7]
            if locus_line_info[5] == "circular":
                genbank_metadata_objects[accession]["is_circular"] = "True"
                contig_information_dict[accession]["is_circular"] = "True"
            elif locus_line_info[5] == "linear":
                genbank_metadata_objects[accession]["is_circular"] = "False"
                contig_information_dict[accession]["is_circular"] = "False"
            else:
                genbank_metadata_objects[accession]["is_circular"] = "Unknown"
                contig_information_dict[accession]["is_circular"] = "Unknown"
        else:
            date_text = locus_line_info[5]

        try:
            record_time = datetime.datetime.strptime(date_text, '%d-%b-%Y')
            if min_date == None:
                min_date = record_time
            elif record_time < min_date:
                min_date = record_time
            if max_date == None:
                max_date = record_time
            elif record_time > max_date:
                max_date = record_time

        except ValueError:
            fasta_file_handle.close()
            exception_string = "Incorrect date format, should be 'DD-MON-YYYY' , attempting to parse the following as a date: %s , the locus line elements: %s " % (date_text, ":".join(locus_line_info))
#            raise ValueError("Incorrect date format, should be 'DD-MON-YYYY' , attempting to parse the following as a date:" + date_text)
            raise ValueError(exception_string)

        genbank_metadata_objects[accession]["external_source_origination_date"] = date_text

        num_metadata_lines = len(metadata_lines)
        metadata_line_counter = 0

        for metadata_line in metadata_lines:
            if metadata_line.startswith("DEFINITION  "):
                definition = metadata_line[12:]
                definition_loop_counter = 1
                if ((metadata_line_counter + definition_loop_counter)<= num_metadata_lines):
                    next_line = metadata_lines[metadata_line_counter + definition_loop_counter]
                    while (next_line.startswith("            ")) and ((metadata_line_counter + definition_loop_counter)<= num_metadata_lines) :
                        definition = "%s %s" % (definition,next_line[12:])
                        definition_loop_counter += 1
                        if ((metadata_line_counter + definition_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + definition_loop_counter]
                        else:
                            break
                genbank_metadata_objects[accession]["definition"] = definition 
                contig_information_dict[accession]["definition"] = definition 
            elif metadata_line.startswith("  ORGANISM  "): 
                organism = metadata_line[12:] 
                if organism not in organism_dict:
                    fasta_file_handle.close()
                    raise ValueError("There is more than one organism represented in these Genbank files, they do not represent single genome. First record's organism is %s , but %s was also found" 
                                     % (str(organism_dict.keys()),organism)) 
            elif metadata_line.startswith("COMMENT     "):
                comment = metadata_line[12:] 
                comment_loop_counter = 1 
                if ((metadata_line_counter + comment_loop_counter)<= num_metadata_lines):
                    next_line = metadata_lines[metadata_line_counter + comment_loop_counter] 
                    while (next_line.startswith("            ")) : 
                        comment = "%s %s" % (comment,next_line[12:]) 
                        comment_loop_counter += 1 
                        if ((metadata_line_counter + comment_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + comment_loop_counter]
                        else:
                            break
#                genome_comment = "%s<%s :: %s> " % (genome_comment,accession,comment)
#                genome_comment_io.write("<%s :: %s> " % (accession,comment))
            elif metadata_line.startswith("REFERENCE   "):
                #PUBLICATION SECTION (long)
                authors = ''
                title = ''
                journal = ''
                pubmed = ''
                consortium = ''
                publication_key = metadata_line

                reference_loop_counter = 1
                if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines): 
                    next_line = metadata_lines[metadata_line_counter + reference_loop_counter] 
                # while (next_line and re.match(r'\s', next_line) and not nextline[0].isalpha()):
                while (next_line and re.match(r'\s', next_line)):
                    publication_key += next_line
                    if next_line.startswith("  AUTHORS   "):
                        authors = next_line[12:] 
                        reference_loop_counter += 1
                        if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + reference_loop_counter] 
                        else:
                            break
                        while (next_line.startswith("            ")) :     
                            authors = "%s %s" % (authors,next_line[12:]) 
                            reference_loop_counter += 1
                            if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines): 
                                next_line = metadata_lines[metadata_line_counter + reference_loop_counter] 
                            else: 
                                break 
                    elif next_line.startswith("  TITLE     "):
                        title = next_line[12:]
                        reference_loop_counter += 1
                        if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                        else:
                            break
                        while (next_line.startswith("            ")) :
                            title = "%s %s" % (title,next_line[12:])
                            reference_loop_counter += 1
                            if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                                next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                            else:
                                break
                    elif next_line.startswith("  JOURNAL   "):
                        journal = next_line[12:]
                        reference_loop_counter += 1
                        if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                        else:
                            break
                        while (next_line.startswith("            ")) :
                            journal = "%s %s" % (journal,next_line[12:])
                            reference_loop_counter += 1
                            if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                                next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                            else:
                                break
                    elif next_line.startswith("   PUBMED   "): 
                        pubmed = next_line[12:] 
                        reference_loop_counter += 1
                        if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                        else:
                            break
                        while (next_line.startswith("            ")) : 
                            pubmed = "%s %s" % (journal,next_line[12:]) 
                            reference_loop_counter += 1
                            if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines): 
                                next_line = metadata_lines[metadata_line_counter + reference_loop_counter] 
                            else: 
                                break 
                    elif next_line.startswith("  CONSRTM   "):
                        consortium = next_line[12:]
                        reference_loop_counter += 1
                        if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines): 
                            next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                        else:
                            break 
                        while (next_line.startswith("            ")) : 
                            consortium = "%s %s" % (journal,next_line[12:]) 
                            reference_loop_counter += 1
                            if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                                next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                            else:
                                break
                    else:
                        reference_loop_counter += 1
                        if ((metadata_line_counter + reference_loop_counter)<= num_metadata_lines):
                            next_line = metadata_lines[metadata_line_counter + reference_loop_counter]
                        else:
                            break
                #Done grabbing reference lines, time to build the reference object.

                pubmed_link = ''
                publication_source = ''
                publication_date = ''
                if pubmed != '':
                    publication_source = "PubMed"
                elif consortium != '':
                    publication_source = consortium
                try:
                    pubmed = int(pubmed)
                except ValueError:
                    pubmed = 0
                if pubmed != 0:
                    pubmed_link = "http://www.ncbi.nlm.nih.gov/pubmed/%s" % str(pubmed)
                if journal != '':
                    potential_date_regex = r'(?<=\().+?(?=\))'
                    potential_dates = re.findall(potential_date_regex, journal)
                    
                    for potential_date in reversed(potential_dates):                        
                        try:
                            record_time = datetime.datetime.strptime(potential_date, '%d-%b-%Y')
                            if now_date > record_time:
                                publication_date = potential_date
                                break
                        except ValueError:
                            try:
                                record_time = datetime.datetime.strptime(potential_date, '%b-%Y')
                                if now_date > record_time:
                                    publication_date = potential_date
                                    break       
                            except ValueError:
                                try:
                                    record_time = datetime.datetime.strptime(potential_date, '%Y')
                                    if now_date > record_time:
                                        publication_date = potential_date
                                        break
                                except ValueError:
                                    next
                publication = [pubmed,publication_source,title,pubmed_link,publication_date,authors,journal]
                genome_publication_dict[publication_key] = publication
                #END OF PUBLICATION SECTION

            metadata_line_counter += 1

        if len(genome_publication_dict) > 0 :
            genome["publications"] = genome_publication_dict.values() 

        ##################################################################################################
        #MAKE SEQUENCE PART INTO CONTIG WITH NO INTERVENING SPACES OR NUMBERS
        ##################################################################################################
        sequence_part = re.sub('[0-9]+', '', sequence_part)
        sequence_part = re.sub('\s+','',sequence_part)
        sequence_part = sequence_part.replace("?","")

        contig_length = len(sequence_part)
        if contig_length == 0:
            fasta_file_handle.close() 
            raise Exception("The genbank record %s does not have any sequence associated with it." % (accession))
            

        ##################################################################################################
        #FEATURE ANNOTATION PORTION - Build up datastructures to be able to build feature containers.
        ##################################################################################################
        #print "GOT TO FEATURE PORTION"
        features_lines = features_part.split("\n") 

        num_feature_lines = len(features_lines)
        features_list = list()

        #break up the features section into individual features.
        for feature_line_counter in range(0,(num_feature_lines)):
            feature_line = features_lines[feature_line_counter]
            if ((not feature_line.startswith("                     ")) and (feature_line.startswith("     ")) and (feature_line[5:7].strip() != "")):
                #Means a new feature:
                #
                current_feature_string = feature_line
                while ((feature_line_counter + 1) < num_feature_lines) and (features_lines[(feature_line_counter + 1)].startswith("                     ")): 
                    feature_line_counter += 1 
                    feature_line = features_lines[feature_line_counter]
                    current_feature_string += " %s" % (feature_line)

                features_list.append(current_feature_string)

            elif ((feature_line_counter + 1) < num_feature_lines): 
                feature_line_counter += 1 
                feature_line = features_lines[feature_line_counter]
        
        #Go through each feature and determine key value pairs, properties and importantly the id to use to group for interfeature_relationships.
        for feature_text in features_list:
#            if("gene" in feature_text):
#                print feature_text
#                sys.exit()
            feature_object = dict()
            #split the feature into the key value pairs. "/" denotes start of a new key value pair.
            feature_key_value_pairs_list = feature_text.split("                     /")
            feature_header = feature_key_value_pairs_list.pop(0)
            if len(feature_header[:5].strip()) != 0:
                continue
            coordinates_info = feature_header[21:] 
            feature_type = feature_header[:21] 
            feature_type = feature_type.strip().replace(" ","_")
            if feature_type not in ['CDS','gene']:
                #skip non core feature types. We currently decided to not include mRNA
                continue
            feature_object["type"] = feature_type

            quality_warnings = list() #list of warnings about the feature. Can do more with this at a later time.
            feature_keys_present_dict = dict() #dict of keys present in the feature

            #Get feature key value pairs
            for feature_key_value_pair in feature_key_value_pairs_list: 
                #the key value pair removing unnecessary white space (including new lines as these often span multiple lines)
                temp_string = re.sub( '\s+', ' ', feature_key_value_pair ).strip() 
                try: 
                    key, value = temp_string.split('=', 1) 
                except Exception, e: 
                    #Does not follow key value pair structure.  This unexpected. Skipping.
                    key = temp_string
                    value = ""

                value = re.sub(r'^"|"$', '', value.strip())
                feature_keys_present_dict[key.strip()] = 1

            coordinates_info = re.sub( '\s+', '', coordinates_info ).strip()
            original_coordinates = coordinates_info
            coordinates_list = list()
            apply_complement_to_all = False
            need_to_reverse_locations = False
            has_odd_coordinates = False
            can_not_process_feature = False
            if coordinates_info.startswith("complement") and coordinates_info.endswith(")"): 
                apply_complement_to_all = True
                need_to_reverse_locations = True
                coordinates_info = coordinates_info[complement_len:-1]
            if coordinates_info.startswith("join") and coordinates_info.endswith(")"):
                coordinates_info = coordinates_info[join_len:-1]
            if coordinates_info.startswith("order") and coordinates_info.endswith(")"):
                coordinates_info = coordinates_info[order_len:-1]
                has_odd_coordinates = True
                temp_warning = "Feature with the text %s has the rare 'order' coordinate. The sequence was joined together because KBase does not allow for a non contiguous resulting sequence with multiple locations for a feature.\n\n" % (feature_text)
#                quality_warnings.append(temp_warning)
                report.write(temp_warning)
                #annotation_metadata_warnings.append(temp_warning)
#                sql_cursor.execute("insert into annotation_metadata_warnings values(:warning)",(temp_warning,))
            coordinates_list = coordinates_info.split(",")
            last_coordinate = 0
            dna_sequence_length = 0
            dna_sequence = ''
            locations = list()#list of location objects
            for coordinates in coordinates_list:
                apply_complement_to_current = False
                if coordinates.startswith("complement") and coordinates.endswith(")"): 
                    apply_complement_to_current = True 
                    coordinates = coordinates[complement_len:-1]
                #Look for and handle odd coordinates
                if (("<" in coordinates) or (">" in coordinates)):
                    has_odd_coordinates = True
                    temp_warning = "Feature with the text %s has a '<' or a '>' in the coordinates.  This means the feature starts or ends beyond the known sequence.\n\n" % (feature_text)
                    #quality_warnings.append(temp_warning)
                    report.write(temp_warning)
                    #annotation_metadata_warnings.append(temp_warning)
#                    sql_cursor.execute("insert into annotation_metadata_warnings values(:warning)",(temp_warning,))
                    coordinates= re.sub('<', '', coordinates)
                    coordinates= re.sub('>', '', coordinates)

                period_count = coordinates.count('.')
                if ((period_count == 2) and (".." in coordinates)):
                    start_pos, end_pos = coordinates.split('..', 1)                    
                elif period_count == 0:
                    start_pos = coordinates
                    end_pos = coordinates
                elif period_count == 1:
                    start_pos, end_pos = coordinates.split('.', 1) 
                    has_odd_coordinates = True
                    temp_warning = "Feature with the text %s has a single period in the original coordinate this indicates that the exact location is unknown but that it is one of the bases between bases %s and %s, inclusive.  Note the entire sequence range has been put into this feature.\n\n" % (feature_text, str(start_pos),str(end_pos))
                    #quality_warnings.append(temp_warning)
                    report.write(temp_warning)
                    #annotation_metadata_warnings.append(temp_warning)
#                    sql_cursor.execute("insert into annotation_metadata_warnings values(:warning)",(temp_warning,))
                elif period_count > 2 :
                    can_not_process_feature = True
                else:
                    can_not_process_feature = True
                if "^" in coordinates:
                    start_pos, end_pos = coordinates.split('^', 1) 
                    has_odd_coordinates = True
                    temp_warning = "Feature with the text %s is between bases.  It points to a site between bases %s and %s, inclusive.  Note the entire sequence range has been put into this feature.\n\n" % (feature_text, str(start_pos),str(end_pos))
                    #quality_warnings.append(temp_warning)
                    report.write(temp_warning)
                    #annotation_metadata_warnings.append(temp_warning)       
#                    sql_cursor.execute("insert into annotation_metadata_warnings values(:warning)",(temp_warning,))

                if not can_not_process_feature:
                    if (represents_int(start_pos) and represents_int(end_pos)):
                        if int(start_pos) > int(end_pos):
                            fasta_file_handle.close() 
                            print "FEATURE TEXT: " + feature_text
                            raise Exception("The genbank record %s has coordinates that are out of order. Start coordinate %s is bigger than End coordinate %s. Should be ascending order." % (accession, str(start_pos), str(end_pos)))

#CANT COUNT ON THEM BEING IN ASCENDING POSITIONAL ORDER
#                    if (int(start_pos) < last_coordinate or int(end_pos) < last_coordinate) and ("trans_splicing" not in feature_keys_present_dict) :
#                        fasta_file_handle.close()
#                        raise Exception("The genbank record %s has coordinates that are out of order. Start coordinate %s and/or End coordinate %s is larger than the previous coordinate %s within this feature. Should be ascending order since this is not a trans_splicing feature." % (accession, str(start_pos), str(end_pos),str(last_coordinate)))

                        if (int(start_pos) > contig_length) or (int(end_pos) > contig_length):
                            fasta_file_handle.close() 
                            raise Exception("The genbank record %s has coordinates (start: %s , end: %s) that are longer than the sequence length %s." % \
                                            (accession,str(start_pos), int(end_pos),str(contig_length)))

                        segment_length = (int(end_pos) - int(start_pos)) + 1
                        dna_sequence_length += segment_length
                        temp_sequence = sequence_part[(int(start_pos)-1):int(end_pos)] 
                        strand = "+"
                        location_start = int(start_pos)
                        if apply_complement_to_current or apply_complement_to_all: 
                            my_dna = Seq(temp_sequence, IUPAC.ambiguous_dna)
                            my_dna = my_dna.reverse_complement()
                            temp_sequence = str(my_dna).upper()      
                            strand = "-"
                            location_start = location_start + (segment_length - 1)
                        if apply_complement_to_all:
                            dna_sequence =  temp_sequence + dna_sequence 
                        else:
                            dna_sequence +=  temp_sequence 

                        locations.append([accession,location_start,strand,segment_length]) 
                    else:
                        #no valid coordinates
                        print "Feature text : {} :".format(feature_text)
                        fasta_file_handle.close() 
                        raise Exception("The genbank record %s contains coordinates that are not valid number(s).  Feature text is : %s" % (accession,feature_text)) 

                    last_coordinate = int(end_pos)

            if has_odd_coordinates:
                    quality_warnings.insert(0,"Note this feature contains some atypical coordinates, see the rest of the warnings for details : %s" % (original_coordinates))
            if can_not_process_feature: 
                #skip source feature types.
                continue
            
            dna_sequence = dna_sequence.upper()

            if len(locations) > 0:
                if need_to_reverse_locations and (len(locations) > 1):
                    locations.reverse()
            feature_object["location"]=locations

            feature_object["dna_sequence_length"] = dna_sequence_length
            feature_object["dna_sequence"] = dna_sequence
            try:
                feature_object["md5"] = hashlib.md5(dna_sequence).hexdigest() 
            except Exception, e:
#                print "THE FEATURE TEXT IS : %s" % (feature_text)
#                print "THE FEATURE SEQUENCE IS : %s : " % (dna_sequence)
#                print "Help %s" % help(dna_sequence)
                raise Exception(e)

            #Need to determine id for the feature : order selected by gene, then locus.
            alias_dict = dict() #contains locus_tag, gene, gene_synonym, dbxref, then value is 1 (old way value is a list of sources).
            inference = ""
            notes = ""
            additional_properties = dict()
            feature_specific_id = None
            feature_id = None
            product = None
            EC_number = None
            pseudo_non_gene = False
            has_protein_id = False
            ontology_terms = dict()
            ontology_terms_not_found = dict() #Term and then count of occurences

            for feature_key_value_pair in feature_key_value_pairs_list:
                #the key value pair removing unnecessary white space (including new lines as these often span multiple lines)
                temp_string = re.sub( '\s+', ' ', feature_key_value_pair ).strip()

                try: 
                    key, value = temp_string.split('=', 1) 
                except Exception, e: 
                    #Does not follow key value pair structure.  This unexpected. Skipping.
                    if temp_string == "pseudo":
                        if feature_type == "gene":
                            feature_object["type"] = "pseudogene"
                        else:
                            pseudo_non_gene = True
                    elif temp_string != "trans_splicing":
                        temp_warning = "%s has the following feature property does not follow the expected key=value format : %s" % (feature_id, temp_string) 
                        quality_warnings.append(temp_warning)
                        #annotation_metadata_warnings.append(temp_warning)
#                        sql_cursor.execute("insert into annotation_metadata_warnings values(:warning)",(temp_warning,))       
                    key = temp_string 
                    value = "" 

                key = key.strip()
                value = re.sub(r'^"|"$', '', value.strip())

                if key == "gene":
                    feature_object["gene"] = value 
                    alias_dict[value]=1 
                    if source.upper() == "ENSEMBL" and feature_type == "gene":
                        if value in feature_ids:
                            raise Exception("More than one feature has the specific feature id of {}.  All feature ids need to be unique.".format(value))
                        else:
                            feature_id = value
                            feature_ids[value] = 1
#Kept lines, for dealing with aliases if keeping track of sources/source field
#                    if value in alias_dict and ("Genbank Gene" not in alias_dict[value]) :
#                        alias_dict[value].append("Genbank Gene")
#                    else:
#                        alias_dict[value]=["Genbank Gene"] 
                elif key == "locus_tag":
                    feature_object["locus_tag"] = value 
                    alias_dict[value]=1 
                    if source.upper() != "ENSEMBL" and ( feature_type == "gene" or feature_type == "CDS" ):
                        if value in feature_ids:
                            raise Exception("More than one feature has the specific feature id of {}.  All feature ids need to be unique.".format(value))
                        else:
                            feature_id = value
#                    if feature_type == "gene":
#                        feature_object["feature_specific_id"] = value
                elif key == "old_locus_tag" or key == "standard_name":
                    alias_dict[value]=1 
                elif key == "gene_synonym":
                    synonyms = value.split(';') 
                    for i in synonyms:
                        i = i.strip()
                        alias_dict[i]=1 
                elif (key == "transcript_id"):
#                    if feature_type == "mRNA":
#                        feature_object["feature_specific_id"] = value 
                    alias_dict[value]=1 
                elif (key == "protein_id"):
#                    if feature_type == "CDS":
#                        feature_object["feature_specific_id"] = value
#                    if feature_type == "CDS":
#                        if value in feature_ids:
#                            raise Exception("More than one feature has the specific feature id of {}.  All feature ids need to be unique.".format(value))
#                        else:
#                            feature_id = value
#                            feature_ids[value] = 1 
                    alias_dict[value]=1 
                    has_protein_id = True

                elif (key == "db_xref"):
                    try:
                        db_xref_source, db_xref_value = value.strip().split(':',1)
                        db_xref_value = db_xref_value.strip()
                        db_xref_source = db_xref_source.strip()
                        if db_xref_source.upper() == "GO" or db_xref_source.upper() == "PO":
                            if exclude_ontologies == 0:
                                ontology_id=value.strip()
                                ontology_source = db_xref_source.upper()
                                if ontology_source == "GO":
                                    ontology_ref = "KBaseOntology/gene_ontology"
                                elif ontology_source == "PO":
                                    ontology_ref = "KBaseOntology/plant_ontology"
                                if ontology_id not in ontology_sources[ontology_source]:
#                                alias_dict[value]=1 
    #                                print ("Term {} was not found in our ontology database. It is likely a deprecated term.".format(ontology_id))
                                    if ontology_id not in ontology_terms_not_found:
                                        ontology_terms_not_found[ontology_id] = 1
                                    else:
                                        ontology_terms_not_found[ontology_id] = ontology_terms_not_found[ontology_id] + 1
                                else:
                                    if(ontology_source not in ontology_terms):
                                        ontology_terms[ontology_source]=dict()
                                    if( ontology_id not in ontology_terms[ontology_source]):
                                        OntologyEvidence=[{"method":"KBase_Genbank_uploader from db_xref field","timestamp":time_string,"method_version":"1.0"}]
                                        OntologyData={"id":ontology_id,"ontology_ref":ontology_ref,
                                                      "term_name":ontology_sources[ontology_source][ontology_id]["name"],
                                                      "term_lineage":[],"evidence":OntologyEvidence}
                                        ontology_terms[ontology_source][ontology_id]=OntologyData
                        else:
                            alias_dict[value]=1 
                    except Exception, e: 
                        alias_dict[value]=1 
#                        db_xref_source = "Unknown"
#                        db_xref_value = value.strip()
#                    if db_xref_value.strip() in alias_dict: 
#                        if (db_xref_source.strip() not in alias_dict[db_xref_value.strip()]) :
#                            alias_dict[db_xref_value.strip()].append(db_xref_source.strip())
#                    else:
#                        alias_dict[db_xref_value.strip()]=[db_xref_source.strip()]
#                elif (key == "note"):
#                    if notes != "":
#                        notes += ";"
#                    notes += value
                elif (key == "translation"):
                    #
                    # TODO
                    #NOTE THIS IS A PLACE WHERE A QUALITY WARNING CHECK CAN BE DONE, 
                    #see if translation is accurate.(codon start (1,2,3) may need to be used)
                    #
                    value = re.sub('\s+','',value)
                    feature_object["protein_translation"] = value
                    feature_object["protein_translation_length"] = len(value)
                elif ((key == "function") and (value is not None) and (value.strip() == "")) :
                    feature_object["function"] = value
                elif (key == "product"):
                    product = value
#                    additional_properties[key] = value
#                elif (key == "trans_splicing"):
#                    feature_object["trans_splicing"] = 1
#                elif (key == "EC_number") and feature_type == "CDS":
#                    EC_number = value
#                else:
#                    if key in additional_properties:
#                        additional_properties[key] =  "%s::%s" % (additional_properties[key],value)
#                    else:
#                        additional_properties[key] = value


#            if len(additional_properties) > 0:
#                feature_object["additional_properties"] = additional_properties
#            if len(notes) > 0:
#                feature_object["notes"] = notes
#            if len(inference) > 0:
#                feature_object["inference"] = inference
            if len(alias_dict) > 0:
                feature_object["aliases"] = alias_dict.keys()
            if ("function" not in feature_object) and (product is not None):
                feature_object["function"] = product

            if feature_type == 'CDS':
                #GET TRANSLATION OF THE CDS.  IF THE GENBANK FILE DOES NOT HAVE IT.  
                coding_dna = Seq(feature_object["dna_sequence"], generic_dna)
                aa_seq = coding_dna.translate(table=genetic_code, to_stop=True)
                aa_trans_seq = str(aa_seq[0:].upper())

                if "protein_translation" in feature_object:
                    if aa_trans_seq != feature_object["protein_translation"].upper():
                        temp_warning = "%s translated amino acid sequence does not match the supplied amino acid sequence.\n\n" % (feature_id) 
                        report.write(temp_warning)
#                        quality_warnings.append(temp_warning) 
#                        sql_cursor.execute("insert into annotation_metadata_warnings values(:warning)",(temp_warning,)) 
                else:
                    if not pseudo_non_gene:
                        raise Exception("Error: CDS with the text : {} has no protein translation".format(feature_text))
#                    if "dna_sequence" in feature_object:
#                        feature_object["protein_translation"] = aa_trans_seq
#                        feature_object["protein_translation_length"] = len(aa_trans_seq)
            
            if pseudo_non_gene:
                if feature_type == "CDS" and has_protein_id:
                    read.write("Feature text : {} is a CDS with pseudo and protein_id.\n\n".format(feature_text))
                    #don not include this feature.
                continue

            if feature_object["type"] in feature_type_counts:
                feature_type_counts[feature_object["type"]] += 1
            else:
                feature_type_counts[feature_object["type"]] = 1     

            if feature_id is None:
                if generate_ids_if_needed == 1:
                    #MAKE AUTOGENERATED ID
                    #MAKING ALL IDS UNIQUE ACROSS THE GENOME.
                    if feature_type not in feature_type_id_counter_dict:
                        feature_type_id_counter_dict[feature_object["type"]] = 1;
                        feature_id = "%s_%s" % (feature_object["type"],str(1)) 
                    else: 
                        feature_type_id_counter_dict[feature_object["type"]] += 1; 
                        feature_id = "%s_%s" % (feature_type,str(feature_type_id_counter_dict[feature_object["type"]]))
                else:
                    #Throw an error informing user they can set the generate ids if needed checkbox
                    #do specific errors so they know where we look for the ids.
                    raise Exception("There was no feature specific id for {}.  \
For gene type we take the id from the locus tag \
(except for Ensembl, then the gene field)\
For CDS type we take the id from the protein id field. \
NOTE IF YOU WANT THIS STILL UPLOADED GO TO THE \
ADVANCED OPTIONS AND CHECK THE\
\"Generate IDs if needed\" checkbox".format(feature_text))

#            feature_object["quality_warnings"] = quality_warnings

#            ############################################
#            #DETERMINE ID TO USE FOR THE FEATURE OBJECT
#            ############################################
#            if feature_type not in features_type_containers_dict:
#                features_type_containers_dict[feature_type] = dict()
#            feature_id = None

#OLD WAY TRIED TO USE ID FROM THE FEATURE, UNIQUENESS ONLY GUARANTEED WITH FEATURE CONTAINER AND NOT ACROSS THE GENOME ANNOTATION
#            if "feature_specific_id" not in feature_object:
#                if "locus_tag" not in feature_object:
#                    if feature_type not in feature_type_id_counter_dict:
#                        feature_type_id_counter_dict[feature_type] = 1;
#                        feature_id = "%s_%s" % (feature_type,str(1))
#                    else:
#                        feature_type_id_counter_dict[feature_type] += 1;
#                        feature_id = "%s_%s" % (feature_type,str(feature_type_id_counter_dict[feature_type]))
#                else:
#                    feature_id = feature_object["locus_tag"]
#            else:
#                feature_id = feature_object["feature_specific_id"]
#            if feature_id in features_type_containers_dict[feature_type]:
#                #Insure that no duplicate ids exist
#                if feature_type not in feature_type_id_counter_dict:
#                    feature_type_id_counter_dict[feature_type] = 1;
#                    feature_id = "%s_%s" % (feature_type,str(1))
#                else: 
#                    feature_type_id_counter_dict[feature_type] += 1;
#                    feature_id = "%s_%s" % (feature_type,str(feature_type_id_counter_dict[feature_type]))
#END OLD WAY


##NEW WAY:  MAKING ALL IDS UNIQUE ACROSS THE GENOME.
#            if feature_type not in feature_type_id_counter_dict:
#                feature_type_id_counter_dict[feature_type] = 1;
#                feature_id = "%s_%s" % (feature_type,str(1))
#            else: 
#                feature_type_id_counter_dict[feature_type] += 1;
#                feature_id = "%s_%s" % (feature_type,str(feature_type_id_counter_dict[feature_type]))
##END NEW WAY
            if len(ontology_terms) > 0:
                feature_object["ontology_terms"]=ontology_terms
            feature_object["id"] = feature_id

            ########################################
            #CLEAN UP UNWANTED FEATURE KEYS
            #######################################
            if "locus_tag" in feature_object: 
                del feature_object["locus_tag"]
            if "gene" in feature_object: 
                del feature_object["gene"]
            if "feature_specific_id" in feature_object: 
                del feature_object["feature_specific_id"]

#            feature_object["quality_warnings"] = quality_warnings

            #MAKE ENTRY INTO THE FEATURE TABLE
#            pickled_feature = cPickle.dumps(feature_object, cPickle.HIGHEST_PROTOCOL) 
#            sql_cursor.execute("insert into features values(:feature_id, :feature_type , :sequence_length, :feature_data)", 
#                               (feature_id, feature_object["type"], feature_object["dna_sequence_length"], sqlite3.Binary(pickled_feature),))
#            print key,value,feature_id,feature_type
#            print feature_object
   
            #Only add genes to the feature object, but maintain dict of locus_tag indices
            if(feature_object["type"]=="gene"):
                list_of_features.append(feature_object)
                last_index = len(list_of_features)-1
                if(feature_object['id'] in locus_tag_index_dict):
                            raise Exception("More than one gene has this locus tag: {}".format(value))
                locus_tag_index_dict[feature_object['id']]=last_index

            if(feature_object["type"]=="CDS"):
                gene_index = locus_tag_index_dict[feature_object['id']]
                list_of_features[gene_index]['protein_translation']=feature_object['protein_translation']
                list_of_features[gene_index]['protein_translation_length']=feature_object['protein_translation_length']
                list_of_features[gene_index]['function']=feature_object['function']
                aliases_dict=dict()
                if('aliases' in list_of_features[gene_index]):
                   for alias in list_of_features[gene_index]['aliases']:
                       aliases_dict[alias]=1
                if('aliases' in feature_object):
                   for alias in feature_object['aliases']:
                       aliases_dict[alias]=1
                list_of_features[gene_index]['aliases']=aliases_dict.keys()
                if("ontology_terms" in feature_object):
                    list_of_features[gene_index]["ontology_terms"]=feature_object["ontology_terms"]
                                                  
#            if(len(list_of_features)==10):
#                sys.exit()
#        for feature_type in feature_type_counts:
#            print "Feature " + feature_type + "  count: " + str(feature_type_counts[feature_type])


        ##################################################################################################
        #SEQUENCE PARSING PORTION  - Write out to Fasta File
        ##################################################################################################

#        print "The len of sequence part is: " + str(len(sequence_part))
#        print "The number from the record: " + genbank_metadata_objects[accession]["number_of_basepairs"]        
#        print "First 100 of sequence part : " + sequence_part[0:100] 
        fasta_file_handle.write(">{}\n".format(accession))
        #write 80 nucleotides per line
        fasta_file_handle.write(insert_newlines(sequence_part,80))
        
    fasta_file_handle.close()
    if min_date == max_date:
        genbank_time_string = min_date.strftime('%d-%b-%Y').upper()
    else:
        genbank_time_string = "%s to %s" %(min_date.strftime('%d-%b-%Y').upper(), max_date.strftime('%d-%b-%Y').upper())

    ##########################################
    #ASSEMBLY CREATION PORTION  - consume Fasta File
    ##########################################

    logger.info("Calling FASTA to Assembly Uploader")
    assembly_reference = "%s/%s_assembly" % (workspace_name,core_genome_name)
    try:
        fasta_working_dir = str(os.getcwd()) + "/temp_fasta_file_dir"

        print "HANDLE SERVICE URL " + handle_service_url
        assembly.upload_assembly(shock_service_url = shock_service_url,
                                 handle_service_url = handle_service_url,
                                 input_directory = fasta_working_dir,
                                 #                  shock_id = args.shock_id,
                                 #                  handle_id = args.handle_id,
                                 #                  input_mapping = args.input_mapping, 
                                 workspace_name = workspace_name,
                                 workspace_service_url = workspace_service_url,
                                 taxon_reference = taxon_id,
                                 assembly_name = "%s_assembly" % (core_genome_name),
                                 source = source_name,
                                 contig_information_dict = contig_information_dict,
                                 date_string = genbank_time_string,
                                 logger = logger)
        shutil.rmtree(fasta_working_dir)
    except Exception, e: 
        logger.exception(e) 
        sys.exit(1) 

    logger.info("Assembly Uploaded")

#    sys.exit(1)

    #Do size check of the features
#    sql_cursor.execute("select sum(length(feature_data)) from features where feature_type = ?", (feature_type,))
#    sql_cursor.execute("select sum(length(feature_data)) from features")
#    for row in sql_cursor:
#        data_length = row[0]

#    if data_length < 900000000:
        #Size is probably ok Try the save
        #Retrieve the features from the sqllite DB
#        sql_cursor.execute("select feature_id, feature_data from features")

#        for row in sql_cursor: 
#            feature_id = row[0]
#            feature_data = cPickle.loads(str(row[1])) 
#            list_of_features.append(feature_data)

#    else:
        #Features too large
        #raising an error for now.
#        raise Exception("This genome can not be saved due to the resulting object being too large for the workspace")

    #Save genome
    #Then Finally store the GenomeAnnotation.                                                                            

    shock_id = None
    handle_id = None
    if shock_id is None:
        shock_info = script_utils.upload_file_to_shock(logger, shock_service_url, input_file_name, token=token)
        shock_id = shock_info["id"]
        handles = script_utils.getHandles(logger, shock_service_url, handle_service_url, [shock_id], [handle_id], token)   
        handle_id = handles[0]

    genome['genbank_handle_ref'] = handle_id
    # setup provenance
    provenance_action = {"script": __file__, "script_ver": "0.1", "description": "features from upload from %s" % (source_name)}
    genome_annotation_provenance = []
    if provenance is not None:
        genome_annotation_provenance = provenance
    genome_annotation_provenance.append(provenance_action)
    genome_object_name = core_genome_name 
    genome['type'] = type 
    if type == "Reference":
        genome['reference_annotation'] = 1
    else:
        genome['reference_annotation'] = 0
    genome['taxon_ref'] = taxon_id
    genome['original_source_file_name'] = source_file_name
    genome['assembly_ref'] =  assembly_reference 
    genome['id'] = genome_object_name
    genome['source'] = source_name
    temp_source_id = locus_name_order[0]
    if len(locus_name_order) > 1:
        temp_source_id += ' (' + str(len(locus_name_order) - 1) + ' more accessions)'
    genome['source_id'] = temp_source_id
    genome['external_source_origination_date'] = genbank_time_string
    genome['features'] = list_of_features
    if release is not None:
        genome['release'] = release
    if len(ontology_terms_not_found) > 0:
        report.write("\nThere were ontologies in the source file that were not found in the onology database.\n\
These are like to be deprecated terms.\n\
Below is a list of the term and the countof the number of features that contained that term:\n")

        for term in ontology_terms_not_found:
            report.write("{} --- {}\n".format(term,str(ontology_terms_not_found[term])))
        report.write("\n")

#    print "Genome id %s" % (genome['id'])
 
    logger.info("Attempting Genome save for %s" % (genome_object_name))
#    while genome_annotation_not_saved:
#        try:
    genome_annotation_info =  ws_client.save_objects({"workspace":workspace_name,
                                                      "objects":[ { "type":"KBaseGenomes.Genome",
                                                                    "data":genome,
                                                                    "name": genome_object_name,
                                                                    "provenance":genome_annotation_provenance,
                                                                    "meta":usermeta
                                                                }]}) 
#            genome_annotation_not_saved = False 
    logger.info("Genome saved for %s" % (genome_object_name))
#        except biokbase.workspace.client.ServerError as err: 
#            raise 

#    if not make_sql_in_memory:
#        os.remove(db_name) 

    logger.info("Conversions completed.")
    report.write("\n\nGENOME AND ASSEMBLY OBJECTS HAVE BEEN SUCESSFULLY SAVED.")


    output_data_ref = "{}/{}".format(workspace_name,genome_object_name)
    reportObj = {
        'objects_created':[{'ref':output_data_ref, 'description':'Assembled contigs'}],
        'text_message':report.getvalue()
    }
    report_kb = KBaseReport(callback_url)
    report_info = report_kb.create({'report':reportObj, 'workspace_name':workspace_name})
    report.close

    return {
        'genome_info': genome_annotation_info[0],
        'report_name': report_info['name'],
        'report_ref': report_info['ref']
    }


















# called only if script is run from command line
if __name__ == "__main__":
    script_details = script_utils.parse_docs(upload_genome.__doc__)    

    import argparse

    parser = argparse.ArgumentParser(prog=__file__, 
                                     description=script_details["Description"],
                                     epilog=script_details["Authors"])
                                     
    parser.add_argument('--shock_service_url', 
                        help=script_details["Args"]["shock_service_url"],
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--handle_service_url', 
                        action='store', type=str, nargs='?', default=None, required=True)
    parser.add_argument('--workspace_name', nargs='?', help='workspace name to populate', required=True)
    parser.add_argument('--taxon_wsname', nargs='?', help='workspace name with taxon in it, assumes the same workspace_service_url', required=False, default='ReferenceTaxons')
#    parser.add_argument('--taxon_names_file', nargs='?', help='file with scientific name to taxon id mapping information in it.', required=False, default="/homes/oakland/jkbaumohl/Genome_Spec_files/Taxonomy/names.dmp")
    parser.add_argument('--taxon_reference', nargs='?', help='ONLY NEEDED IF PERSON IS DOING A CUSTOM TAXON NOT REPRESENTED IN THE NCBI TAXONOMY TREE', required=False)
    parser.add_argument('--workspace_service_url', action='store', type=str, nargs='?', required=True) 

    parser.add_argument('--object_name', 
                        help="genbank file", 
                        nargs='?', required=False)
    parser.add_argument('--source', 
                        help="data source : examples Refseq, Genbank, Pythozyme, Gramene, etc", 
                        nargs='?', required=False, default="Genbank") 
    parser.add_argument('--type', 
                        help="data source : examples Reference, Representative, User Upload", 
                        nargs='?', required=False, default="User upload") 
    parser.add_argument('--release', 
                        help="Release or version of the data.  Example Ensembl release 30", 
                        nargs='?', required=False) 
    parser.add_argument('--genetic_code', 
                        help="genetic code for the genome, normally determined by taxon information. Will override taxon supplied genetic code if supplied. Defaults to 1", 
                        nargs='?', type=int, required=False)
    parser.add_argument('--generate_ids_if_needed', 
                        help="If the fields used for ID determination are not present the uploader will fail by default. If generate_ids_id_needed is 1 then it will generate IDs (Feature_AutoincrementNumber format)", 
                        nargs='?', type=int, required=False)
    parser.add_argument('--exclude_ontologies', 
                        help="Some larger genomes may not fit in the 1 GB limit, one way to increase likelihood they will fit is to exclude the ontologies.", 
                        nargs='?', type=int, required=False)
    parser.add_argument('--input_directory', 
                        help="directory the genbank file is in", 
                        action='store', type=str, nargs='?', required=True)

    args, unknown = parser.parse_known_args()

    logger = script_utils.stderrlogger(__file__)

    logger.debug(args)

    try:
        obj_name = upload_genome(shock_service_url = args.shock_service_url,
                                 handle_service_url = args.handle_service_url, 
                                 input_directory = args.input_directory, 
                                 workspace_name = args.workspace_name,
                                 workspace_service_url = args.workspace_service_url,
                                 taxon_wsname = args.taxon_wsname,
                                 taxon_reference = args.taxon_reference,
                                 taxon_lookup_obj_name='taxon_lookup',
                                 core_genome_name = args.object_name,
                                 source = args.source,
                                 release = args.release,
                                 type = args.type,
                                 genetic_code = args.genetic_code,
                                 generate_ids_if_needed = args.generate_ids_if_needed,
                                 logger = logger)
    except Exception, e:
        logger.exception(e)
        sys.exit(1)

    sys.exit(0)


