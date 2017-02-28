#!/usr/bin/env python
from __future__ import print_function

#GFF3 format
#http://www.sequenceontology.org/gff3.shtml
#http://gmod.org/wiki/GFF3

# Standard imports
import sys,os,time,datetime,re
import gzip,shutil,subprocess,copy
import itertools,hashlib,logging,collections
from pprint import pprint

# Defining pythonic stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# 3rd party imports
import simplejson
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

codon_table = CodonTable.ambiguous_generic_by_name["Standard"]

# KBase imports
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil

def convert_ftr_object(old_ftr,contig):
    new_ftr = dict()
    new_ftr["id"] = old_ftr["ID"]

    #GFF use 1-based integers
    substr_start = old_ftr["start"]-1 
    substr_end = old_ftr["end"]
    if(old_ftr["strand"] == "-"):
        substr_start = old_ftr["end"]-1
        substr_end = old_ftr["start"]

    dna_sequence = Seq(contig[old_ftr["start"]-1:old_ftr["end"]], IUPAC.ambiguous_dna)

    #reverse complement
    if(old_ftr["strand"] == "-"):
        dna_sequence = dna_sequence.reverse_complement()
        old_ftr["start"]=old_ftr["end"]

    new_ftr["dna_sequence"]=str(dna_sequence).upper()
    new_ftr["dna_sequence_length"]=len(dna_sequence)
    new_ftr["md5"]=hashlib.md5(str(dna_sequence)).hexdigest()
    new_ftr["location"] = [[old_ftr["contig"],old_ftr["start"],old_ftr["strand"],len(dna_sequence)]]
    return new_ftr

def upload_genome(shock_service_url=None,
                  handle_service_url=None,
                  workspace_service_url=None,
                  callback_url=None,

                  input_gff_file=None,
                  input_fasta_file=None,
                  
                  workspace_name=None,
                  core_genome_name=None,

                  genome_type=None,
                  scientific_name=None,
                  taxonomy=None,
                  taxon_reference=None,
                  source=None):

    assembly_ref = None
    gff_handle_ref = None
    time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))

    ###########################################
    #Create logger
    ###########################################
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.INFO)
    
    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stderr)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)
    ###########################################
    #End logger creation
    ###########################################

    ##########################################
    #Reading in Fasta file, Code taken from https://www.biostars.org/p/710/
    ##########################################
    logger.info("Reading FASTA file.") 

    assembly = {"contigs":{},"dna_size":0,"gc_content":0,"md5":[],"base_counts":{}}
    contig_seq_start=0

#    input_file_handle = gzip.open(input_fasta_file,'rb')
    input_file_handle = open(input_fasta_file,'rb')
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in itertools.groupby(input_file_handle, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())

        try:
            fasta_header,fasta_description = header.split(' ',1)
        except:
            fasta_header = header
            fasta_description = None

        #Handle record
        seq=seq.upper()

        #Build contig objects for Assembly
        seq_count = dict(collections.Counter(seq))

        #to delete at end, but required for now
        contig_dict = {"sequence":seq}

        Ncount = 0
        if "N" in seq_count:
            Ncount = seq_count["N"]
        contig_dict["Ncount"] = Ncount

        for character in seq_count:
            if character in assembly["base_counts"]:
                assembly["base_counts"][character]+= seq_count[character]
            else:
                assembly["base_counts"][character] = seq_count[character]

        contig_seq_length = len(seq)
        assembly["dna_size"]+=contig_seq_length

        contig_gc_length = seq.count("G")
        contig_gc_length+= seq.count("C")
        contig_dict["gc_content"] = float("{0:.2f}".format(float(contig_gc_length)/float(contig_seq_length)))
        assembly["gc_content"]+=contig_gc_length

        contig_dict["contig_id"] = fasta_header
        contig_dict["name"] = fasta_header
        contig_dict["length"] = contig_seq_length
        contig_dict["md5"] = hashlib.md5(seq).hexdigest()
        assembly["md5"].append(contig_dict["md5"])

        if fasta_description is not None: 
            contig_dict["description"] = fasta_description

        contig_dict["is_circular"] = "Unknown"
        contig_dict["start_position"] = contig_seq_start
        contig_dict["num_bytes"] = sys.getsizeof(contig_dict["sequence"])

        assembly["contigs"][fasta_header]=contig_dict
        
        #used for start of next sequence and total gc_content
        contig_seq_start+=contig_seq_length

    assembly["gc_content"] = float("{0:.2f}".format(float(assembly["gc_content"])/float(contig_seq_start)))
    assembly["md5"] = hashlib.md5(",".join(assembly["md5"])).hexdigest()
    assembly["assembly_id"] = core_genome_name+"_assembly"
    assembly["name"] = scientific_name
    assembly["external_source"] = source
    assembly["external_source_id"] = os.path.basename(input_fasta_file)
    assembly["external_source_origination_date"] = str(os.stat(input_fasta_file).st_ctime)
    assembly["num_contigs"] = len(assembly["contigs"].keys())
    assembly["type"] = "Unknown"
    assembly["notes"] = "Note MD5s are generated from uppercasing the sequences" 

    if taxon_reference is not None:
        assembly["taxon_ref"] = taxon_reference

    logger.info("Reading GFF file.")

    header = list()
    feature_list = dict()
    original_CDS_count=dict()
    original_feature_ids=dict()

#    gff_file_handle = gzip.open(input_gff_file, 'rb')
    gff_file_handle = open(input_gff_file, 'rb')
    current_line = gff_file_handle.readline()
    gff_object = dict()
    while ( current_line != '' ):
        current_line=current_line.strip()
        
        if(current_line.startswith("##") or current_line.startswith("#!")):
            header.append(current_line)
            if('headers' not in gff_object):
                gff_object['headers']=list()
            gff_object['headers'].append(current_line)
        else:
            if('features' not in gff_object):
                gff_object['features']=list()

            contig_id, source_id, feature_type, start, end, score, strand, phase, attributes = current_line.split('\t')
            attributes_dict=dict()
            for attribute in attributes.split(";"):
                if(attribute == "" or "=" not in attribute):
                    continue
                key, value = attribute.split("=",1)
                attributes_dict[key]=value            

            #ID should be transferred from Name or Parent
            old_id=None
            for key in ("ID","PACid","pacid"):
                if(key in attributes_dict):
                    old_id=attributes_dict[key]
                    break
            if(old_id is None):
                eprint("Cannot find unique ID, PACid, or pacid in GFF attributes: "+attributes)
                continue

            if("Name" in attributes_dict):
                attributes_dict["ID"]=attributes_dict["Name"]
            else:
                attributes_dict["ID"]=original_feature_ids[attributes_dict["Parent"]]+"."+feature_type

                #if CDS have to increment
                if(feature_type == "CDS"):
                    if(attributes_dict["ID"] not in original_CDS_count):
                        original_CDS_count[attributes_dict["ID"]]=1
                    else:
                        original_CDS_count[attributes_dict["ID"]]+=1

                    attributes_dict["ID"]+="."+str(original_CDS_count[attributes_dict["ID"]])

            #Update parent
            if("Parent" in attributes_dict):
                attributes_dict["Parent"]=original_feature_ids[attributes_dict["Parent"]]

            original_feature_ids[old_id]=attributes_dict["ID"]

            #recreate line for GFF
            partial_line, attributes = current_line.rsplit('\t',1)
            new_line = partial_line + "\t" + ";".join(key+"="+attributes_dict[key] for key in attributes_dict.keys())
            gff_object['features'].append(new_line)

            if(contig_id not in assembly["contigs"]):
                logger.warn("Missing contig: "+contig_id)

            if(contig_id not in feature_list):
                feature_list[contig_id]=list()

            feature = {'type':feature_type,'start':int(start),'end':int(end),'score':score,'strand':strand,'phase':phase}
            for attribute in attributes.split(";"):
                if(attribute == "" or "=" not in attribute):
                    continue
                key, value = attribute.split("=",1)
                feature[key]=value

            #Append contig identifier
            feature["contig"]=contig_id
            feature_list[contig_id].append(feature)

        current_line = gff_file_handle.readline()
    gff_file_handle.close()

    #Writing updated lines to gff_file_handle
    input_gff_file = input_gff_file.replace("gene","edited_gene")
    gff_file_handle = gzip.open(input_gff_file, 'wb')
    if('headers' in gff_object):
        gff_file_handle.write("\n".join(gff_object["headers"]))
    gff_file_handle.write("\n".join(gff_object["features"]))
    gff_file_handle.close()

    #New code inserted to better handle feature identifiers
    #Start by extracting and group them first
    features_identifiers_dict = dict()
    features_identifiers_list = list()
    features_identifiers_count = dict()
    features_parents_dict = dict()
    features_name_id_dict = dict()
    CDS_count=dict()
    for contig in sorted(feature_list):
        for feature in feature_list[contig]:
            #We're only considering gene, mRNA, and CDS for brevity's sake
            if(feature["type"] not in ("gene","mRNA","CDS")):
                continue

            #gene and mRNA always have name, CDS do not
            if("Name" not in feature):
                feature["Name"]=None

            #Update parent following name/id switch
            if("Parent" in feature and feature["Parent"] in features_name_id_dict):
                feature["Parent"]=features_name_id_dict[feature["Parent"]]

            #ID should be transferred to Name, but need to maintain parent
            if(feature["Name"] is not None):
                features_name_id_dict[feature["ID"]]=feature["Name"]
                feature["ID"]=feature["Name"]
            else:
                feature["ID"]=feature["Parent"]+"."+feature["type"]
                #if CDS have to increment
                if(feature["type"] == "CDS"):
                    if(feature["ID"] not in CDS_count):
                        CDS_count[feature["ID"]]=1
                    else:
                        CDS_count[feature["ID"]]+=1

                    feature["ID"]+="."+str(CDS_count[feature["ID"]])

            #Collect
            if(feature["type"] == "gene"):
                features_identifiers_dict[feature["ID"]]=dict()
            if(feature["type"] == "mRNA"):
                features_identifiers_dict[feature["Parent"]][feature["ID"]]=dict()
                features_parents_dict[feature["ID"]]=feature["Parent"]
            if(feature["type"] == "CDS"):
                features_identifiers_dict[features_parents_dict[feature["Parent"]]][feature["Parent"]][feature["ID"]]=1

            features_identifiers_list.append(feature)
            features_identifiers_count[feature["ID"]]=len(features_identifiers_list)-1

    updated_features_identifiers_dict = dict()
    updated_features_list = list()
    updated_features_identifiers_count = dict()
    updated_features_parents_dict = dict()
    updated_CDS_count = dict()
    for gene in sorted(features_identifiers_dict):

        #retrieve original object
        gene_ftr = features_identifiers_list[features_identifiers_count[gene]]

        #store gene
        updated_features_identifiers_dict[gene_ftr["ID"]]=dict()
        updated_features_list.append(gene_ftr)
        updated_features_identifiers_count[gene_ftr["ID"]]=len(updated_features_list)-1

        for mRNA in sorted(features_identifiers_dict[gene], key=lambda x: features_identifiers_count[x]):
            #retrieve feature
            mRNA_ftr = features_identifiers_list[features_identifiers_count[mRNA]]

            if("PAC" in mRNA[0:3]):
                if("Name" in mRNA_ftr):
                    mRNA_ftr["ID"]=mRNA_ftr["Name"]

            updated_features_identifiers_dict[gene_ftr["ID"]][mRNA_ftr["ID"]]=dict()
            updated_features_parents_dict[mRNA_ftr["ID"]]=mRNA_ftr["Parent"]

            updated_features_list.append(mRNA_ftr)
            updated_features_identifiers_count[mRNA_ftr["ID"]]=len(updated_features_list)-1

            for CDS in sorted(features_identifiers_dict[gene][mRNA],key=lambda x: features_identifiers_count[x]):
                #retrieve feature
                CDS_ftr = features_identifiers_list[features_identifiers_count[CDS]]

                if("PAC" in CDS[0:3]):
                    CDS_ftr["ID"]=mRNA_ftr["ID"]+".CDS"

                    if(CDS_ftr["ID"] not in updated_CDS_count):
                        updated_CDS_count[CDS_ftr["ID"]]=1
                    else:
                        updated_CDS_count[CDS_ftr["ID"]]+=1

                    CDS_ftr["ID"]+="."+str(updated_CDS_count[CDS_ftr["ID"]])
                    CDS_ftr["Parent"]=mRNA_ftr["ID"]

                updated_features_identifiers_dict[gene_ftr["ID"]][mRNA_ftr["ID"]][CDS_ftr["ID"]]=1
                updated_features_parents_dict[CDS_ftr["ID"]]=CDS_ftr["Parent"]

                updated_features_list.append(CDS_ftr)
                updated_features_identifiers_count[CDS_ftr["ID"]]=len(updated_features_list)-1

    genome_features_list = list()
    genome_mrnas_list = list()
    genome_cdss_list = list()
    for gene in sorted(updated_features_identifiers_dict):
        #retrieve updated object
        gene_ftr = updated_features_list[updated_features_identifiers_count[gene]]

        gene_object = convert_ftr_object(gene_ftr,assembly["contigs"][gene_ftr["contig"]]["sequence"])
        gene_object["type"]="gene"

        #New terms, TODO, move to end of gene loop
        gene_object["cdss"]=list()
        gene_object["mrnas"]=list()

        #use function of longest CDS for gene
        longest_protein_length=0
        longest_protein_sequence=""
        for mRNA in sorted(updated_features_identifiers_dict[gene], key=lambda x: updated_features_identifiers_count[x]):
            #retrieve updated object
            mRNA_ftr = updated_features_list[updated_features_identifiers_count[mRNA]]
            
            feature_object = convert_ftr_object(mRNA_ftr,assembly["contigs"][mRNA_ftr["contig"]]["sequence"])
            feature_object['parent_gene']=gene_object['id']
            
            mrna_object=copy.deepcopy(feature_object)
            cds_object=copy.deepcopy(feature_object)

            cds_object['id']=mrna_object['id']+".CDS"
            mrna_object['cds']=cds_object['id']

            cds_object['parent_mrna']=mrna_object['id']

            del mrna_object["dna_sequence"]
            del mrna_object["dna_sequence_length"]

            cds_object["ontology_terms"]=dict()

            gene_object["mrnas"].append(mrna_object["id"])
            gene_object["cdss"].append(cds_object["id"])

            #CDS aggregation needs to be done in order to build protein sequence and list of locations
            CDS_list = sorted(updated_features_identifiers_dict[gene][mRNA], key=lambda x: updated_features_identifiers_count[x])

            dna_sequence = ""
            locations = list()

            #collect phases, and lengths of exons
            #right now, this is only for the purpose of error reporting
            phases = list()
            exons = list()

            for CDS in (CDS_list):
                #retrieve updated partial CDS
                add_ftr = updated_features_list[updated_features_identifiers_count[CDS]]
                phases.append(add_ftr["phase"])

                add_ftr_obj = convert_ftr_object(add_ftr,assembly["contigs"][add_ftr["contig"]]["sequence"])
                exons.append(len(add_ftr_obj["dna_sequence"]))

                #Remove base(s) according to phase, but only for first CDS
                if(CDS == CDS_list[0] and int(add_ftr["phase"]) != 0):
                    logger.info("Adjusting phase for first CDS: "+CDS)
                    add_ftr_obj["dna_sequence"] = add_ftr_obj["dna_sequence"][int(add_ftr["phase"]):]

                dna_sequence+=add_ftr_obj["dna_sequence"]
                locations.append(add_ftr_obj["location"][0])

            #translate sequence
            dna_sequence_obj = Seq(dna_sequence, IUPAC.ambiguous_dna)
            rna_sequence = dna_sequence_obj.transcribe()

            #Incomplete gene model with no start codon
            #Translate as is
            if str(rna_sequence.upper())[:3] not in codon_table.start_codons:
                logger.info("Missing start codon for "+feature_object["id"]+" Assuming incomplete gene model.")
                #temp_seq = 'AUG'+str(rna_sequence.upper())[3:]
                #rna_sequence = Seq(temp_seq, IUPAC.ambiguous_dna)

            #You should never have this problem, needs to be reported rather than "fixed"
            codon_count = len(str(rna_sequence)) % 3
            if codon_count != 0:
                logger.info("Number of bases for RNA sequence for "+feature_object["id"]+" is not divisible by 3. The resulting protein may well be mis-translated.")
                #temp_seq = str(rna_sequence.upper())+"N"
                #if codon_count == 1:
                #    temp_seq+="N"
                #new_codon_count=len(temp_seq) % 3
                #rna_sequence = Seq(temp_seq, IUPAC.ambiguous_dna)

            protein_sequence = Seq("")
            try:
                protein_sequence = rna_sequence.translate() #cds=True)
            except CodonTable.TranslationError as te:
                logger.info("TranslationError for: "+feature_object["id"],phases,exons," : "+str(te))

            cds_object["protein_translation"] = str(protein_sequence).upper()
            cds_object["protein_translation_length"]=len(cds_object["protein_translation"])
            cds_object["md5"] = hashlib.md5(cds_object["protein_translation"]).hexdigest()

            if(cds_object["protein_translation_length"] > longest_protein_length):
                longest_protein_length=cds_object["protein_translation_length"]
                longest_protein_sequence=cds_object["protein_translation"]
                    
            del cds_object["dna_sequence"]
            del cds_object["dna_sequence_length"]
            if("aliases" not in cds_object):
                cds_object["aliases"]=list()
            if("function" not in cds_object):
                cds_object["function"]=""

            #End of mRNA loop
            genome_mrnas_list.append(mrna_object)
            genome_cdss_list.append(cds_object)

        #End of gene loop
        gene_object["ontology_terms"]=dict()
        gene_object["protein_translation"]=longest_protein_sequence
        gene_object["protein_translation_length"]=longest_protein_length
        genome_features_list.append(gene_object)

    #remove sequences before loading
    for contig in assembly["contigs"]:
        del assembly["contigs"][contig]["sequence"]

#    assembly_string = simplejson.dumps(assembly, sort_keys=True, indent=4, ensure_ascii=False)
#    assembly_file = open("Bulk_Phytozome_Upload/"+assembly["name"]+'.json', 'w+')
#    assembly_file.write(assembly_string)
#    assembly_file.close()
    
    if(assembly_ref == None):
        #Upload FASTA to shock
        #Need to gunzip file first
        gunzipped_fasta_file=input_fasta_file
#        gunzipped_fasta_file=input_fasta_file[0:-3]
#        with gzip.open(input_fasta_file, 'rb') as f_in:
#            with open(gunzipped_fasta_file, 'wb') as f_out:
#                shutil.copyfileobj(f_in, f_out)

        token = os.environ.get('KB_AUTH_TOKEN') 

        logger.info("Attempting Assembly save for %s" % (assembly["assembly_id"]))
        aUtil = AssemblyUtil(callback_url)
        assembly_ref =  aUtil.save_assembly_from_fasta({'file':{'path':gunzipped_fasta_file, 'assembly_name':assembly['assembly_id']},
                                                        'workspace_name':workspace_name, 'assembly_name':assembly['assembly_id']})
        logger.info("Assembly saved for %s" % (assembly["name"]))
        
        #Remove gunzipped file
        #os.remove(input_fasta_file[0:-3])

    genome = dict()
    genome["id"]=core_genome_name
    genome["scientific_name"]=scientific_name
    genome["assembly_ref"]=assembly_ref
    genome["features"]=genome_features_list
    genome["cdss"]=genome_cdss_list
    genome["mrnas"]=genome_mrnas_list
    genome["source"]=source
    genome["domain"]="Eukaryota"
    genome["genetic_code"]=1
    genome["gc_content"]=assembly["gc_content"]
    genome["dna_size"]=assembly["dna_size"]

    if taxon_reference is not None:
        genome["taxon_ref"] = taxon_reference
        genome["taxonomy"]=taxonomy

    UserMeta = dict()
    UserMeta['Taxonomy']=taxonomy
    UserMeta['Source']=source
    UserMeta['Domain']="Eukaryota"
    UserMeta['Source ID']=core_genome_name
    UserMeta['Name']=scientific_name
    UserMeta['Genetic code']=1;

    UserMeta['GC content']=assembly["gc_content"]
    UserMeta['Size']=assembly["dna_size"]
    UserMeta['Number contigs']=assembly['num_contigs']

    #id_source_version_array = core_genome_name.split("_")
    #version = "_".join(id_source_version_array[2:])
    #UserMeta['Version']=version
    #UserMeta['url']='';

    dfUtil = DataFileUtil(callback_url)

    if(gff_handle_ref == None):
        token = os.environ.get('KB_AUTH_TOKEN') 
        file_upload = dfUtil.file_to_shock({'file_path' : input_gff_file,'make_handle': 1,'pack' : "gzip"})
        gff_handle_ref=file_upload['handle']['hid']

    genome['gff_handle_ref'] = gff_handle_ref

#    genome_string = simplejson.dumps(genome, sort_keys=True, indent=4, ensure_ascii=False)
#    genome_file = open("Bulk_Phytozome_Upload/"+core_genome_name+'.json', 'w+')
#    genome_file.write(genome_string)
#    genome_file.close()

    logger.info("Attempting Genome save for %s" % (core_genome_name))
    workspace_id = dfUtil.ws_name_to_id(workspace_name)
    genome_info =  dfUtil.save_objects({"id":workspace_id, "objects":[ {"name": core_genome_name, "type": "KBaseGenomes.Genome", "data": genome} ]})
    logger.info("Genome saved for %s" % (core_genome_name))

    return { 'genome_info': genome_info[0] }
#        'report_name': report_info['name'],
#        'report_ref': report_info['ref']

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog=__file__)

    parser.add_argument('--input_gff_file', nargs='?', help='GFF file', required=True)
    parser.add_argument('--input_fasta_file', nargs='?', help='FASTA file', required=True)
#    parser.add_argument('--input_annotation_file', nargs='?', help='Annotations file', required=False)
#    parser.add_argument('--input_protein_file', nargs='?', help='Proteins file', required=False)

    parser.add_argument('--name', nargs='?', help='Taxon', required=True)
    parser.add_argument('--version', help="Version of the gene models in the data.  Example: TAIR10", nargs='?', required=False) 
    parser.add_argument('--ws_name', nargs='?', help='Workspace to save data in', required=True)

    parser.add_argument('--taxon_ref', nargs='?', help='Taxon Workspace', required=False, default='ReferenceTaxons')
    parser.add_argument('--source', help="data source : examples Refseq, Genbank, Pythozyme, Gramene, etc", nargs='?', required=False, default="Unknown") 

    handle_service_url="https://kbase.us/services/handle_service/"
    shock_service_url="https://kbase.us/services/shock-api/"
    workspace_service_url="https://kbase.us/services/ws/"

    args, unknown = parser.parse_known_args()

    #Check files
    for check_file in (args.input_fasta_file,args.input_gff_file):
        if not os.path.isfile(check_file):
            print(check_file+" not a recognizable file")
            sys.exit(0)

        if(check_file[-3:len(check_file)] != '.gz'):
            print("{0} is not a gzipped file".format(check_file))
            subprocess.check_output(["gzip","-f",check_file])
            check_file+=".gz"

    ws_client = biokbase.workspace.client.Workspace(workspace_service_url)

    #Taxon reference
    tax_id=0
    taxon_object_name = "unknown_taxon"
    taxon_info = ws_client.get_objects([{"workspace":"ReferenceTaxons", 
                                         "name": taxon_object_name}])[0]

    taxon_ref = "%s/%s/%s" % (taxon_info['info'][6], taxon_info['info'][0], taxon_info['info'][4])
    display_sc_name = taxon_info['data']['scientific_name']
    taxonomy = taxon_info['data']['scientific_lineage']

    genome_type="User"
    upload_genome(input_gff_file=args.input_gff_file,input_fasta_file=args.input_fasta_file,workspace_name=args.ws_name,
                  shock_service_url=shock_service_url,handle_service_url=handle_service_url,workspace_service_url=workspace_service_url,
                  taxon_reference=taxon_ref,taxonomy=taxonomy,source=args.source,core_genome_name=args.name,genome_type=genome_type,scientific_name=display_sc_name)
    sys.exit(0)
