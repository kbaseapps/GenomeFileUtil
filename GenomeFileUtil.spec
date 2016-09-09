/*
** A KBase module: GenomeFileUtil
**
** This module contains methods for converting KBase Data Objects to common bioinformatics file formats
** 
*/
module GenomeFileUtil {

    /* A boolean - 0 for false, 1 for true.
       @range (0, 1)
    */
    typedef int boolean;

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects2([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string data_obj_name;
    typedef string data_obj_ref;


    /*
    ** characteristics of genomes, files, and methods
    */
    typedef string sequence;
    typedef string feature_id;
    typedef string genome_id;
    typedef string path_type;
    typedef string pattern_type;
    typedef string log_msg;
    typedef string true_false;


    typedef structure {
        string path;
        string shock_id;
        string ftp_url;
    } File;

    /* */
    typedef structure {
        File file;

        string genome_name;
        string workspace_name;

        string source;
        string taxon_wsname;

    } GenbankToGenomeParams;

    typedef structure {
        string genome_ref;
    } GenomeSaveResult;

    funcdef genbank_to_genome(GenbankToGenomeParams params)
                returns (GenomeSaveResult result) authentication required;




    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
    } GenomeToGFFParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        File gff_file;
        boolean from_cache;
    } GenomeToGFFResult;

    funcdef genome_to_gff(GenomeToGFFParams params)
                returns (GenomeToGFFResult result) authentication required;



    /* TranslateNucToProtSeq()
    */
    typedef structure {
        string  nuc_seq;
        string  genetic_code;
    } TranslateNucToProtSeq_Params;

    typedef structure {
        string  prot_seq;
    } TranslateNucToProtSeq_Output;

    funcdef TranslateNucToProtSeq (TranslateNucToProtSeq_Params params)  returns (TranslateNucToProtSeq_Output) authentication required;


    /* GenomeToFASTA()
    */
    typedef structure {
        data_obj_ref   genome_ref;
        path_type      file;
        path_type      dir;
        list<log_msg>  console;
        list<log_msg>  invalid_msgs;
        string         residue_type;
        string         feature_type;
        pattern_type   record_id_pattern;
        pattern_type   record_desc_pattern;
        string         case;
        int            linewrap;
    } GenomeToFASTA_Params;

    typedef structure {
        path_type         fasta_file_path;
        list<feature_id>  feature_ids;
    } GenomeToFASTA_Output;

    funcdef GenomeToFASTA (GenomeToFASTA_Params params)  returns (GenomeToFASTA_Output) authentication required;


    /* GenomeSetToFASTA()
    */
    typedef structure {
        data_obj_ref   genomeSet_ref;
        path_type      file;
        path_type      dir;
        list<log_msg>  console;
        list<log_msg>  invalid_msgs;
        string         residue_type;
        string         feature_type;
        pattern_type   record_id_pattern;
        pattern_type   record_desc_pattern;
        string         case;
        int            linewrap;
        true_false     merge_fasta_files;
    } GenomeSetToFASTA_Params;

    typedef structure {
        list<path_type>                       fasta_file_path_list;
        mapping<genome_id, list<feature_id>>  feature_ids_by_genome_id;
    } GenomeSetToFASTA_Output;

    funcdef GenomeSetToFASTA (GenomeSetToFASTA_Params params)  returns (GenomeSetToFASTA_Output) authentication required;


    /* FeatureSetToFASTA()
    */
    typedef structure {
        data_obj_ref   featureSet_ref;
        path_type      file;
        path_type      dir;
        list<log_msg>  console;
        list<log_msg>  invalid_msgs;
        string         residue_type;
        string         feature_type;
        pattern_type   record_id_pattern;
        pattern_type   record_desc_pattern;
        string         case;
        int            linewrap;
    } FeatureSetToFASTA_Params;

    typedef structure {
        path_type                                fasta_file_path;
        mapping<data_obj_ref, list<feature_id>>  feature_ids_by_genome_ref;
    } FeatureSetToFASTA_Output;

    funcdef FeatureSetToFASTA (FeatureSetToFASTA_Params params)  returns (FeatureSetToFASTA_Output) authentication required;

};
