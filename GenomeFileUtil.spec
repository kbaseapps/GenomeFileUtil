/*

*/
#include <KBaseGenomes.spec>
#include <workspace.spec>

module GenomeFileUtil {

    /* A boolean - 0 for false, 1 for true.
       @range (0, 1)
    */
    typedef int boolean;

    typedef structure {
        string path;
        string shock_id;
        string ftp_url;
    } File;

    typedef mapping<string, string> usermeta;

    /*
    genome_name - becomes the name of the object
    workspace_name - the name of the workspace it gets saved to.
    source - Source of the file typically something like RefSeq or Ensembl
    taxon_ws_name - where the reference taxons are : ReferenceTaxons
    taxon_id - if defined, will try to link the Genome to the specified
        taxonomy id in lieu of performing the lookup during upload
    release - Release or version number of the data
          per example Ensembl has numbered releases of all their data: Release 31
    generate_ids_if_needed - If field used for feature id is not there,
          generate ids (default behavior is raising an exception)
    genetic_code - Genetic code of organism. Overwrites determined GC from
          taxon object
    scientific_name - will be used to set the scientific name of the genome
        and link to a taxon
    generate_missing_genes - If the file has CDS or mRNA with no corresponding
        gene, generate a spoofed gene.
    use_existing_assembly - Supply an existing assembly reference

    */
    typedef structure {
        File file;

        string genome_name;
        string workspace_name;

        string source;
        string taxon_wsname;
        string taxon_id;

        string release;
        string generate_ids_if_needed;
        int    genetic_code;
        string scientific_name;
        usermeta metadata;
        boolean generate_missing_genes;
        string use_existing_assembly;
    } GenbankToGenomeParams;

    typedef structure {
        string genome_ref;
    } GenomeSaveResult;

    funcdef genbank_to_genome(GenbankToGenomeParams params)
                returns (GenomeSaveResult result) authentication required;

    typedef structure {
        File file;
        string genome_name;

        string source;
        string taxon_wsname;
        string taxon_id;

        string release;
        string generate_ids_if_needed;
        int    genetic_code;
        string scientific_name;
        usermeta metadata;
        boolean generate_missing_genes;
        string use_existing_assembly;
    } GenbankToGenomeInput;

    typedef structure {
        int workspace_id;
        list<GenbankToGenomeInput> inputs;
    } GenbanksToGenomesParams;

    /* Results for the genbanks_to_genomes function.
        results - the results of the save operation in the same order as the input.
    */
    typedef structure {
        list<GenomeSaveResult> results;
    } GenomeSaveResults;

    funcdef genbank_to_genome(GenbanksToGenomesParams params)
                returns (GenomeSaveResults results) authentication required;

    /*
        is_gtf - optional flag switching export to GTF format (default is 0,
            which means GFF)
        target_dir - optional target directory to create file in (default is
            temporary folder with name 'gff_<timestamp>' created in scratch)
    */
    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
        boolean is_gtf;
        string target_dir;
    } GenomeToGFFParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        string file_path;
        boolean from_cache;
    } GenomeToGFFResult;

    funcdef genome_to_gff(GenomeToGFFParams params)
                returns (GenomeToGFFResult result) authentication required;

    /*
        is_gtf - optional flag switching export to GTF format (default is 0,
            which means GFF)
        target_dir - optional target directory to create file in (default is
            temporary folder with name 'gff_<timestamp>' created in scratch)
    */

    typedef structure {
        string metagenome_ref;
        list <string> ref_path_to_genome;
        boolean is_gtf;
        string target_dir;
    } MetagenomeToGFFParams;

    typedef structure {
        string file_path;
        boolean from_cache;
    } MetagenomeToGFFResult;

    typedef structure {
        string metagenome_ref;
    } MetagenomeSaveResult;

    funcdef metagenome_to_gff(MetagenomeToGFFParams params)
                returns (MetagenomeToGFFResult result) authentication required;

    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
    } GenomeToGenbankParams;

    typedef structure {
        string file_path;
    } GBFile;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        GBFile genbank_file;
        boolean from_cache;
    } GenomeToGenbankResult;

    funcdef genome_to_genbank(GenomeToGenbankParams params)
                returns (GenomeToGenbankResult result) authentication required;


    typedef structure {
        string file_path;
    } FASTAResult;

    /*
        Produce a FASTA file with the nucleotide sequences of features in a genome.

        string genome_ref: reference to a genome object
        list<string> feature_lists: Optional, which features lists (features, mrnas, cdss, non_coding_features) to provide sequences. Defaults to "features".
        list<string> filter_ids: Optional, if provided only return sequences for matching features.
        boolean include_functions: Optional, add function to header line. Defaults to True.
        boolean include_aliases: Optional, add aliases to header line. Defaults to True.
    */
    typedef structure {
        string genome_ref;
        list<string> feature_lists;
        list<string> filter_ids;
        boolean include_functions;
        boolean include_aliases;
    } GenomeFeaturesToFastaParams;

    funcdef genome_features_to_fasta(GenomeFeaturesToFastaParams params)
                returns (FASTAResult result) authentication required;

    /*
        Produce a FASTA file with the protein sequences of CDSs in a genome.

        string genome_ref: reference to a genome object
        list<string> filter_ids: Optional, if provided only return sequences for matching features.
        boolean include_functions: Optional, add function to header line. Defaults to True.
        boolean include_aliases: Optional, add aliases to header line. Defaults to True.
    */

    typedef structure {
        string genome_ref;
        list<string> filter_ids;
        boolean include_functions;
        boolean include_aliases;
    } GenomeProteinToFastaParams;

    funcdef genome_proteins_to_fasta(GenomeProteinToFastaParams params)
                returns (FASTAResult result) authentication required;


    /*  input and output structure functions for standard downloaders */
    typedef structure {
        string input_ref;
    } ExportParams;

    typedef structure {
        string shock_id;
    } ExportOutput;

    funcdef export_genome_as_genbank(ExportParams params)
                returns (ExportOutput output) authentication required;

    funcdef export_genome_as_gff(ExportParams params)
                returns (ExportOutput output) authentication required;

    funcdef export_genome_features_protein_to_fasta(ExportParams params)
                returns (ExportOutput output) authentication required;

    funcdef export_metagenome_as_gff(ExportParams params)
                returns (ExportOutput output) authentication required;

    /*
    genome_name - becomes the name of the object
    workspace_name - the name of the workspace it gets saved to.
    source - Source of the file typically something like RefSeq or Ensembl
    taxon_ws_name - where the reference taxons are : ReferenceTaxons
    taxon_id - if defined, will try to link the Genome to the specified
        taxonomy id in lieu of performing the lookup during upload
    release - Release or version number of the data
          per example Ensembl has numbered releases of all their data: Release 31
    genetic_code - Genetic code of organism. Overwrites determined GC from
          taxon object
    scientific_name - will be used to set the scientific name of the genome
        and link to a taxon
    generate_missing_genes - If the file has CDS or mRNA with no corresponding
        gene, generate a spoofed gene. Off by default
    existing_assembly_ref - a KBase assembly upa, to associate the genome with.
        Avoids saving a new assembly when specified.
    */
    typedef structure {
        File fasta_file;
        File gff_file;

        string genome_name;
        string workspace_name;

        string source;
        string taxon_wsname;
        string taxon_id;
        string release;
        int    genetic_code;
        string scientific_name;
        usermeta metadata;
        boolean generate_missing_genes;
        string existing_assembly_ref;
    } FastaGFFToGenomeParams;

    funcdef fasta_gff_to_genome(FastaGFFToGenomeParams params)
                returns (GenomeSaveResult returnVal) authentication required;

    /* As above but returns the genome instead */
    funcdef fasta_gff_to_genome_json(FastaGFFToGenomeParams params)
                returns (UnspecifiedObject genome) authentication required;

    /*
    genome_name - becomes the name of the object
    workspace_name - the name of the workspace it gets saved to.
    source - Source of the file typically something like RefSeq or Ensembl
    taxon_ws_name - where the reference taxons are : ReferenceTaxons
    taxon_id - if defined, will try to link the Genome to the specified
        taxonomy id in lieu of performing the lookup during upload
    release - Release or version number of the data
          per example Ensembl has numbered releases of all their data: Release 31
    genetic_code - Genetic code of organism. Overwrites determined GC from
          taxon object
    scientific_name - will be used to set the scientific name of the genome
        and link to a taxon
    generate_missing_genes - If the file has CDS or mRNA with no corresponding
        gene, generate a spoofed gene. Off by default
    existing_assembly_ref - a KBase assembly upa, to associate the metagenome with.
        Avoids saving a new assembly when specified.
    */

    typedef structure {
        File fasta_file;
        File gff_file;

        string genome_name;
        string workspace_name;

        string source;
        usermeta metadata;
        boolean generate_missing_genes;
        string existing_assembly_ref;
    } FastaGFFToMetagenomeParams;

    funcdef fasta_gff_to_metagenome(FastaGFFToMetagenomeParams params)
                returns (MetagenomeSaveResult returnVal) authentication required;

    typedef structure {
        string workspace;
        string name;
        KBaseGenomes.Genome data;
        boolean hidden;
        boolean upgrade;
    } SaveOneGenomeParams;

    typedef structure {
        Workspace.object_info info;
    } SaveGenomeResult;

    funcdef save_one_genome(SaveOneGenomeParams params)
                returns (SaveGenomeResult returnVal) authentication required;

    /*
    gff_file - object containing path to gff_file
    ws_ref - input Assembly or Genome reference

    genome_name - becomes the name of the object
    workspace_name - the name of the workspace it gets saved to.
    source - Source of the file typically something like RefSeq or Ensembl
    taxon_ws_name - where the reference taxons are : ReferenceTaxons
    taxon_id - if defined, will try to link the Genome to the specified
        taxonomy id in lieu of performing the lookup during upload
    release - Release or version number of the data
          per example Ensembl has numbered releases of all their data: Release 31
    genetic_code - Genetic code of organism. Overwrites determined GC from
          taxon object
    scientific_name - will be used to set the scientific name of the genome
        and link to a taxon
    metadata - any user input metadata
    generate_missing_genes - If the file has CDS or mRNA with no corresponding
        gene, generate a spoofed gene. Off by default
    */

    typedef structure {
        string ws_ref;
        File gff_file;

        string genome_name;
        string workspace_name;

        string source;
        string taxon_wsname;
        string taxon_id;
        string release;
        int    genetic_code;
        string scientific_name;
        usermeta metadata;
        boolean generate_missing_genes;

    } WsObjGFFToGenomeParams;

    /*
    This function takes in a workspace object of type KBaseGenomes.Genome or KBaseGenomeAnnotations.Assembly and a gff file and produces a KBaseGenomes.Genome reanotated according to the the input gff file.
    */
    funcdef ws_obj_gff_to_genome(WsObjGFFToGenomeParams params)
                returns (GenomeSaveResult returnVal) authentication required;

    /*
    gff_file - object containing path to gff_file
    ws_ref - input Assembly or AnnotatedMetagenomeAssembly reference

    genome_name - becomes the name of the object
    workspace_name - the name of the workspace it gets saved to.
    source - Source of the file typically something like RefSeq or Ensembl

    genetic_code - Genetic code of organism. Overwrites determined GC from
          taxon object
    metadata - any user input metadata
    generate_missing_genes - If the file has CDS or mRNA with no corresponding
        gene, generate a spoofed gene. Off by default
    */

    typedef structure {
        string ws_ref;
        File gff_file;

        string genome_name;
        string workspace_name;

        string source;
        usermeta metadata;
        boolean generate_missing_genes;

    } WsObjGFFToMetagenomeParams;

    /*
    This function takes in a workspace object of type KBaseMetagenomes.AnnotatedMetagenomeAssembly or KBaseGenomeAnnotations.Assembly and a gff file and produces a KBaseMetagenomes.AnnotatedMetagenomeAssembly reanotated according to the the input gff file.
    */
    funcdef ws_obj_gff_to_metagenome(WsObjGFFToMetagenomeParams params)
                returns (MetagenomeSaveResult returnVal) authentication required;

    /*
    Parameters for the update_taxon_assignments function.
    Fields:
        workspace_id: a workspace UPA of a Genome object
        taxon_assignments: an optional mapping of assignments to add or replace. This will perform a
            merge on the existing assignments. Any new assignments are added, while any existing
            assignments are replaced.
        remove_assignments: an optional list of assignment names to remove.

    @optional taxon_assignments remove_assignments
    */
    typedef structure {
        int workspace_id;
        int object_id;
        mapping<string, string> taxon_assignments;
        list<string> remove_assignments;
    } UpdateTaxonAssignmentsParams;

    /*
    Result of the update_taxon_assignments function.
    Fields:
        ws_obj_ref: a workspace UPA of a Genome object
    */
    typedef structure {
        string ws_obj_ref;
    } UpdateTaxonAssignmentsResult;

    /*
    Add, replace, or remove taxon assignments for a Genome object.
    */
    funcdef update_taxon_assignments(UpdateTaxonAssignmentsParams params)
        returns (UpdateTaxonAssignmentsResult returnVal) authentication required;

};
