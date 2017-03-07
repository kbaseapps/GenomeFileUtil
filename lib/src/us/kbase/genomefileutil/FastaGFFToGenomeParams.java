
package us.kbase.genomefileutil;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: FastaGFFToGenomeParams</p>
 * <pre>
 * genome_name - becomes the name of the object
 * workspace_name - the name of the workspace it gets saved to.
 * source - Source of the file typically something like RefSeq or Ensembl
 * taxon_ws_name - where the reference taxons are : ReferenceTaxons
 * taxon_reference - if defined, will try to link the Genome to the specified
 * taxonomy object insteas of performing the lookup during upload
 * release - Release or version number of the data 
 *   per example Ensembl has numbered releases of all their data: Release 31
 * generate_ids_if_needed - If field used for feature id is not there, 
 *   generate ids (default behavior is raising an exception)
 * genetic_code - Genetic code of organism. Overwrites determined GC from 
 *   taxon object
 * type - Reference, Representative or User upload
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_file",
    "gff_file",
    "genome_name",
    "workspace_name",
    "source",
    "taxon_wsname",
    "taxon_reference",
    "release",
    "genetic_code",
    "type",
    "metadata"
})
public class FastaGFFToGenomeParams {

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("fasta_file")
    private File fastaFile;
    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("gff_file")
    private File gffFile;
    @JsonProperty("genome_name")
    private java.lang.String genomeName;
    @JsonProperty("workspace_name")
    private java.lang.String workspaceName;
    @JsonProperty("source")
    private java.lang.String source;
    @JsonProperty("taxon_wsname")
    private java.lang.String taxonWsname;
    @JsonProperty("taxon_reference")
    private java.lang.String taxonReference;
    @JsonProperty("release")
    private java.lang.String release;
    @JsonProperty("genetic_code")
    private Long geneticCode;
    @JsonProperty("type")
    private java.lang.String type;
    @JsonProperty("metadata")
    private Map<String, String> metadata;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("fasta_file")
    public File getFastaFile() {
        return fastaFile;
    }

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("fasta_file")
    public void setFastaFile(File fastaFile) {
        this.fastaFile = fastaFile;
    }

    public FastaGFFToGenomeParams withFastaFile(File fastaFile) {
        this.fastaFile = fastaFile;
        return this;
    }

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("gff_file")
    public File getGffFile() {
        return gffFile;
    }

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("gff_file")
    public void setGffFile(File gffFile) {
        this.gffFile = gffFile;
    }

    public FastaGFFToGenomeParams withGffFile(File gffFile) {
        this.gffFile = gffFile;
        return this;
    }

    @JsonProperty("genome_name")
    public java.lang.String getGenomeName() {
        return genomeName;
    }

    @JsonProperty("genome_name")
    public void setGenomeName(java.lang.String genomeName) {
        this.genomeName = genomeName;
    }

    public FastaGFFToGenomeParams withGenomeName(java.lang.String genomeName) {
        this.genomeName = genomeName;
        return this;
    }

    @JsonProperty("workspace_name")
    public java.lang.String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(java.lang.String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public FastaGFFToGenomeParams withWorkspaceName(java.lang.String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("source")
    public java.lang.String getSource() {
        return source;
    }

    @JsonProperty("source")
    public void setSource(java.lang.String source) {
        this.source = source;
    }

    public FastaGFFToGenomeParams withSource(java.lang.String source) {
        this.source = source;
        return this;
    }

    @JsonProperty("taxon_wsname")
    public java.lang.String getTaxonWsname() {
        return taxonWsname;
    }

    @JsonProperty("taxon_wsname")
    public void setTaxonWsname(java.lang.String taxonWsname) {
        this.taxonWsname = taxonWsname;
    }

    public FastaGFFToGenomeParams withTaxonWsname(java.lang.String taxonWsname) {
        this.taxonWsname = taxonWsname;
        return this;
    }

    @JsonProperty("taxon_reference")
    public java.lang.String getTaxonReference() {
        return taxonReference;
    }

    @JsonProperty("taxon_reference")
    public void setTaxonReference(java.lang.String taxonReference) {
        this.taxonReference = taxonReference;
    }

    public FastaGFFToGenomeParams withTaxonReference(java.lang.String taxonReference) {
        this.taxonReference = taxonReference;
        return this;
    }

    @JsonProperty("release")
    public java.lang.String getRelease() {
        return release;
    }

    @JsonProperty("release")
    public void setRelease(java.lang.String release) {
        this.release = release;
    }

    public FastaGFFToGenomeParams withRelease(java.lang.String release) {
        this.release = release;
        return this;
    }

    @JsonProperty("genetic_code")
    public Long getGeneticCode() {
        return geneticCode;
    }

    @JsonProperty("genetic_code")
    public void setGeneticCode(Long geneticCode) {
        this.geneticCode = geneticCode;
    }

    public FastaGFFToGenomeParams withGeneticCode(Long geneticCode) {
        this.geneticCode = geneticCode;
        return this;
    }

    @JsonProperty("type")
    public java.lang.String getType() {
        return type;
    }

    @JsonProperty("type")
    public void setType(java.lang.String type) {
        this.type = type;
    }

    public FastaGFFToGenomeParams withType(java.lang.String type) {
        this.type = type;
        return this;
    }

    @JsonProperty("metadata")
    public Map<String, String> getMetadata() {
        return metadata;
    }

    @JsonProperty("metadata")
    public void setMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
    }

    public FastaGFFToGenomeParams withMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((((((("FastaGFFToGenomeParams"+" [fastaFile=")+ fastaFile)+", gffFile=")+ gffFile)+", genomeName=")+ genomeName)+", workspaceName=")+ workspaceName)+", source=")+ source)+", taxonWsname=")+ taxonWsname)+", taxonReference=")+ taxonReference)+", release=")+ release)+", geneticCode=")+ geneticCode)+", type=")+ type)+", metadata=")+ metadata)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
