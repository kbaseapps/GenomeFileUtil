"""
GenomeAnnotation to GenBank file conversion.
"""
__author__ = 'Marcin Joachimiak <mjoachimiak@lbl.gov>'
__date__ = '6/10/16'

# Stdlib
from collections import OrderedDict
import os
import string
import cStringIO as StringIO

from Bio import SeqIO

# Local
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI as GenomeAnnotationAPI_local
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI

STD_PREFIX = " " * 21


class GenomeToGenbank(object):

    def __init__(self, sdk_config):
        self.cfg = sdk_config

    def validate_params(self, params):
        if 'genome_ref' not in params:
            raise ValueError('required "genome_ref" field was not defined')

    def export(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome genbank handle reference
        getGenomeOptions = {
            'genomes': [{
                'ref': params['genome_ref']
            }],
            'included_fields': ['genbank_handle_ref'],
            'ignore_errors': 0  # if we can't find the genome, throw an error
        }
        if 'ref_path_to_genome' in params:
            getGenomeOptions['genomes'][0]['ref_path_to_genome'] = params['ref_path_to_genome']

        api = GenomeAnnotationAPI(self.cfg.callbackURL)
        genome_data = api.get_genome_v1(getGenomeOptions)['genomes'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split('-')[0] != 'KBaseGenomes.Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        # 4) build the genbank file and return it
        print('not cached, building file...')
        result = self.build_genbank_file(getGenomeOptions, "KBase_derived_" + info[1] + ".gbff")
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result

    def export_original_genbank(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome genbank handle reference
        getGenomeOptions = {
            'genomes':[{
                'ref': params['genome_ref']
            }],
            'included_fields':['genbank_handle_ref'],
            'ignore_errors':0 # if we can't find the genome, throw an error
        }
        if 'ref_path_to_genome' in params:
            getGenomeOptions['genomes'][0]['ref_path_to_genome'] = params['ref_path_to_genome']

        api = GenomeAnnotationAPI(self.cfg.callbackURL)
        genome_data = api.get_genome_v1(getGenomeOptions)['genomes'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split('-')[0] != 'KBaseGenomes.Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        # 4) if the genbank handle is there, get it and return
        print('checking if genbank file is cached...')
        result = self.get_genbank_handle(data)
        return result

    def get_genbank_handle(self, data):
        if 'genbank_handle_ref' not in data:
            return None
        if data['genbank_handle_ref'] is None:
            return None

        print('pulling cached genbank file from Shock: '+str(data['genbank_handle_ref']))
        dfu = DataFileUtil(self.cfg.callbackURL)
        file = dfu.shock_to_file({
                            'handle_id': data['genbank_handle_ref'],
                            'file_path': self.cfg.sharedFolder,
                            'unpack': 'unpack'
                        })
        return {
            'genbank_file': {
                'file_path': file['file_path']
            }
        }

    def build_genbank_file(self, getGenomeOptions, output_filename):
        ref = getGenomeOptions["genomes"][0]["ref"]
        services = {
            "workspace_service_url": self.cfg.workspaceURL,
            "shock_service_url": self.cfg.shockURL,
            "handle_service_url": self.cfg.handleURL
        }

        genome = GenomeAnnotationAPI_local(services=services,
                                           token=os.environ['KB_AUTH_TOKEN'],
                                           ref=ref)
        g = GenbankAnnotations(genome, self.cfg.callbackURL)
        g.to_file(output_filename)

        return {
            'genbank_file': {
                'file_path': str(output_filename)
            }
        }


class GenbankAnnotations(object):
    def __init__(self, genome=None, callbackURL=''):
        self._contents = StringIO.StringIO()
        self._ga = genome
        print('downloading assembly')
        au = AssemblyUtil(callbackURL)
        assembly_info = self._ga.get_assembly().get_info()
        assembly_ref = str(assembly_info['workspace_id']) + '/' + str(assembly_info['object_id']) + '/' + str(assembly_info['version'])
        print('Assembly reference = ' + assembly_ref)
        assembly_file_path = au.get_assembly_as_fasta({'ref': assembly_ref})['path']

        print('extracting taxonomy information')
        self._taxa = self._ga.get_taxon()
        self._tax_lineage = self._taxa.get_scientific_lineage()

        print('assembling feature and protein data')
        self._genome_name = str(self._ga.get_id())
        self._proteins = self._ga.get_proteins()
        self._features = self._ga.get_features()

        print('writing file')
        # read in fasta file to build the contig index
        contigs = {}
        for record in SeqIO.parse(assembly_file_path, "fasta"):
            # SeqRecord(seq=Seq('TTAT...', SingleLetterAlphabet()),
            #           id='gi|113968346|ref|NC_008321.1|',
            #           name='gi|113968346|ref|NC_008321.1|',
            #           description='gi|113968346|ref|NC_008321.1| Shewanella sp. MR-4 chromosome, complete genome',
            #           dbxrefs=[])
            contigs[record.id] = {
                'length': len(record.seq),
                'sequence': str(record.seq)
            }
        self._contigs = contigs
        contig_length_dict = dict()
        for contig_id in contigs:
            contig_length_dict[contig_id] = contigs[contig_id]["length"]
        del contigs
        # longest contigs first
        contigs_tuples = sorted(contig_length_dict.items(), key=lambda x: x[1],
                                reverse=True)
        # print("Contig tuples : " + str(contigs_tuples))

        # organize features by location
        feature_ids_by_region = self._ga.get_feature_ids(group_by="region")["by_region"]
        # print('FEATURE IDS BY REGION :: ' + str(feature_ids_by_region))

        for cid, length in contigs_tuples:
            self._add_contig_header(cid)
            if cid not in feature_ids_by_region:
                continue
            # get the features from both strands, sort them by the starting loc
            features = sorted(feature_ids_by_region[cid].get("+", {}).items() +
                              feature_ids_by_region[cid].get("-", {}).items(),
                              key=lambda x: x[0].split("_")[0])
            for loc, fid_list in features:
                for fid in fid_list:
                    self._add_feature(fid)

    # formats a string into lines of given length (first line can be different)
    @staticmethod
    def _format_protein_sequence(s):
        if not s:
            raise ValueError("Protein sequence not valid")

        s = s.replace("\"", "")
        out = []
        start = 0
        last = 0
        firstChunk = 43
        chunkSize = 57
        size = len(s)

        if size < firstChunk:
            return s + "\n"

        # create the first line, which is shorter
        out.append(s[start:firstChunk] + "\n")
        start += firstChunk

        while start + chunkSize < size:
            out.append(STD_PREFIX + s[start:start + chunkSize] + "\n")
            start += chunkSize

        remainder = s[start:]
        # check to see if the last line has sequence to the end, if so split
        if len(remainder) == chunkSize:
            out.append(STD_PREFIX + remainder[:-3] + "\n")
            out.append(STD_PREFIX + remainder[-3:])
        else:
            out.append(STD_PREFIX + remainder)

        formattedProtein = "".join(out)
        if len(formattedProtein.strip()) == 0:
            formattedProtein = ""

        return formattedProtein

    @staticmethod
    def _format_dna_sequence(s, charnum, linenum):
        out = StringIO.StringIO()
        # e.g., "        1 tctcgcagag ttcttttttg tattaacaaa cccaaaaccc atagaattta atgaacccaa\n"
        # start at position 1 of the sequence
        out.write("        1 ")
        index = 0
        size = len(s)
        end = 0
        while end < size:
            for n in xrange(linenum / charnum):
                # compute the boundary of the chunk
                end = index + charnum
                # if this runs past the end of the overall sequence length, take the remainder and exit
                if end > size:
                    out.write(s[index:])
                    index = size
                    break
                else:
                    out.write("{} ".format(s[index:index + charnum]))
                    index += charnum
            # add a line break
            out.write("\n")
            if index >= size:
                break
            # if we haven't reached the end, write the current index starting from 1
            indexString = str(index + 1)
            out.write("{}{} ".format(" " * (9 - len(indexString)), indexString))
        return out.getvalue()

    def _sort_feature_ids(self, feature_ids=None):
        def rank_type(t):
            type_sorting = {"gene": 1, "CDS": 2}
            if t in type_sorting:
                return type_sorting[t]
            else:
                return 3

        # first sort by location
        location_sorted = sorted(feature_ids, cmp=lambda x,y: cmp(self._features[x]["feature_locations"][0]["start"],
                                                                  self._features[y]["feature_locations"][0]["start"]))

        # now sort by type
        type_sorted = sorted(location_sorted, cmp=lambda x,y: cmp(rank_type(self._features[x]["feature_type"]),
                                                                  rank_type(self._features[y]["feature_type"])))

        return type_sorted

    def _add_feature(self, feature_id=None):
        # may want to have different behavior depending on feature type
        self._add_feature_data(feature_id)

    def _add_feature_location(self, feature_locations=None):
        added = 0
        complement = False
        isJoin = False
        numloc = len(feature_locations)
        for n in xrange(numloc):
            now4 = feature_locations[n]
            if added == 0 and now4['strand'] == "-":
                self._contents.write("complement(")
                complement = True

            if len(feature_locations) > 1:
                if added == 0:
                    self._contents.write("join(")
                isJoin = True

            if not complement:
                start = now4['start']
                stop = now4['start'] + now4['length'] - 1
            else:
                start = now4['start'] - now4['length'] + 1
                stop = now4['start']

            self._contents.write("{}..{}".format(start, stop))

            if numloc > 0 and n < (numloc - 1):
                self._contents.write(",")

            added += 1

        tail = "\n"

        if complement:
            tail = ")" + tail
        if isJoin:
            tail = ")" + tail

        self._contents.write(tail)

    def _add_feature_function(self, function=None):
        # line length is 80 characters, - 21 for the prefix, -1 for the newline
        size = 58
        # max length of first line
        first_size = 46

        if not function:
            return

        if len(function) < first_size:
            self._contents.write("{}/function=\"{}\"\n".format(STD_PREFIX, function))
            return

        lines = []
        position = min(first_size, len(function) - 1)
        while position > 0 and function[position - 1] != " ":
            position -= 1
        if position == 0:
            position = min(first_size, len(function) - 1)
        first_line = function[:position]
        lines.append(first_line + "\n")
        remainder = function[position:].rstrip()

        if len(remainder) >= size:
            tokens = remainder.split()
            current_line = ""
            for t in tokens:
                current_token = t

                if len(current_line) + len(current_token) + 1 >= size:
                    lines.append(current_line.rstrip() + "\n")
                    current_line = ""

                    # split up all tokens larger than a line and append
                    while len(current_token) > size:
                        lines.append(current_token[:size] + "\n")
                        current_token = current_token[size:]

                current_line += current_token + " "

            if len(current_line) > 0:
                lines.append(current_line)
            lines[-1] = lines[-1].rstrip()
        elif len(remainder) > 0:
            lines.append(remainder)

        format_function = STD_PREFIX.join(lines)
        self._contents.write("{}/function=\"{}\"\n".format(STD_PREFIX, format_function))

    def _add_feature_aliases(self, aliases=None):
        if aliases is not None and len(aliases) > 0:
            # TODO handle different alias cases and types
            for id in aliases:
                for source in aliases[id]:
                    if source.find("Genbank ") == -1:
                        self._contents.write("{}/db_xref=\"{}:{}\"\n".format(STD_PREFIX, source, id))
                    else:
                        key = source.replace("Genbank ", "")
                        self._contents.write("{}/{}=\"{}\"\n".format(STD_PREFIX, key, id))

    def _add_contig_header(self, contig_id=None):
        if contig_id not in self._contigs:
            raise ValueError('ContigID (' + str(contig_id) + ') was not found in the assembly.  This may be a corrupted genome.')
        contig_id_field_length = 16
        bp_field_length = 11
        temp_contig_id = contig_id
        if len(contig_id) > contig_id_field_length:
            temp_contig_id = contig_id[-16:]
        contig_spaces_padding = contig_id_field_length - len(temp_contig_id)
        bp_field_padding = bp_field_length - len(str(self._contigs[contig_id]["length"]))
        self._contents.write("LOCUS{}{}{} {}{} bp    DNA{}\n".format(" " * 7,
                                                                     temp_contig_id,
                                                                     " " * contig_spaces_padding,
                                                                     " " * bp_field_padding,
                                                                     self._contigs[contig_id]["length"],
                                                                     " " * 27))
        sn = self._taxa.get_scientific_name()
        self._contents.write("DEFINITION  {} genome.\n".format(sn))
        self._contents.write("SOURCE      {}\n".format(sn))
        self._contents.write("  ORGANISM  {}\n".format(sn))

        if self._tax_lineage:
            chunkSize = 68
            lineage = []
            prefix = " " * 12
            currentLine = ""
            for ancestor in self._tax_lineage:
                if len(ancestor + ";") >= chunkSize:
                    remainder = chunkSize - len(currentLine)
                    lineage.append(prefix + currentLine + ancestor[:remainder])
                    currentLine = ""
                    start = remainder
                    while start < len(ancestor):
                        lineage.append(prefix + ancestor[start:start + chunkSize])
                        start += chunkSize
                elif len(currentLine + ancestor + ";") >= chunkSize:
                    lineage.append(prefix + currentLine)
                    currentLine = ancestor + ";"
                else:
                    currentLine += ancestor + ";"
            if len(lineage) > 0:
                if prefix + currentLine != lineage[-1]:
                    lineage.append(prefix + currentLine)
                if lineage[-1][-1] == ';':
                    lineage[-1] = lineage[-1][:-1]
                self._contents.write("{}.\n".format("\n".join(lineage)))

        self._contents.write("COMMENT{}Exported from the DOE KnowledgeBase.\n".format(" " * 5))
        self._contents.write("FEATURES{}Location/Qualifiers\n".format(" " * 13))
        self._contents.write("{}source{}1..{}\n".format(" " * 5, " " * 10, self._contigs[contig_id]["length"]))
        self._contents.write("{}/organism=\"{}\"\n".format(STD_PREFIX, self._taxa.get_scientific_name()))
        self._contents.write("{}/mol_type=\"DNA\"\n".format(STD_PREFIX))

    def _add_contig_sequence(self, contig_id=None):
        self._contents.write("ORIGIN\n")
        self._contents.write(self._format_dna_sequence(self._contigs[contig_id]['sequence'], 10, 60))
        self._contents.write("//\n")

    def _add_feature_data(self, feature_id=None):
        data = self._features[feature_id]
        self._contents.write("{}{}{}".format(" " * 5, data["feature_type"],
                                             " " * (16 - len(data["feature_type"]))))
        self._add_feature_location(data["feature_locations"])
        self._contents.write("{}/kbase_id=\"{}\"\n".format(STD_PREFIX, feature_id))
        self._add_feature_function(data['feature_function'])
        self._add_feature_aliases(data['feature_aliases'])

        if data["feature_type"] == "CDS":
            self._contents.write("{}/protein_id=\"{}\"\n".format(STD_PREFIX, feature_id))

            if feature_id in self._proteins:
                protein_translation = self._proteins[feature_id]['protein_amino_acid_sequence']
                if not protein_translation:
                    protein_translation_final = ""
                else:
                    protein_translation_final = self._format_protein_sequence(protein_translation)

                self._contents.write("{}/translation=\"{}\"\n".format(STD_PREFIX, protein_translation_final))

        #     gene            3631..5899
        #             /gene="NAC001"
        #             /locus_tag="AT1G01010"
        #             /gene_synonym="ANAC001; NAC domain containing protein 1;
        #             NAC001; T25K16.1; T25K16_1"
        #             /db_xref="GeneID:839580"
        #
        #     CDS             complement(join(6915..7069,7157..7232,7384..7450,
        # 7564..7649,7762..7835,7942..7987,8236..8325,8417..8464,
        # 8571..8666))
        # /gene="ARV1"
        # /locus_tag="AT1G01020"
        # /gene_synonym="T25K16.2; T25K16_2"
        # /inference="Similar to DNA sequence:INSD:AY758070.1"
        # /inference="Similar to RNA sequence,
        # EST:INSD:EH846835.1,INSD:AU227271.1,INSD:EG512767.1,
        # INSD:EG452794.1,INSD:EL102675.1,INSD:EG512768.1,
        # INSD:EH959539.1"
        # /note="ARV1; CONTAINS InterPro DOMAIN/s: Arv1-like protein
        # (InterPro:IPR007290); BEST Arabidopsis thaliana protein
        # match is: Arv1-like protein (TAIR:AT4G01510.1); Has 311
        # Blast hits to 311 proteins in 154 species: Archae - 0;
        # Bacteria - 0; Metazoa - 110; Fungi - 115; Plants - 42;
        # Viruses - 0; Other Eukaryotes - 44 (source: NCBI BLink)."
        # /codon_start=1
        # /product="protein Arv1"
        # /protein_id="NP_171610.2"
        # /db_xref="GI:79332834"
        # /db_xref="GeneID:839569"
        # /db_xref="TAIR:AT1G01020"
        # /translation="MAASEHRCVGCGFRVKSLFIQYSPGNIRLMKCGNCKEVADEYIE
        # CERMIIFIDLILHRPKVYRHVLYNAINPATVNIQHLLWKLVFAYLLLDCYRSLLLRKS
        # DEESSFSDSPVLLSIKVLIGVLSANAAFIISFAIATKGLLNEVSRRREIMLGIFISSY
        # FKIFLLAMLVWEFPMSVIFFVDILLLTSNSMALKVMTESTMTRCIAVCLIAHLIRFLV
        # GQIFEPTIFLIQIGSLLQYMSYFFRIV"
        #
        # mRNA            complement(join(5928..6263,6437..7069,7157..7232,
        # 7384..7450,7564..7649,7762..7835,7942..7987,8236..8325,
        # 8417..8464,8571..8737))
        # /gene="ARV1"
        # /locus_tag="AT1G01020"
        # /gene_synonym="T25K16.2; T25K16_2"
        # /product="protein Arv1"
        # /inference="Similar to DNA sequence:INSD:AY758070.1"
        # /inference="Similar to RNA sequence,
        # EST:INSD:EH846835.1,INSD:AU227271.1,INSD:EG512767.1,
        # INSD:EG452794.1,INSD:EL102675.1,INSD:EG512768.1,
        # INSD:EH959539.1"
        # /transcript_id="NM_099984.5"
        # /db_xref="GI:240253989"

    def to_file(self, outfile=None):
        valid_chars = "-_.(){0}{1}".format(string.ascii_letters, string.digits)
        filename_chars = []

        if outfile is None:
            for character in self._genome_name:
                if character in valid_chars:
                    filename_chars.append(character)
                else:
                    filename_chars.append("_")

                    outfile = "".join(filename_chars) + ".gbff"

        with open(outfile, 'w') as f:
            f.write(self._contents.getvalue())
