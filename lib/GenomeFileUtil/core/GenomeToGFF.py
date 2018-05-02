from collections import defaultdict
import csv
import os
import time
import urllib

from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeUtils import get_start, get_end
import traceback


class GenomeToGFF:
    '''
    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
        int is_gtf;
    } GenomeToGFFParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        File file_path;
        boolean from_cache;
    } GenomeToGFFResult;

    funcdef genome_to_gff(GenomeToGFFParams params)
                returns (GenomeToGFFResult result) authentication required;
    '''

    def __init__(self, sdk_config):
        self.cfg = sdk_config
        self.dfu = DataFileUtil(self.cfg.callbackURL)
        self.child_dict = {}
        self.transcript_counter = defaultdict(int)

    def export(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome info
        genome_data = self.dfu.get_objects({
            'object_refs': [params['genome_ref']]
        })['data'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split(".")[1].split('-')[0] != 'Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        is_gtf = params.get('is_gtf', 0)

        target_dir = params.get('target_dir')
        if not target_dir:
            target_dir = os.path.join(self.cfg.sharedFolder, "gff_" + str(int(time.time() * 1000)))
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        # 4) if the GFF handle is there, get it and return
        if is_gtf != 1:
            print('checking if GFF file is cached...')
            result = self.get_gff_handle(data, target_dir)
            if result is not None:
                result['from_cache'] = 1
                return result
            print('not cached, building file...')

        # 5) otherwise, build the GFF/GTF file and return it
        result = self.build_gff_file(data, target_dir, info[1], is_gtf == 1)
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result

    def get_gff_handle(self, data, output_dir):

        if 'gff_handle_ref' not in data:
            return None
        if data['gff_handle_ref'] is None:
            return None

        print('pulling cached GFF file from Shock: '+str(data['gff_handle_ref']))
        file_ret = self.dfu.shock_to_file(
            {'handle_id': data['gff_handle_ref'],
             'file_path': output_dir,
             'unpack': 'unpack'})
        return {'file_path': file_ret['file_path']}

    def build_gff_file(self, genome_data, output_dir, output_filename, is_gtf):
        def feature_sort(feat):
            order = ('gene', 'mRNA', 'CDS')
            if feat.get('children'):
                priority = 0
            elif feat['type'] not in order:
                priority = len(order)
            else:
                priority = order.index(feat['type'])
            return get_start(self.get_common_location(
                feat['location'])), priority

        gff_header = ['seqname', 'source', 'type', 'start', 'end', 'score',
                      'strand', 'frame', 'attribute']

        # create the file
        file_ext = ".gtf" if is_gtf else ".gff"
        out_file_path = os.path.join(output_dir, output_filename + file_ext)
        print('Creating file: ' + str(out_file_path))

        """There is two ways of printing, if a feature has a parent_gene, it 
        will be printed breadth first when it's parent parent gene is printed.
        if not, it needs to be added to the features_by_contig to be printed"""
        # sort every feature in the feat_arrays into a dict by contig
        features_by_contig = defaultdict(list)
        for feature in genome_data['features'] + genome_data.get(
                'non_coding_features', []):
            # type is not present in new gene array
            if 'type' not in feature:
                feature['type'] = 'gene'
            features_by_contig[feature['location'][0][0]].append(feature)

        for mrna in genome_data.get('mrnas', []):
            mrna['type'] = 'mRNA'
            if mrna.get('parent_gene'):
                self.child_dict[mrna['id']] = mrna
            else:
                features_by_contig[mrna['location'][0][0]].append(mrna)

        for cds in genome_data.get('cdss', []):
            cds['type'] = 'CDS'
            if cds.get('parent_gene'):
                self.child_dict[cds['id']] = cds
            else:
                features_by_contig[cds['location'][0][0]].append(cds)

        file_handle = open(out_file_path, 'w')
        writer = csv.DictWriter(file_handle, gff_header, delimiter="\t",
                                escapechar='\\', quotechar="'")
        for contig in genome_data.get('contig_ids', features_by_contig.keys()):
            file_handle.write("##sequence-region {}\n".format(contig))
            features_by_contig[contig].sort(key=feature_sort)
            for feature in features_by_contig[contig]:
                writer.writerows(self.make_feature_group(feature, is_gtf))

        return {'file_path': out_file_path}

    def make_feature_group(self, feature, is_gtf):
        # RNA types make exons if they have compound locations
        if feature['type'] in {'RNA', 'mRNA', 'tRNA', 'rRNA', 'misc_RNA', 'transcript'}:
            loc = self.get_common_location(feature['location'])
            lines = [self.make_feature(loc, feature, is_gtf)]
            for i, loc in enumerate(feature['location']):
                exon = {'id': "{}_exon_{}".format(feature['id'], i + 1),
                        'parent_gene': feature.get('parent_gene', ""),
                        'parent_mrna': feature['id']}
                lines.append(self.make_feature(loc, exon, is_gtf))
        # other types duplicate the feature
        else:
            lines = [self.make_feature(loc, feature, is_gtf)
                     for loc in feature['location']]

        #if this is a gene with mRNAs, make the mrna (and subfeatures)
        if feature.get('mrnas', False):
            for mrna_id in feature['mrnas']:
                lines += self.make_feature_group(self.child_dict[mrna_id], is_gtf)
        # if no mrnas are present in a gene and there are CDS, make them here
        elif feature.get('cdss', False):
            for cds_id in feature['cdss']:
                lines += self.make_feature_group(self.child_dict[cds_id], is_gtf)
        # if this is a mrna with a child CDS, make it here
        elif feature.get('cds', False):
            lines += self.make_feature_group(self.child_dict[feature['cds']], is_gtf)

        return lines

    def make_feature(self, location, in_feature, is_gtf):
        """Make a single feature line for the file"""
        try:
            out_feature = {
                'seqname': location[0],
                'source': 'KBase',
                'type': in_feature.get('type', 'exon'),
                'start': str(get_start(location)),
                'end': str(get_end(location)),
                'score': '.',
                'strand': location[2],
                'frame': '0',
            }
            if is_gtf:
                out_feature['attribute'] = self.gen_gtf_attr(in_feature)
            else:
                out_feature['attribute'] = self.gen_gff_attr(in_feature)
        except Exception as e:
            traceback.print_exc()
            raise Exception('Unable to parse {}:{}'.format(in_feature, e))
        return out_feature

    def gen_gtf_attr(self, feature):
        """Makes the attribute line for a feature in gtf style"""
        if feature.get('type') == 'gene':
            return 'gene_id "{}"; transcript_id ""'.format(feature['id'])

        if "parent" in feature:
            feature['parent_gene'] = feature['parent']

        return 'gene_id "{}"; transcript_id "{}"'.format(
            feature.get('parent_gene', feature['id']),
            feature.get('parent_mrna', feature['id']))

    @staticmethod
    def gen_gff_attr(feature):
        """Makes the attribute line for a feature in gff style"""
        def _one_attr(key, val):
            return '{}={}'.format(key, urllib.quote(val, " /:"))

        # don't add an attribute that could be 0 without refactor
        for key in ('parent_gene', 'parent_mrna'):
            if key in feature:
                feature['parent'] = feature[key]
        attr_keys = (('id', 'ID'), ('parent', 'Parent'), ('note', 'note'))
        attrs = [_one_attr(pair[1], feature[pair[0]])
                 for pair in attr_keys if feature.get(pair[0])]
        attrs.extend([_one_attr('db_xref', '{}:{}'.format(*x))
                     for x in feature.get('db_xref', [])])
        attrs.extend([_one_attr(pair[0], pair[1])
                      for pair in feature.get('aliases', [''])
                      if isinstance(pair, list)])
        if feature.get('functional_descriptions'):
            attrs.append(_one_attr('function', ";".join(
                feature['functional_descriptions'])))
        if feature.get('functions'):
            attrs.append(_one_attr('product', ";".join(feature['functions'])))
        elif feature.get('function'):
            attrs.append(_one_attr('product', feature['function']))
        for ont in feature.get('ontology_terms', []):
            attrs.extend([_one_attr('Ontology_terms', x)
                          for x in feature['ontology_terms'][ont]])

        if 'inference_data' in feature:
            attrs.extend([_one_attr(
                'inference', ":".join([x[y] for y in ('category', 'type', 'evidence') if x[y]]))
                for x in feature['inference_data']])
        if 'trans_splicing' in feature.get('flags', []):
            attrs.append(_one_attr("exception", "trans-splicing"))
        return "; ".join(attrs)

    def get_common_location(self, location_array):
        """Merges a compound location array into an overall location"""
        contig = location_array[0][0]
        strand = location_array[0][2]
        min_pos = min([get_start(loc) for loc in location_array])
        max_pos = max([get_end(loc) for loc in location_array])
        common_length = max_pos - min_pos + 1
        common_start = min_pos if strand == '+' else max_pos
        return [contig, common_start, strand, common_length]

    @staticmethod
    def validate_params(params):
        if 'genome_ref' not in params:
            raise ValueError('required "genome_ref" field was not defined')
