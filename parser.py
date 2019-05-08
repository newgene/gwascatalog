import os
import unicodedata
from collections import defaultdict

from csv import DictReader
from biothings.utils.dataload import dict_sweep, open_anyfile


def str2float(item):
    """Convert string type to float type
    """
    if item:
        return float(item)
    else:
        return item

def load_data(data_folder):

    input_file = os.path.join(data_folder,"gwas_catalog_v1.0.2-associations_e96_r2019-04-21.tsv")
    assert os.path.exists(input_file), "Can't find input file '%s'" % input_file
    with open_anyfile(input_file) as in_f:

        # Remove duplicated lines if any
        header = next(in_f).strip().split('\t')
        lines = set(list(in_f))
        reader = DictReader(lines, fieldnames=header, delimiter='\t')

        results = defaultdict(list)
        for row in reader:
            variant = {}

            # Use gDNA as variant identifier
            variant['_id'] = row['SNPS']
            variant['gwascatalog'] = {}
            variant['gwascatalog']['pubmed'] = int(row['PUBMEDID'])
            variant['gwascatalog']['date_added'] = row['DATE ADDED TO CATALOG']
            variant['gwascatalog']['study'] = row['STUDY']
            variant['gwascatalog']['trait'] = row['DISEASE/TRAIT']
            variant['gwascatalog']['region'] = row['REGION']
            variant['gwascatalog']['chrom'] = row['CHR_ID']
            variant['gwascatalog']['position'] = row['CHR_POS']
            variant['gwascatalog']['genes'] = row['REPORTED GENE(S)'].split(',')
            variant['gwascatalog']['context'] = row['CONTEXT']
            variant['gwascatalog']['intergenic'] = row['INTERGENIC']
            variant['gwascatalog']['risk_allele_freq'] = str2float(row['RISK ALLELE FREQUENCY'])
            variant['gwascatalog']['p_val'] = str2float(row['P-VALUE'])
            variant['gwascatalog']['p_val_mlog'] = str2float(row['PVALUE_MLOG'])
            variant['gwascatalog']['platform'] = row['PLATFORM [SNPS PASSING QC]']
            variant['gwascatalog']['study_accession'] = row['STUDY ACCESSION']
            variant['gwascatalog']['mapped_trait'] = row['MAPPED_TRAIT'].split(',')
            variant['gwascatalog']['mapped_trait_efo'] = row['MAPPED_TRAIT_URI'].split(',')
            results[variant['_id']].append(variant)

        # Merge duplications
        for v in results.values():
            if len(v) == 1:
                yield v[0]
            else:
                yield {
                    '_id': v[0]['_id'],
                    'gwascatalog': [i['gwascatalog'] for i in v]
                }