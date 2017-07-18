#!/usr/bin/env python3

import sys
import os
import csv
import warnings
import logging
from pprint import pprint

from pkg_resources import resource_string
#foo_config = resource_string(__name__, 'foo.conf')

import pandas as pd
from Bio import SeqIO

aa_set = set('GPAVLIMCFYWHKRQNEDST')  # aminoacid one-letter code
dn_dir = os.path.dirname(__file__)
db_dir = os.path.abspath(os.path.join(dn_dir, 'db'))
template_dir = os.path.abspath(os.path.join(dn_dir, 'templates'))


def parse_drm():
    '''Parse drug resistance mutations listed in files db/*Variation.txt'''

    df_list = []
    for gene, drm_file_name in [('protease', 'masterComments_PI.txt'),
                                ('RT', 'masterComments_RTI.txt'),
                                ('integrase', 'masterComments_INI.txt')]:
        drm_file = os.path.join(db_dir, drm_file_name)
        try:
            d1 = pd.read_table(drm_file, header=0, names=['pos', 'mut', 'category',
                               'commented', 'comment'])
        except:  # integrase does not have commented column
            d1 = pd.read_table(drm_file, header=0, names=['pos', 'mut', 'category',
                               'comment'])
        gs = pd.Series([gene] * len(d1))
        d1['gene'] = gs
        df_list.append(d1)
    df = pd.concat(df_list)
    return df


def write_header(handle, subtype_file=None):
    '''Write header to a file in markdown format'''
    md_header = 'Drug resistance mutations detected by NGS sequencing'
    mc = len(md_header)
    md_header += '\n' + '=' * len(md_header) + '\n\n'
    md_header +='Subtype inference with blast\n'
    md_header +='----------------------------\n'
    md_header += '|     subtype     | support [%] |\n'
    md_header += '|:{:-^15}:|:{:-^11}:|\n'.format('', '', '')
    if subtype_file:
        with open(subtype_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for mtype, freq in spamreader:
                int_freq = int(round(100 * float(freq), 0))
                if int_freq >= 1:
                    md_header += '|{: ^17}|{: ^13}|\n'.format(mtype, int_freq)
    md_header += '''

Parsing mutations
-----------------

The list of mutations was downloaded from HIVdb and includes:

- xyz positions on protease
- zyx positions on RT
- abc positions on integrase.
'''
    print(md_header, file=handle)


def aa_unpack(mut_string):
    if not mut_string.startswith('NOT'):
        return set(mut_string)
    else:
        return aa_set - set(mut_string.split()[1])


def parse_merged(mer_file):
    '''This is done by hand because it was too complicated to achieve
    this functionality with panda alone'''

    with open(mer_file) as csvfile:
        reader = csv.DictReader(csvfile)
        mdf = pd.DataFrame(columns=reader.fieldnames)
        for row in reader:
            row['freq'] = float(row['freq'])
            row['pos'] = int(row['pos'])
            if row['mut_x'] in aa_unpack(row['mut_y']):
                mdf = mdf.append(row, ignore_index=True)
            elif row['mut_x'] == '-' and 'd' in aa_unpack(row['mut_y']):
                mdf = mdf.append(row, ignore_index=True)
            elif row['mut_x'] and row['mut_y'] == '':
                row['category'] = 'unannotated'
                mdf = mdf.append(row, ignore_index=True)
    return mdf


def main(mut_file='annotated_mutations.csv', subtype_file='subtype_evidence.csv'):
    '''What does the main do?'''

    import subprocess

    rh = open('report.md', 'w')
    write_header(rh, subtype_file)

    #cov_info = get_coverage_info()

    logging.info('Parsing DRM from database')
    resistance_mutations = parse_drm()
    logging.debug('Shape is: %s', str(resistance_mutations.shape))

    logging.info('Reading mutations from %s' % mut_file)
    mutation_detected = pd.read_csv(mut_file)
    logging.debug('Shape is: %s' % str(mutation_detected.shape))

    mpd = pd.merge(mutation_detected, resistance_mutations, how='left',
                   on=['gene', 'pos'])
    mpd['pos'] = mpd['pos'].astype(int)
    mpd.to_csv(path_or_buf='merged_muts_drm_annotated.csv', index=False)
    logging.debug('Shape of raw merged is: %s' % str(mpd.shape))

    # too complicated with panda, do it by hand
    drms = parse_merged('merged_muts_drm_annotated.csv')
    logging.debug('Shape of merged is: %s' % str(drms.shape))
    #os.remove('merged_muts.csv')

    drms.drop(['commented', 'mut_y'], axis=1, inplace=True)
    drms.rename(columns={'mut_x': 'mut'}, inplace=True)

    cols = ['gene', 'pos', 'mut', 'freq', 'category']
    drms = drms[cols]
    drms.sort_values(cols[:3], inplace=True)
    drms.to_csv('annotated_DRM.csv', index=False)

    for gene in ['protease', 'RT', 'integrase']:
        gene_muts = drms[drms.gene == gene]
        if gene_muts.shape[0] == 0:
            logging.info('No mutations on %s' % gene)
            print('No mutations on %s' % gene, file=rh)
            print('-' * (len(gene) + 16), file=rh)
            print(file=rh)
            continue

        grouped = gene_muts.groupby(['pos', 'mut'])
        print('%s' % gene, file=rh)
        print('-'*len(gene), file=rh)
        print('| position | mutation | frequency [%] |      category      |',
              file=rh)
        print('|:{:-^8}:|:{:-^8}:|:{:-^13}:|:{:-^18}:|'.format('', '', '', ''),
              file=rh)
        for name, group in grouped:
            # same pos mut tuple must give same annotation, probably redundant
            #assert group['category'].nunique() == 1, group['category'].describe()
            mut_cat = group['category'].unique()[0]
            int_freq = int(round(100 * group['freq'].sum(), 0))
            print('|{: ^10}|{: ^10}|{: ^15}|{: ^20}|'.format(int(name[0]),
                                                     name[1],
                                                     int_freq,
                                                     mut_cat), file=rh)
        print('\n', file=rh)
    rh.close()
    # convert to PDF with pandoc
    tmpl_file = os.path.abspath(os.path.join(template_dir, 'template.tex'))
    os.symlink(tmpl_file, 'template.tex')
    pand_cml = 'pandoc --template=template.tex report.md -o report.pdf'
    logging.info('Run %s' % pand_cml)
    subprocess.call(pand_cml, shell=True)


if __name__ == '__main__':
    main()
