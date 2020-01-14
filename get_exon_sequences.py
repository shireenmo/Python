import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

## This code will extract exon sequences into a fasta file based on feature coordinates

### indica
df_indica_coordinates =  pd.read_csv('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/processing_results/exon_coordinates_indica.txt',sep ='\t')
dicts = df_indica_coordinates.to_dict('records')#Convert the DataFrame to a dictionary.

ofile = open('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/processing_results/exon_seq_indica.fasta', "w") #Output
for record in SeqIO.parse('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/Databases/Oryza_indica.ASM465v1.dna.toplevel.fa','fasta'):
    def get_seq(start, end):# Extract sequence from start to end, take away one from start position because counting in python starts from zero
        sequence = record.seq[int(start) - 1:int(end)]
        return sequence
    for index in range(len(dicts)):

        if str(dicts[index]['chr']) == str(record.id) :  # If same Chromosome number
            get_seq(dicts[index]['start'], dicts[index]['end'])
            HEader = dicts[index]['id'] + '_' + str(dicts[index]['chr']) + '_' + str(dicts[index]['start']) + '_' + str(dicts[index]['end']) + '_' + dicts[index]['strand']  # Header includes ID_chr number_start_end_strand
            ofile.write(">" + str(HEader) + "\n" + str(get_seq(dicts[index]['start'], dicts[index]['end'])) + "\n")  # Write fasta of full exon sequences

ofile.close()

### N22
df_N22_coordinates =  pd.read_csv('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/processing_results/exon_coordinates_N22.txt',sep ='\t')
dicts = df_N22_coordinates.to_dict('records')#Convert the DataFrame to a dictionary.

ofile = open('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/processing_results/exon_seq_N22.fasta', "w") #Output
for record in SeqIO.parse('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/Databases/oryza_aus-toplevel-20180831.fa','fasta'):
    def get_seq(start, end):# Extract sequence from start to end, take away one from start position because counting in python starts from zero
        sequence = record.seq[int(start) - 1:int(end)]
        return sequence
    for index in range(len(dicts)):

        if str(dicts[index]['chr']) == str(record.id) :  # If same Chromosome number
            get_seq(dicts[index]['start'], dicts[index]['end'])
            HEader = dicts[index]['id'] + '_' + str(dicts[index]['chr']) + '_' + str(dicts[index]['start']) + '_' + str(dicts[index]['end']) + '_' + dicts[index]['strand']  # Header includes ID_chr number_start_end_strand
            ofile.write(">" + str(HEader) + "\n" + str(get_seq(dicts[index]['start'], dicts[index]['end'])) + "\n")  # Write fasta of full exon sequences
ofile.close()
