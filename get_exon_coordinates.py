from BCBio import GFF
import pandas as pd

"""Code will read the relevant portions of the indica GFF files to get exon coordinates 
of possible TFs from preliminary files manually created from BLAST for N22 and BGI"""
# Indica
in_file = "/Users/shiree/Desktop/rice-tf-db/get_exon_coordinates/Oryza_indica.ASM465v1.44.gff3"
in_handle = open(in_file)
id = []
chr = []
start = []
end = []
strand = []
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        for exon_feature in feature.sub_features:
            for _feature in exon_feature.sub_features:
                if _feature.type == 'exon':
                    split_id = str(feature.id).split(':')
                    id.append(split_id[1] + '-PA')# gene id
                    chr.append(rec.id)# Chromosome number
                    location = str(_feature.location).replace('[', ':').replace('](', ':').replace(')', ':')
                    split_loc = location.split(':')
                    start.append(int(split_loc[1]) + 1)# start pos
                    end.append(split_loc[2])  # end pos
                    strand.append(split_loc[3]) # strand

# Add all Indica exon_coordinates in a dataframe
df_indica = pd.DataFrame()
df_indica['id'] = id
df_indica['chr'] = chr
df_indica['start'] = start
df_indica['end'] = end
df_indica['strand'] = strand
in_handle.close()

df_Possible_indica = pd.read_csv('/Users/shiree/Desktop/rice-tf-db/input_data/map_MSU_to_BGIIndica.tsv',sep ='\t')#TF IDs of mapped BGI Indica
df_Possible_indica = df_Possible_indica.drop_duplicates(subset= df_Possible_indica.columns[1], keep="last")
# take exon_coordinates of possible TFs
df_indica = df_indica.loc[df_indica.iloc[:,0].isin(df_Possible_indica.iloc[:,1])]
df_indica.to_csv('/Users/shiree/Desktop/rice-tf-db/get_exon_coordinates/exon_coordinates_indica.txt', sep='\t', index=False, header=True)

#N22
in_file = "/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/Databases/oryza_aus_core_3_87_1.gff"
in_handle = open(in_file)
id = []
chr = []
start = []
end = []
strand = []
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        for exon_feature in feature.sub_features:
            for _feature in exon_feature.sub_features:
                if _feature.type == 'exon':
                    id.append(exon_feature.id)# gene id
                    chr.append(rec.id)# Chromosome number
                    location = str(_feature.location).replace('[', ':').replace('](', ':').replace(')', ':')
                    split_loc = location.split(':')
                    start.append(int(split_loc[1]) + 1)# start pos
                    end.append(split_loc[2])  # end pos
                    strand.append(split_loc[3]) # strand

# Add all N22 exon coordinates in a dataframe
df_N22 = pd.DataFrame()
df_N22['id'] = id
df_N22['chr'] = chr
df_N22['start'] = start
df_N22['end'] = end
df_N22['strand'] = strand
in_handle.close()

df_Possible_N22 = pd.read_csv('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/processing_results/Possible_N22_tfs_but_some_will_not_be.txt',sep ='\t')#TF IDs of mapped N22
df_Possible_N22 = df_Possible_N22.drop_duplicates(subset= df_Possible_N22.columns[0], keep="last")
# Take exon_coordinates of possible TFs
df_N22 = df_N22.loc[df_N22.iloc[:,0].isin(df_Possible_N22.iloc[:,0])]
df_N22.to_csv('/Users/shiree/NewDropBox/Dropbox/liv/PycharmProjects/rice-tf-db/processing_results/exon_coordinates_N22.txt', sep='\t', index=False, header=True)

