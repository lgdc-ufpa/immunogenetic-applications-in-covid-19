import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from plots import *
from chromossome import chromossome_sequence_to_dataframe
from chromossome import chromossome_dataframe_order
from styles import plt_figure_size
from tqdm import tqdm
from alleles import *


##############
# Load Dataset
# status: DONE
##############
df = pd.read_csv('input/gene_link-orphanet_link-ensembl_id-ensembl_chromossome.csv', sep='>')

#################
# Analyse Dataset
# status: DONE

# There are 308 genes in 308 different chromossome positions
# Each gene has its unique ensembl code
#################
# print('\n', 10 * '=', 'DATASET', 10 * '=')

# print('\n', 10 * '=', 'HEAD', 10 * '=')

# print(df.head())

# print('\n', 10 * '=', 'INFO', 10 * '=')

# print(df.info())

# print('\n', 10 * '=', 'DESCRIBE', 10 * '=')

# print(df.describe())


################
# Question 01
# The genes are next to each other consideringchromossome position?
# status: DONE

# https://medlineplus.gov/genetics/understanding/howgeneswork/genelocation/
################

###################################
# Extract Chromossomes from Dataset
# Status: DONE
###################################
chromossomes = df['Localização Cromossômica']


####################################
# Generate Dataframe of Chromossomes
# Status: DONE
####################################
df_chromossomes = chromossome_sequence_to_dataframe(chromossomes)
df_chromossomes_ordered = chromossome_dataframe_order(df_chromossomes, x=True)
print(df_chromossomes_ordered)

############################
# Anayse df_chromossome Dataset
# Status: DONE

# May have correlaction between chromossome
# Chromossomes 1, 2, 11, 5, 17, 6, and 19 compounds 50.82 % of the 291
# Chromossomes 1, 2, 11, 5, 17, 6, 19, and 16 compounds 52.88 % of the 308

# Maybe can be correlaction between chromossomes
# 1 and 2 (strong),
# 5 and 6 (strong),
# 17 and 19 (not too strong)

# 11 may have not correlaction with either chromossomes
############################
# print('\n', 10 * '=', 'CHROMOSSOME DATASET', 10 * '=')

# print('\n', 10 * '=', 'HEAD', 10 * '=')

# print(df_chromossomes.head())

# print('\n', 10 * '=', 'INFO', 10 * '=')

# print(df_chromossomes.info())

# print('\n', 10 * '=', 'DESCRIBE', 10 * '=')

# print(df_chromossomes.describe())

# print('\n', 10 * '=', 'VALUE COUNTS', 10 * '=')

# print(df_chromossomes.value_counts())


############################
# Anayse df_chromossome['chromossome'] Dataset
# Status: DONE

# 291 genes of a total of 308 have chromossome position specified

# The 291 genes appears in all 22 chromossomes

# The X chromossome appears 17 times: 	  17/291 == 05.84% | 17/308 == 5.51%
# The X chromossome is lest frequent then chromossome 1, 2, 11, 5, and 17

# The chromossome 1 is the most frequent: 34/291 == 11.68% | 34/308 == 11,03%
# The chromossome 2 succeeds 1: 		  23/291 == 07.90% | 23/308 == 07.46%
# The chromossome 11   ''    2: 		  22/291 == 07.56% | 22/308 == 07.14%
# The chromossome 5    ''   11:			  19/291 == 06.52% | 19/308 == 06.16%
# The chromossome 17   ''    5:			  18/291 == 06.18% | 18/308 == 05.84%
# The chromossome 6    ''   17:			  16/291 == 05.49% | 16/308 == 05.19%
# The chromossome 19   ''    6:			  16/291 == 05.49% | 16/308 == 05.19%
# The chromossome 16   ''	19:			  15/291 == 05.15  | 15/308 == 04.87%

# Chromossomes 1, 2, 11, 5, 17, 6, and 19 compounds 50.82 % of the 291
# Chromossomes 1, 2, 11, 5, 17, 6, 19, and 16 compounds 52.88 % of the 308

# Maybe can be correlaction between chromossomes
# 1 and 2 (strong),
# 5 and 6 (strong),
# 17 and 19 (not too strong)

# 11 may have not correlaction with either chromossomes
############################

# print('\n', 10 * '=', 'CHROMOSSOME[CHROMOSSOME] DATASET - CHROMOSSOME POSITION', 10 * '=')

# print('\n', 10 * '=', 'HEAD', 10 * '=')

# print(df_chromossomes['chromossome'].head(10))

# print('\n', 10 * '=', 'DESCRIBE', 10 * '=')

# print(df_chromossomes['chromossome'].describe())

# print('\n', 10 * '=', 'VALUE COUNTS', 10 * '=')

# print(df_chromossomes['chromossome'].value_counts())

bol1 = df_chromossomes['arm_2'] == 'q'
bol2 = df_chromossomes['arm_2'] == 'p'
condition = (bol1 | bol2)
df_chromossomes_arm_2 = df_chromossomes[condition][['chromossome', 'arm_2']]
group_chromossome_arm_1 = df_chromossomes[['chromossome', 'arm_1']].groupby('chromossome')

# print(group_chromossome_arm_1.count())

group_arm_1_chromossome = df_chromossomes[['chromossome', 'arm_1']].groupby('arm_1')

# print(group_arm_1_chromossome.count())


group_chromossome_arm_2 = df_chromossomes[['chromossome', 'arm_2']].groupby('chromossome')

# print(group_chromossome_arm_2.count())

group_arm_2_chromossome = df_chromossomes[['chromossome', 'arm_2']].groupby('arm_2')

# print(group_arm_2_chromossome.count())


###############################################
# Data Visualization of Dataset
# Status: DOING
###############################################

#########################
# 01 - Distribution plots
#########################

# ok
# distplot(df=df_chromossomes_ordered, column='chromossome')
# distplot(df=df_chromossomes, column='arm_1')
# sns.displot(df_chromossomes['arm_2'].dropna())
# distplot(df_chromossomes[['arm_1', 'region_1']].sort_values(by='region_1'), column='region_1')
# distplot(df_chromossomes[['arm_1', 'region_1']].sort_values(by='region_1'), column='arm_1')

# ugly
# jointplot(df=df_chromossomes, x='chromossome', y='arm_1')
# jointplot(df=df_chromossomes_arm_2, x='chromossome', y='arm_2')


########################
# 02 - Categorical plots
########################
# ok
# countplot(df_chromossomes_ordered, x='chromossome', color='gray')

# issue: order all dataframe, and not only the chromossome column
# countplot(df_chromossomes_ordered, x='chromossome', hue='arm_1')

# countplot(df_chromossomes, x='chromossome', hue='arm_1')
# countplot(df_chromossomes, x='chromossome', hue='arm_2')
# countplot(df_chromossomes, x='arm_1')
# countplot(df_chromossomes.sort_values(by='region_1'), x='region_1')
# countplot(df_chromossomes.sort_values(by='band_1'), x='band_1')
# countplot(df_chromossomes.sort_values(by='subband_1'), x='subband_1')
# countplot(df_chromossomes, x='arm_2')
# countplot(df_chromossomes.sort_values(by='region_2'), x='region_2')
# countplot(df_chromossomes.sort_values(by='band_2'), x='band_2')
# countplot(df_chromossomes.sort_values(by='subband_2'), x='subband_2')

######################
# 03 - Matricial plots
######################

# Issue: This graph will be very good if all dataframe
# were be ordered by chromossomes
# heatmap(df_chromossomes.groupby('chromossome').count())

# bad
# heatmap(group_arm_1_chromossome.count())

# Issue: order by chromossomes
# heatmap(group_chromossome_arm_1.count())

# bad
# heatmap(group_arm_2_chromossome.count())

# Issue: ordered by chromossomes
# heatmap(group_chromossome_arm_2.count())

# Issue: This graph will be very good if all dataframe
# were be ordered by chromossomes
# clustermap(df_chromossomes.groupby('chromossome').count(), annot=True)

#######################
# 04 - Regression plots
# Don't make sense
#######################


#######################
# 05 - pair grids
# Needs non categorical column values
#######################
plt.show()


##########################
# Question 02
# What genes are important?
# Status: DONE

# Analysing just the gene incidences in chromossomes,
# the 34 genes of chromossome 01 are important and should be
# analyzed and, if some pattern is observed, try to undestand it
# and verify if this patter occurs in other genes of chromossomes who has
# high frequency too
##########################


######################################
# Data analysis
# df3 = chromossome and gene dataframe
# 		ordered by chromossome
# Status: DONE
######################################

df2 = pd.DataFrame()

df2['gene'] = df['Código no Ensembl']

df2['chromossome'] = df_chromossomes['chromossome']

df3 = pd.DataFrame()

li_chrm = ['1', '2', '5', '6', '11', '17', '19']

# order by chromossome
for chrm in li_chrm:
	bol_chrm = df2['chromossome'].apply(lambda x: x == chrm)
	df3 = pd.concat([df3, df2[df2[['chromossome']].apply(lambda x: x == chrm).dropna()['chromossome']]], axis=0)


#######################################
# Question 03
# There is something in common in genes
# present in chromossome 01?

# Status: DOING
#######################################


#########################
# Data Anaylis
# Status: DONE

# Genes of chromossome 01
#########################

df4 = df3[df3['chromossome'] == '1'].reset_index()[['gene']]

##########################
# Data Anaylis
# Status: DOING

# Get genes (only) informations

# There is only 266 of a total of 308 with patogenic variants
##########################

##########
# Read Dataset
# Genes with patogenic variants
# Status: DONE
##########
df5 = pd.DataFrame()

relative_path_genes = os.path.join("input", "table_with_genes_with_patogenic_variants")

path_genes_with_patogenic_variants = os.path.join(os.getcwd(), relative_path_genes)

for gene in os.listdir(path_genes_with_patogenic_variants):
	path_gene = os.path.join(path_genes_with_patogenic_variants, gene)
	df_temp = pd.read_csv(path_gene, sep=',')
	df_temp['gene'] = gene.split('.')[0]
	df5 = pd.concat([df5, df_temp], axis=0)

del df5['Unnamed: 0']


#############################
# Number of variants per gene
# Status: DONE
#############################

df6 = df5[['gene', 'Variant ID']]

df7 = df6.groupby('gene').count().sort_values(by='Variant ID', ascending=False).reset_index()

df7.columns = ['gene', 'variants']

bol_patogenic_genes_chr_1 = [True if gene in np.array(df4['gene']) else False for gene in df7['gene']]

df8 = df7[bol_patogenic_genes_chr_1].reset_index()[['gene', 'variants']]

n_genes_patogenic_chr_1 = len([gene for gene in df7['gene'] if gene in np.array(df4['gene'])])

n_genes_chr_1 = len(df4['gene'])

print(f"{n_genes_patogenic_chr_1} ({100 * (n_genes_patogenic_chr_1 / n_genes_chr_1):.2f}%) \
	genes (of {n_genes_chr_1}) in chromossome 01 has patogenic variants")

df9 = df8.iloc[0:3]

df10 = df8.iloc[3:]

df7.describe()
df8.describe()
df10.describe()

########################
# Data Visualization of number of patogenic genes per gene in chr 1
# Status: DOING
########################

#########################
# 01 - Distribution plots
#########################
# distplot(df7, column='variants')
# distplot(df8, column='variants')
# distplot(df10, column='variants')

# kdeplot(df7, column='variants')
# kdeplot(df8, column='variants')
# kdeplot(df10, column='variants')

########################
# 02 - Categorical plots
########################
# boxplot(df7)
# boxplot(df8)
# boxplot(df10)

#######
# Some more data analysis
# Status: DOING
#######

n_variants_total = len(df6.groupby('Variant ID').count().reset_index()['Variant ID'].unique())

n_genes_with_variants = len(df6['gene'].unique())

print(df6['Variant ID'].value_counts())

print(f"There is a total of {n_variants_total} unique variants in all {n_genes_with_variants} genes with variants")

df_nGenes_x_variant = df6.groupby('Variant ID').count().reset_index().sort_values(by='gene', ascending=False)

variant_most_frequency_in_genes = df_nGenes_x_variant.iloc[0, 0]

genes_from_most_frequent_variant = df6[df6['Variant ID'] == variant_most_frequency_in_genes]['gene'].unique()

print("Analysing the genes from the most frequent variant")

print(pd.DataFrame(genes_from_most_frequent_variant).describe())

print(f"{variant_most_frequency_in_genes} is the variant most frequency in genes")

print(f"{len(genes_from_most_frequent_variant)} is/are the number of genes from the most frequent variant")

print(f"{genes_from_most_frequent_variant} is/are the genes from the most frequent variant")

print("Analysing the df_nGenes_x_variant")

print(df_nGenes_x_variant.describe())

# distplot(df_nGenes_x_variant, column='gene', bins=10)

# barplot(df_nGenes_x_variant, 'Variant ID', 'gene')


############
# Question 04
# What is the chromossome of the most frequent variant?
############

gene_from_most_frequent_variant = df5[df5['Variant ID'] == variant_most_frequency_in_genes]['gene'].iloc[0]

chromossome_from_most_frequent_variant = df2[df2['gene'] == gene_from_most_frequent_variant].iloc[0, 1]

print(f"Chromossome {chromossome_from_most_frequent_variant} contains the gene {gene_from_most_frequent_variant} wich contains the most frequent variant {variant_most_frequency_in_genes}")

###########
# Question 05
# What are the chromossomes wich contains the most 25 % frequent variants of the genes?
###########

percent_25 = int(len(df_nGenes_x_variant) / 4)

df_nGenes_x_variant.iloc[0:percent_25]

top_variants = (np.array(df_nGenes_x_variant.iloc[0:percent_25]['Variant ID'].unique()))

df_gen_var_nVar = pd.merge(left=df6[['gene', 'Variant ID']], right=df_nGenes_x_variant.iloc[0:percent_25], how='inner', on='Variant ID')

df_gen_var_nVar.columns = ['gene', 'Variant ID', 'n_var_in_genes']

df_gen_var_nVar.describe()

df_final_01 = pd.merge(left=df2, right=df_gen_var_nVar, how='inner', on='gene').sort_values(by='n_var_in_genes', ascending=False)

del df_final_01['n_var_in_genes']

# df_final_01.drop_duplicates()
# df_final_01.drop_duplicates(subset=['gene'])
# df_final_01.drop_duplicates(subset=['gene']) # Top variants per gene

#############
# Data Analysis
# Analyse df_final_01
# Status: DOING
#############

# df_chr_variant = chromossome_dataframe_order_with_all_columns(df_final_01, x=True)

n_variants_in_general = df_final_01.drop_duplicates(subset=['Variant ID']).shape[0]

dff2 = df_final_01.drop_duplicates(subset=['Variant ID'])[['chromossome', 'gene', 'Variant ID']]

dff2.columns = ['chromossome', 'gene', 'variant']

print(f"There is {n_variants_in_general} variants in whole chromossomes and genes")

# countplot(dff2, x='chromossome')

li = [str(i) for i in range(1, 23)]

li.append('X')

dff3 = pd.DataFrame()

for chrm in li:
	bol_temp = dff2['chromossome'] == chrm
	dff3 = pd.concat([dff3, dff2[bol_temp]])

# dff3.describe()

chrm_unique, chrm_top, chrm_freq = dff3.describe()['chromossome'].iloc[1:]
gene_unique, gene_top, gene_freq = dff3.describe()['gene'].iloc[1:]
vari_unique, vari_top, vari_freq = dff3.describe()['variant'].iloc[1:]

print(">>Chromossome")
print(f"There is/are {chrm_unique} unique chromossomes")
print(f"The most frequent chromossome {chrm_top} appears {chrm_freq} times")

print(">>Gene")
print(f"There is/are {gene_unique} unique genes")
print(f"The most frequent gene {gene_top} appears {gene_freq} times")

print(">>Variant")
print(f"There is/are {vari_unique} unique variants")
print(f"The most frequent variant {vari_top} appears {vari_freq} times")

# countplot(dff3, x='chromossome')
#################
# Question 06
# Which are the patogenic variants?
# Status: DONE
#################

clinic_signal_uniqued = df5['Clin. Sig.'].unique()

print(clinic_signal_uniqued)

df5[df5['Clin. Sig.'] == 'pathogenic']['Variant ID']

dff4 = df5[df5['Clin. Sig.'] == 'pathogenic'].drop_duplicates(subset=['Variant ID'])

print('\n>>>Pathogenic Variants')
print(dff4[['Variant ID']])

# countplot(dff3, x='chromossome')

###############
# Question 07
# What are the genes that the pathogenic variants are present
# Status: DONE
###############

c = list(dff4.columns)

c[0] = 'variant'

dff4.columns = c

dff5 = pd.merge(left=dff3, right=dff4, how='inner', on='variant').iloc[:, :3]

pathogenic_genes = dff4[['gene']].drop_duplicates()

n_pathogenic_genes = pathogenic_genes.shape[0]

print(f"There are a total of {n_pathogenic_genes} genes with pathogenic variants")

# countplot(dff3, x='chromossome')

############
# Data Visualization
# Chromossome x number of variants
# Status: DONE
############

dff6 = dff3.groupby(by='chromossome').count()[['variant']].reset_index().sort_values(by='variant', ascending=False)
dff6.columns = ['chromossome', 'n_variants']
# barplot(dff6, x='chromossome', y='variant')

################
# Data Visualization
# Chromossome x number of genes
# Status: DONE
################

dff7 = dff3.drop_duplicates(subset=['gene'])[['chromossome']]

# countplot(df=dff7, x='chromossome')

################
# Data Analysis
# Alele informations
# Status: DONE
################

print("Alele informations")

df5[['Alleles']].describe()

dff8 = df5[['gene', 'Variant ID', 'Alleles']].drop_duplicates()

dff8.describe()

n_alleles = dff8['Alleles'].unique().shape[0]

print(f"There is a total of {n_alleles}")

# distplot(dff8[['Alleles']], column='Alleles')

dff9 = df5[['Location', 'gene', 'Variant ID', 'Alleles']].drop_duplicates()

dff9.describe()

dff10 = df5[df5['meta_lr_class'] == 'damaging']

print('Damaging Variants')

dff10[['Variant ID']]

dff11 = dff10[['Location', 'gene', 'Variant ID', 'Alleles', 'Clin. Sig.', 'MetaLR']]
# dff11 = df5[['Location', 'gene', 'Variant ID', 'Alleles', 'Clin. Sig.', 'MetaLR']]

li = [x for x in range(1, 23)]

li.append('X')

dff12 = pd.DataFrame()

for chrm in li:
	b = dff11['Location'].apply(lambda x: str(chrm) == x.split(':')[0])
	list_dataframe = [dff12, dff11[b]]
	dff12 = pd.concat(list_dataframe, axis=0)

dff12[['Location']] = pd.DataFrame(dff12['Location'].apply(lambda x: x.split(':')[0]))

dff12[['MetaLR']] = pd.DataFrame(dff12['MetaLR'].apply(lambda x: float(x)))
# dff12[['MetaLR']] = dff12[['MetaLR']].apply(lambda x: pd.to_numeric(x, errors='coerce'))


dff12.describe()

# countplot(dff12, x='Alleles')

# countplot(dff12, x='Location')

plt_figure_size(20, 8)

boxplot(dff12, x='Alleles', y='MetaLR')

# plt_figure_size(20, 8)

# boxplot(dff12, x='Location', y='MetaLR')

# plt_figure_size(20, 8)

# barplot(dff12, x='Alleles', y='MetaLR')

# plt_figure_size(20, 8)

# barplot(dff12, x='Location', y='MetaLR')

dff13 = pd.pivot_table(dff12, values='MetaLR', index=['Location'], columns=['Alleles'])

plt_figure_size(20, 8)

heatmap(dff13)

# plt_figure_size(20, 8)

##########################
# Data Anaylis
# Status: TODO

# Get genes informations about populations
##########################

relative_path_populations = os.path.join("input", "table_with_frequency_of_patogenic_variants_in_populations")

path_genes_with_populations = os.path.join(os.getcwd(), relative_path_populations)

dfpop = pd.DataFrame()

for pop in tqdm(os.listdir(path_genes_with_populations)):
	path_pop = os.path.join(relative_path_populations, pop)
	dftemp = pd.read_csv(path_pop)
	gene = pop.split('_')[1]
	variant = pop.split('_')[-1].split('.')[0]
	dftemp['gene'] = gene
	dftemp['variant'] = variant
	del dftemp['Unnamed: 1']
	dfpop = pd.concat([dfpop, dftemp])


###########
# Data Analysis
# Split alleles frequency and count
# Status: DONE
###########
# https://stackoverflow.com/questions/8569201/get-the-string-within-brackets-in-python

# TODO IN FUTURE
col_allele = 'Allele: frequency (count)'
alleles_colection = dfpop[col_allele].apply(alleles_split)
df_alleles_colection = pd.DataFrame(alleles_colection)

for alleles in alleles_colection:
	print(alleles)
	for allele in alleles:
		print(allele)
		a, freq, count = allele
		print(a, freq, count)
	# print(list(map(lambda x: len(x) == 1, alleles)))
	break

for line in dfpop.itertuples(index=False):
	print(line.Population)	
	break

dfpop3 = pd.read_csv(os.path.join('input', 'output_3.csv'), '>')
dfpopdata = pd.read_csv(os.path.join('input', 'data.csv'))


##############
# Data Visualization including populations
# What is the frequency of the variants in the populations
# Status: DOING
##############

dfpopdataexome = dfpopdata[dfpopdata['source'] == 'gnomAD exomes']
dfpopdataexomeall = dfpopdataexome[dfpopdataexome['population'] == 'all']
dfpopdataexomenotall = dfpopdataexome[dfpopdataexome['population'] != 'all']
dfpopdataexomenotallpathogenic = dfpopdataexomenotall['variant']

dfpopdatagenome = dfpopdata[dfpopdata['source'] == 'gnomAD genomes r3.0']

dfpopdata3 = pd.pivot_table(data=dfpopdata, values='frequency', index='variant', columns='population')
dfpopdata4 = pd.pivot_table(data=dfpopdata, values='frequency', index=['source', 'gene', 'variant', 'allele'], columns=['population'])

distplot(dfpopdata, column='frequency')
distplot(dfpopdataexome, column='count')