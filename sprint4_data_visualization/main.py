import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from plots import distplot, jointplot
from plots import countplot
from plots import heatmap, clustermap
from chromossome import chromossome_sequence_to_dataframe
from chromossome import chromossome_dataframe_order

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
distplot(df=df_chromossomes_ordered, column='chromossome')
distplot(df=df_chromossomes, column='arm_1')
sns.displot(df_chromossomes['arm_2'].dropna())
distplot(df_chromossomes[['arm_1', 'region_1']].sort_values(by='region_1'), column='region_1')
distplot(df_chromossomes[['arm_1', 'region_1']].sort_values(by='region_1'), column='arm_1')

# ugly
jointplot(df=df_chromossomes, x='chromossome', y='arm_1')
jointplot(df=df_chromossomes_arm_2, x='chromossome', y='arm_2')


########################
# 02 - Categorical plots
########################
# ok
countplot(df_chromossomes_ordered, x='chromossome', color='gray')

# issue: order all dataframe, and not only the chromossome column
# countplot(df_chromossomes_ordered, x='chromossome', hue='arm_1')

countplot(df_chromossomes, x='chromossome', hue='arm_1')
countplot(df_chromossomes, x='chromossome', hue='arm_2')
countplot(df_chromossomes, x='arm_1')
countplot(df_chromossomes.sort_values(by='region_1'), x='region_1')
countplot(df_chromossomes.sort_values(by='band_1'), x='band_1')
countplot(df_chromossomes.sort_values(by='subband_1'), x='subband_1')
countplot(df_chromossomes, x='arm_2')
countplot(df_chromossomes.sort_values(by='region_2'), x='region_2')
countplot(df_chromossomes.sort_values(by='band_2'), x='band_2')
countplot(df_chromossomes.sort_values(by='subband_2'), x='subband_2')

######################
# 03 - Matricial plots
######################

# Issue: This graph will be very good if all dataframe
# were be ordered by chromossomes
heatmap(df_chromossomes.groupby('chromossome').count())

# bad
heatmap(group_arm_1_chromossome.count())

# Issue: order by chromossomes
heatmap(group_chromossome_arm_1.count())

# bad
heatmap(group_arm_2_chromossome.count())

# Issue: ordered by chromossomes
heatmap(group_chromossome_arm_2.count())

# Issue: This graph will be very good if all dataframe
# were be ordered by chromossomes
clustermap(df_chromossomes.groupby('chromossome').count(), annot=True)

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
# Status: DOING
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
