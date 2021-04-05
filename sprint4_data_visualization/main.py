import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from plots import distplot, jointplot
from plots import countplot
from plots import heatmap, clustermap
from chromossome import chromossome_sequence_to_dataframe


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
# status: DOING

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

print(group_chromossome_arm_1.count())

group_arm_1_chromossome = df_chromossomes[['chromossome', 'arm_1']].groupby('arm_1')

print(group_arm_1_chromossome.count())


group_chromossome_arm_2 = df_chromossomes[['chromossome', 'arm_2']].groupby('chromossome')

print(group_chromossome_arm_2.count())

group_arm_2_chromossome = df_chromossomes[['chromossome', 'arm_2']].groupby('arm_2')

print(group_arm_2_chromossome.count())


###############################################
# Data Visualization of Dataset
# Status: DOING
###############################################

#########################
# 01 - Distribution plots
#########################
distplot(df=df_chromossomes, column='chromossome')
distplot(df=df_chromossomes, column='arm_1')
sns.displot(df_chromossomes['arm_2'].dropna())

jointplot(df=df_chromossomes, x='chromossome', y='arm_1')
jointplot(df=df_chromossomes_arm_2, x='chromossome', y='arm_2')


########################
# 02 - Categorical plots
########################
countplot(df_chromossomes, x='arm_1', y='arm_2', hue='chromossome')
countplot(df_chromossomes, x='chromossome')
countplot(df_chromossomes, x='chromossome', hue='arm_1')
countplot(df_chromossomes, x='chromossome', hue='arm_2')
countplot(df_chromossomes, x='arm_1')
countplot(df_chromossomes, x='region_1')
countplot(df_chromossomes, x='band_1')
countplot(df_chromossomes, x='subband_1')
countplot(df_chromossomes, x='arm_2')
countplot(df_chromossomes, x='region_2')
countplot(df_chromossomes, x='band_2')
countplot(df_chromossomes, x='subband_2')

######################
# 03 - Matricial plots
######################
heatmap(df_chromossomes.groupby('chromossome').count())
heatmap(group_arm_1_chromossome.count())
heatmap(group_chromossome_arm_1.count())
heatmap(group_arm_2_chromossome.count())
heatmap(group_chromossome_arm_2.count())

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