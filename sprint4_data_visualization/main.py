import pandas as pd
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
print('\n', 10 * '=', 'DATASET', 10 * '=')

print('\n', 10 * '=', 'HEAD', 10 * '=')

print(df.head())

print('\n', 10 * '=', 'INFO', 10 * '=')

print(df.info())

print('\n', 10 * '=', 'DESCRIBE', 10 * '=')

print(df.describe())


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
print('\n', 10 * '=', 'CHROMOSSOME DATASET', 10 * '=')

print('\n', 10 * '=', 'HEAD', 10 * '=')

print(df_chromossomes.head())

print('\n', 10 * '=', 'INFO', 10 * '=')

print(df_chromossomes.info())

print('\n', 10 * '=', 'DESCRIBE', 10 * '=')

print(df_chromossomes.describe())

print('\n', 10 * '=', 'VALUE COUNTS', 10 * '=')

print(df_chromossomes.value_counts())


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

print('\n', 10 * '=', 'CHROMOSSOME[CHROMOSSOME] DATASET - CHROMOSSOME POSITION', 10 * '=')

print('\n', 10 * '=', 'HEAD', 10 * '=')

print(df_chromossomes['chromossome'].head(10))

print('\n', 10 * '=', 'DESCRIBE', 10 * '=')

print(df_chromossomes['chromossome'].describe())

print('\n', 10 * '=', 'VALUE COUNTS', 10 * '=')

print(df_chromossomes['chromossome'].value_counts())
