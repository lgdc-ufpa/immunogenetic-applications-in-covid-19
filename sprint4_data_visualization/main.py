import os
import numpy as np
import pandas as pd
from typing import List
from tqdm import tqdm
from modules import split_alelles_from_frequencies
"""
Sprint 4.1: analyze tables with genes
with patogenic variant in populations
"""



# Get output folders of sprint 03
output_directories_of_sprint3: list = \
os.listdir('input')

# Get folder of gene variants
directory_genes: str = [d for d in \
	output_directories_of_sprint3 \
	if 'genes' in d][0]

# Get folder of frequency of pathogenic
# variants in population
directory_variants: str = [d for d in \
output_directories_of_sprint3 \
if 'freq' in d][0]

# form the abs_path_variants folder
abs_path_genes = os.path.join(\
	os.getcwd(), 'input', directory_genes)

# form the abs_path_variants folder
abs_path_variants = os.path.join(\
	os.getcwd(), 'input', directory_variants)

# Get the csv files of variants
li_csv_genes = os.listdir(abs_path_genes)

# Get the csv files of freq
li_csv_variants = os.listdir(abs_path_variants)
# print(li_csv_variants.__len__())

# Create the output dataframe
index = [i for i in range(len(li_csv_variants))]
columns = ['gene', 'variant']
df_output = pd.DataFrame(index=index, columns=columns)

# Create output csv file
with open('output_3.csv', 'w') as f:
	line = 'gene>variant>population>alelle>frequency>count\n'
	f.write(line)

# Iterate over csv files
for csv_gene_file in tqdm(li_csv_genes):	

	# form the absolute path to the
	# current csv gene file
	abs_path_to_current_csv_gene = os.path.join(\
		abs_path_genes, csv_gene_file)

	# read the current csv file
	# examples: ENSG00000159189.csv
	df = pd.read_csv(abs_path_to_current_csv_gene)

	# get the variants	
	try:
		variants: pd.core.series.Series = \
		df['Variant ID']
		# print(variants)
	except Exception as e:
		print(f'>>> {abs_path_to_current_csv_gene}')
		raise e
	else:
		pass
	finally:
		pass
	

	# get the id_ensembl of the gene
	id_ensembl: str = csv_gene_file.replace('.csv', '')
	# print(id_ensembl)

	# get the frequency of each variant in each
	# population already indexed in output df
	for variant in variants:


		# Get the list of csv variant files
		# which matches with the current variant
		# print('variant: ', variant)		
		li_variant_csv: list = [v \
		for v in li_csv_variants \
		if variant in v and not '~lock' in v]
		# print(li_variant_csv)

		# Iterate over the csv files wich matched		
		for v_csv in li_variant_csv:
			# print('variant csv file: ', v_csv)

			# form the absolute path to the
			# current csv variant file
			abs_path_to_current_csv_variant = os.path.\
			join(abs_path_variants, v_csv)

			# read the current variant file
			# examples:gnomadexomes_ENSG00000 (...)
			# 159189_rs200206736.csv
			df_variants = pd.read_csv(\
				abs_path_to_current_csv_variant)

			# print(df_variants)

			# Get the populations
			populations: pd.core.series.Series = \
			df_variants['Population']
			# print('population[0]: ', populations[0],\
				# df_variants['Population'][0])

			# print(populations)

			# Get the frequencies in populations
			frequencies: pd.core.series.Series = \
			df_variants['Allele: frequency (count)']

			# print('OK', csv_gene_file, abs_path_to_current_csv_variant)

			# Iterate over the populations
			# where the current variant is present					

			df_lines = pd.DataFrame()

			for population, frequency in\
			 zip(populations, frequencies):
				# print(population)				

				# Storage informations
				# in output dataframe

				# print(csv_gene_file, id_ensembl, variant,\
				#  v_csv, population, frequency)


				# df_line = pd.DataFrame({
				# 'gene': [id_ensembl], \
				# 'variant': [variant], \
				# population: [frequency]})

				# print('df_line: ', df_line)

				# df_output.append(df_line, \
				# 	ignore_index=True)

				# print('df_output: ', df_output)

				li_gen_var_pop_alel_fre: List[str] = split_alelles_from_frequencies(
					id_ensembl,
					variant,
					population,
					frequency,
					debug=False)
				
				with open('output_3.csv', 'a+') as f:					
					# line = '{}>{}>{}>{}\n'.format(
					# 	id_ensembl, \
					# 	variant, \
					# 	population, \
					# 	frequency)

					for line_gen_var_pop_alel_fre in li_gen_var_pop_alel_fre:
						f.write(line_gen_var_pop_alel_fre + '\n')

		# break


	# break

# df_output.to_csv('output.csv', sep='>', index=False)	