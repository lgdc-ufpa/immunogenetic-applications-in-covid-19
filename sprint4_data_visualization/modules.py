from typing import List

def split_alelles_from_frequencies(\
	gene, variant, population,\
	alleles_and_frequencies: str,\
	debug: bool = True):
	"""Split Allele Frequency per variant
	:df_alelle_freq: pandas dataframe (or serie)
	:return -> line to save in the csv output

	line == gene | variant | population | alelle | frequency
	"""

	# Cleaning trash data
	population = population.replace(" ", "")

	li_gen_var_pop_alel_fre: List[str] = []
	di: dict = {}

	for item in alleles_and_frequencies.split(sep=' '):

		if ":" in item:
			alelle_name: str = item.replace(":", "")

			di[alelle_name] = []

		elif not item == "":
			# if debug:
			# 	print(item)

			di[alelle_name].append(item)

	# return di

	for key_alelle, value_frequency in di.items():
		# print(gene, variant, population, key_alelle, value_frequency)
		line:str = None

		if len(value_frequency) == 1:

			line = f'{gene}>{variant}>{population}>{key_alelle}>{value_frequency[0]}>0'
			if debug:
				print(line)

		else:
			count_alelle = value_frequency[1].replace('(', '').replace(')', '')
			line = f'{gene}>{variant}>{population}>{key_alelle}>{value_frequency[0]}>{count_alelle}'
			if debug:
				print(line)

		li_gen_var_pop_alel_fre.append(line)

	if debug:
		print("--------")

	return li_gen_var_pop_alel_fre