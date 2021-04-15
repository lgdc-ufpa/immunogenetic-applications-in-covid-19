import re
from typing import List, Tuple


def allele_regex(alleles: str):
	"""Alleles regex

	Arguments
	---------

	alleles: string wich contains alleles, frequency
		and counting. Egg: C: 0.9999960193 (251211) G: 3.9807e-06 (1)

	-> Return: regex.findall(alleles)
	"""
	# 
	regex = re.compile(r"(\w+)[: ] ([A-Za-z0-9_]+\.[A-Za-z0-9_]+) \(([A-Za-z0-9_]+)\)")	
	# return regex.search(alleles).groups()
	return regex.findall(alleles)


# alleles = "C: 0.9999960193 (251211) G: 3.9807e-06 (1)"
# match = allele_regex(alleles)
# print(match)


def alleles_split(alleles: str) -> List[Tuple]:
	"""Split alleles in [(allele, frequency, counting)]
	"""
	t1 = map(lambda x: x.replace('(', ''), alleles.split(')'))
	t2 = map(lambda x: x.replace(':', ''), t1)
	t3 = map(lambda x: x.strip(), t2)
	t4 = filter(lambda x: x != '', t3)	
	t5: List[Tuple] = list(map(lambda x: tuple(x.split(' ')), t4))	
	return t5
		
	# for elem in t4:
	# 	allele, freq, count = elem.split(' ')
	# 	print(allele, freq, count)

# alleles = "C: 0.9999960193 (251211) G: 3.9807e-06 (1)"
# li_alleles = alleles_split(alleles)
# print(li_alleles)