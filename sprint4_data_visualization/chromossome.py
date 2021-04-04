import re
from typing import Tuple, List, Iterable
import pandas as pd


def chromossome_regex(chromossome_str: str) -> Tuple[str, List[str], List[str]]:
	"""chromossome regex function

	Arguments
	---------

	chromossome_str : chromossome

	-> return : (chromossome, [arm, region, band, subband], [arm2, region2, band2, subband2])
	
	Examples
	--------

	chromossome1 = chromossome_regex('19p13.2')
	chromossome2 = chromossome_regex('19p13.16')
	chromossome3 = chromossome_regex('6q21')
	chromossome4 = chromossome_regex('Xq25')
	chromossome5 = chromossome_regex('Xp22.11-p21.3')
	chromossome6 = chromossome_regex('1q31.3-q32.1')

	print(chromossome1)
	print(chromossome2)
	print(chromossome3)
	print(chromossome4)
	print(chromossome5)
	print(chromossome6)

	Output
	------

	('19', ['p', '1', '3', '2'], [None, None, None, None])
	('19', ['p', '1', '3', '16'], [None, None, None, None])
	('6', ['q', '2', '1', None], [None, None, None, None])
	(None, ['q', '2', '5', None], [None, None, None, None])
	(None, ['p', '2', '2', '11'], ['p', '2', '1', '3'])
	('1', ['q', '3', '1', '3'], ['q', '3', '2', '1'])

	"""
	regex = re.compile(r'(\d+)?(\w)(\d)(\d)\.?(\d+)?-?(\w)?(\d)?(\d)?\.?(\d+)?')	
	matching_object = regex.search(chromossome_str)
	chromossome, arm_length1, region1, band1, subband1, arm_length2, region2, band2, subband2 = matching_object.groups()
	arm_1 = [arm_length1, region1, band1, subband1]
	arm_2 = [arm_length2, region2, band2, subband2]
	return chromossome, arm_1, arm_2

# chromossome1 = chromossome_regex('19p13.2')
# chromossome2 = chromossome_regex('19p13.16')
# chromossome3 = chromossome_regex('6q21')
# chromossome4 = chromossome_regex('Xq25')
# chromossome5 = chromossome_regex('Xp22.11-p21.3')
# chromossome6 = chromossome_regex('1q31.3-q32.1')

# print(chromossome1)
# print(chromossome2)
# print(chromossome3)
# print(chromossome4)
# print(chromossome5)
# print(chromossome6)


def chromossome_sequence_to_dataframe(chromossomes: Iterable) -> pd.DataFrame:
	"""Generate chromossome dataframe
	
	Arguments
	---------

	chromossome_sequence : chromossome sequence

	-> return : dataframe of chromossomes

	"""
	di_chromossome = {}
	di_chromossome['chromossome'] = []
	di_chromossome['arm_1'] = []
	di_chromossome['region_1'] = []
	di_chromossome['band_1'] = []
	di_chromossome['subband_1'] = []
	di_chromossome['arm_2'] = []
	di_chromossome['region_2'] = []
	di_chromossome['band_2'] = []
	di_chromossome['subband_2'] = []

	for chrmssm in chromossomes:	
		chromossome, arm_1, arm_2 = chromossome_regex(chrmssm)		
		di_chromossome['chromossome'].append(chromossome)
		di_chromossome['arm_1'].append(arm_1[0])
		di_chromossome['region_1'].append(arm_1[1])
		di_chromossome['band_1'].append(arm_1[2])
		di_chromossome['subband_1'].append(arm_1[3])
		di_chromossome['arm_2'].append(arm_2[0])
		di_chromossome['region_2'].append(arm_2[1])
		di_chromossome['band_2'].append(arm_2[2])
		di_chromossome['subband_2'].append(arm_2[3])

	return pd.DataFrame(di_chromossome)