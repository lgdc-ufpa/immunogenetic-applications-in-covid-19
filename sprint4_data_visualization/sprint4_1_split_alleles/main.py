li_alelles = ["G: 1.000 (4378) A: 0.000", \
"CCC: 0.9999534898 (64499)   C: 1.55034e-05 (1)   CC: 3.10068e-05 (2)"]

def calc_alelle_frequency() -> dict:
	"""Split Allele Frequency per variant
	:return -> dict
	"""
	di: dict = {}	

	for elem in li_alelles:
		alelle_name: str = None

		for item in elem.split(sep=' '):

			if ":" in item:
				alelle_name: str = item.replace(":", "")

				di[alelle_name] = []

			elif not item == "":
				di[alelle_name].append(item)

	return di

di_alelle_frequency = calc_alelle_frequency()
print(di_alelle_frequency)
