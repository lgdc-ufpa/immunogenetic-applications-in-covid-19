from bs4 import BeautifulSoup
from urllib.parse import urljoin
import requests
from lxml import etree
from lxml import html
from tqdm import tqdm


with open('gene_link-orphanet_link-ensembl_id-ensembl_chromossome.csv', "w") as f:
    f.write('Gene>Link orphanet>Link Ensembl>Código no Ensembl>Localização Cromossômica>\n')
    f.close()

# soup = BeautifulSoup(open("links_genes.html", "r"), 'html.parser')
soup = BeautifulSoup(open("links_genes.html", "r"), 'html.parser')

list_li = soup.find_all('li')

for li in tqdm(list_li):
	gene_name = li.text
	# print(gene_name)
	# link_orphanet = li.find('a', href=True).get('href')
	base = "https://www.orpha.net/consor/cgi-bin/"
	link_orphanet =  urljoin(base, li.find('a')['href'].replace(' ', '%20'))
	# print(link_orphanet)

	req_gene = requests.get(link_orphanet)
	# tree = etree.parse(req_gene, etree.HTMLParser())
	# tree = html.fromstring(req_gene.content)
	# trs = tree.xpath('//*[@id="ContentType"]/div[2]/ul/li[9]')

	soup2 = BeautifulSoup(req_gene.text, 'html.parser')
	chromossome_region = soup2.find_all('strong')[4].text
	ommin_link = soup2.find_all('strong')[5]

	id_ensembl = soup2.find_all('strong')[9].text
	link_ensembl = soup2.find_all('strong')[9].find('a')['href'].replace(' ', '%20')
	
	with open('gene_link-orphanet_link-ensembl_id-ensembl_chromossome.csv', "a") as f:
		f.write(f'{gene_name}>{link_orphanet}>{link_ensembl}>{id_ensembl}>{chromossome_region}\n')
		f.close()