# see https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf

import requests
import sys

# find proteome
BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/proteomes/'
TOOL_ENDPOINT = '/uploadlists/'

query = 'Mus_musculus' # read config

payload = {
    'query': query,
    'sort': 'score',
    'format': 'list',
    }

proteome_res = requests.get(BASE + KB_ENDPOINT, params=payload, stream=True)
proteome_res.raise_for_status() # throw an error for bad status code

if proteome_res.text.count('\n') > 0:
    proteome = proteome_res.text.split('\n')[0] # get first proteome result
# else throw error
