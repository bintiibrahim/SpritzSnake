# script to get download ensembl references
# assumes all paths exist (bash script in gui)

import yaml
import subprocess
import sys

with open("config.yaml", 'r') as stream:
    data = yaml.safe_load(stream)

species = data["species"][0]
release = data["release"][0]
path = sys.argv[0] # Mus_musculus.GRCm38

gfa1 = "ftp://ftp.ensembl.org/pub/release-" + release + "//fasta/" + species + "/dna/" + path + ".dna.primary_assembly.fa.gz"
gfa2 = "ftp://ftp.ensembl.org/pub/release-" + release + "//fasta/" + species + "/dna/" + path + ".dna.toplevel.fa.gz"

# if gfa1 path exists, download
# else gfa2

gff = "ftp://ftp.ensembl.org/pub/release-" + release + "/gff3/" + species + "/" + path + "." + release + ".gff3.gz"
# download gff

pep = "ftp://ftp.ensembl.org/pub/release-" + release + "//fasta/" + species + "/pep/" + path + ".pep.all.fa.gz"
# download pep

vcf = "ftp://ftp.ensembl.org/pub/release-" + release + "/variation/vcf/" + species + "/" + species + ".vcf.gz"
# download vcf

# gunzip
