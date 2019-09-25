# script to get download ensembl references
# assumes all paths exist (bash script in gui)

import yaml
import subprocess
import sys

with open("config.yaml", 'r') as stream:
    data = yaml.safe_load(stream)

species = data["species"][0].lower()
release = data["release"][0]

path = sys.argv[1] # Mus_musculus.GRCm38

primary = "ftp://ftp.ensembl.org/pub/release-" + release + "//fasta/" + species + "/dna/" + path + ".dna.primary_assembly.fa.gz"
toplevel = "ftp://ftp.ensembl.org/pub/release-" + release + "//fasta/" + species + "/dna/" + path + ".dna.toplevel.fa.gz"

# download gfa
if subprocess.check_output(['./validate.sh', primary]) == b'true\n':
    subprocess.check_call(["wget", "-P", "data/ensembl/", "-P", "data/ensembl/", primary])
else:
    subprocess.check_call(["wget", "-P", "data/ensembl/", toplevel])

# download gff
gff = "ftp://ftp.ensembl.org/pub/release-" + release + "/gff3/" + species + "/" + path + "." + release + ".gff3.gz"
subprocess.check_call(["wget", "-P", "data/ensembl/", gff])

# download pep
pep = "ftp://ftp.ensembl.org/pub/release-" + release + "//fasta/" + species + "/pep/" + path + ".pep.all.fa.gz"
subprocess.check_call(["wget", "-P", "data/ensembl/", pep])

# download vcf
vcf = "ftp://ftp.ensembl.org/pub/release-" + release + "/variation/vcf/" + species + "/" + species + ".vcf.gz" # edit, incosistent naming convention /vcf/Mus_musculus.vcf.gz or /vcf/mus_musculus.vcf.gz
subprocess.check_call(["wget", "-P", "data/ensembl/", vcf])
