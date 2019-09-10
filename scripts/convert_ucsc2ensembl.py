import yaml

with open("config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)

species = data["species"][0]
version = data["genome"][0]

ucsc=open("data/ensembl/" + species + ".vcf")

ucsc2ensembl={}
for line in open("ChromosomeMappings/" + version + "_UCSC2ensembl.txt"):
    linesplit=line.strip().split("\t")
    if len(linesplit) <= 1: continue
    ucsc2ensembl[linesplit[0]] = linesplit[1]

with open("data/ensembl/" + species + ".ensembl.vcf","w") as ensembl:
    for line in ucsc:
        if line.startswith("#"):
            ensembl.write(line)
        splitline = line.split("\t")
        if len(splitline) > 1 and splitline[0] in ucsc2ensembl:
            splitline[0] = ucsc2ensembl[splitline[0]]
            ensembl.write("\t".join(splitline))
