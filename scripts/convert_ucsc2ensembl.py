ucsc=open("data/ensembl/Mus_musculus.vcf")

ucsc2ensembl={}
for line in open("ChromosomeMappings/GRCm38_UCSC2ensembl.txt"):
    linesplit=line.strip().split("\t")
    if len(linesplit) <= 1: continue
    ucsc2ensembl[linesplit[0]] = linesplit[1]

with open("data/ensembl/Mus_musculus.ensembl.vcf","w") as ensembl:
    for line in ucsc:
        if line.startswith("#"):
            ensembl.write(line)
        splitline = line.split("\t")
        if len(splitline) > 1 and splitline[0] in ucsc2ensembl:
            splitline[0] = ucsc2ensembl[splitline[0]]
            ensembl.write("\t".join(splitline))
