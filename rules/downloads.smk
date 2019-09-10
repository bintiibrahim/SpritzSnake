REF=config["species"][0] + "." + config["genome"][0]

rule download_ensembl_references:
    output:
        gfa="data/ensembl/" + REF + ".dna.primary_assembly.fa",
        gff="data/ensembl/" + REF + "." + config["ensembl"][0] + ".gff3",
        pfa="data/ensembl/" + REF + ".pep.all.fa",
        vcf="data/ensembl/" + config["species"][0] + ".vcf",
        vcfidx="data/ensembl/" + config["species"][0] + ".vcf.idx"
    log: "data/ensembl/downloads.log"
    shell:
        """
            (wget -O - ftp://ftp.ensembl.org/pub/release-97//fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz |
            gunzip -c > {output.gfa} &&
            wget -O - ftp://ftp.ensembl.org/pub/release-97/gff3/mus_musculus/Mus_musculus.GRCm38.97.gff3.gz |
            gunzip -c > {output.gff} &&
            wget -O - ftp://ftp.ensembl.org/pub/release-97//fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz |
            gunzip -c > {output.pfa} &&
            wget -O - ftp://ftp.ensembl.org/pub/release-97/variation/vcf/mus_musculus/mus_musculus.vcf.gz |
            gunzip -c > {output.vcf} &&
            gatk IndexFeatureFile -F {output.vcf}) 2> {log}
        """

rule download_chromosome_mappings:
    output: "ChromosomeMappings/" + config["genome"][0] + "_UCSC2ensembl.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.fa"
    output: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule convert_ucsc2ensembl:
    input:
        "data/ensembl/" + config["species"][0] + ".vcf",
        "ChromosomeMappings/" + config["genome"][0] + "_UCSC2ensembl.txt"
    output:
        "data/ensembl/" + config["species"][0] + ".ensembl.vcf",
    script:
        "../scripts/convert_ucsc2ensembl.py"

rule index_ucsc2ensembl:
    input: "data/ensembl/" + config["species"][0] + ".ensembl.vcf"
    output: "data/ensembl/" + config["species"][0] + ".ensembl.vcf.idx"
    shell: "gatk IndexFeatureFile -F {input}"

rule filter_gff3:
    input: "data/ensembl/" + REF + "." + config["ensembl"][0] + ".gff3"
    output: "data/ensembl/202122.gff3"
    shell: "grep \"^#\|20\|^21\|^22\" \"data/ensembl/" + REF + "." + config["ensembl"][0] + ".gff3\" > \"data/ensembl/202122.gff3\""

rule filter_fa:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.fa"
    output: "data/ensembl/202122.fa"
    script: "../scripts/filter_fasta.py"

rule download_sras:
    output:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"), # constrain wildcards, so it doesn't soak up SRR######.trim_1.fastq
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
    log: "{dir}/{sra}.log"
    threads: 4
    shell:
        "fasterq-dump --progress --threads {threads} --split-files --outdir {wildcards.dir} {wildcards.sra} 2> {log}"

rule compress_fastqs:
     input:
         temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"),
         temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
     output:
         "{dir}/{sra,[A-Z0-9]+}_1.fastq.gz",
         "{dir}/{sra,[A-Z0-9]+}_2.fastq.gz"
     shell:
         "gzip {input}"
