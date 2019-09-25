REF=config["species"][0] + "." + config["genome"][0]

rule download_ensembl_references:
    input: REF
    output:
        gfagz=temp("data/ensembl/" + REF + ".dna.primary_assembly.fa.gz"),
        gffgz=temp("data/ensembl/" + REF + "." + config["ensembl"][0] + ".gff3.gz"),
        pfagz=temp("data/ensembl/" + REF + ".pep.all.fa.gz"),
        vcfgz=temp("data/ensembl/" + config["species"][0] + ".vcf.gz"),
        gfa="data/ensembl/" + REF + ".dna.primary_assembly.fa",
        gff="data/ensembl/" + REF + "." + config["ensembl"][0] + ".gff3",
        pfa="data/ensembl/" + REF + ".pep.all.fa",
        vcf="data/ensembl/" + config["species"][0] + ".vcf",
        vcfidx="data/ensembl/" + config["species"][0] + ".vcf.idx"
    log: "data/ensembl/downloads.log"
    shell:
        """
            python scripts/download_ensembl.py {input} &&
            gunzip {output.gfagz} > {output.gfa} &&
            gunzip {output.gffgz} > {output.gff} &&
            gunzip {output.pfagz} > {output.pfa} &&
            gunzip {output.vcfgz} > {output.vcf} &&
            gatk IndexFeatureFile -F {output.vcf}) 2> {log}
        """

rule download_chromosome_mappings:
    output: "ChromosomeMappings/" + config["genome"][0] + "_UCSC2ensembl.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.fa"
    output: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule dict_fa:
    input: "data/ensembl/" + config["species"][0] + "." + config["genome"][0] + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/" + config["species"][0] + "." + config["genome"][0] + ".dna.primary_assembly.karyotypic.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule tmpdir:
    output: temp(directory("tmp"))
    shell: "mkdir tmp"

rule convert_ucsc2ensembl:
    input:
        "data/ensembl/" + config["species"][0] + ".vcf",
        "ChromosomeMappings/" + config["genome"][0] + "_UCSC2ensembl.txt",
        tmp=directory("tmp"),
        fa="data/ensembl/" + config["species"][0] + "." + config["genome"][0] + ".dna.primary_assembly.karyotypic.fa",
        dict="data/ensembl/" + config["species"][0] + "." + config["genome"][0] + ".dna.primary_assembly.karyotypic.dict",
    output:
        ensVcf=temp("data/ensembl/" + config["species"][0] + ".orig.ensembl.vcf"),
        dictVcf="data/ensembl/" + config["species"][0] + ".ensembl.vcf",
    shell:
        "python scripts/convert_ucsc2ensembl.py && "
        "gatk UpdateVCFSequenceDictionary -R {input.fa} --sequence-dictionary {input.dict} -V {output.ensVcf} --output {output.dictVcf} --tmp-dir {input.tmp}"

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
