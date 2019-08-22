# rule download_ensembl_references:
#     output:
#         gfa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#         gff="data/ensembl/Homo_sapiens.GRCh38.81.gff3",
#         pfa="data/ensembl/Homo_sapiens.GRCh38.pep.all.fa",
#         vcf="data/ensembl/" + SPECIES + ".vcf",
#         vcfidx="data/ensembl/common_all_20170710.vcf.idx"
#     log: "data/ensembl/downloads.log"
#     shell:
#         "(wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | "
#         "gunzip -c > {output.gfa} && "
#         "wget -O - ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz | "
#         "gunzip -c > {output.gff} && "
#         "wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz | "
#         "gunzip -c > {output.pfa} && "
#         "wget -O - ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz | "
#         "gunzip -c > {output.vcf} && "
#         "gatk IndexFeatureFile -F {output.vcf}) 2> {log}"
#
# def check_ensembel_references(wildcards):
#     # if files exist in folder

rule download_ensembl_references:
    output:
        gfa="data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.fa",
        gff="data/ensembl/" + SPECIES + "." + GENEMODEL_VERSION + ".gff3",
        pfa="data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".pep.all.fa",
        vcf="data/ensembl/" + SPECIES + ".vcf",
        vcfidx="data/ensembl/" + SPECIES + ".vcf.idx"
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
    output: "ChromosomeMappings/" + GENOME_VERSION + "_UCSC2ensembl.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.fa"
    output: "data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule dict_fa:
    input: "data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule tmpdir:
    output: temp(directory("tmp"))
    shell: "mkdir tmp"

rule convert_ucsc2ensembl:
    input:
        "data/ensembl/" + SPECIES + ".vcf",
        "ChromosomeMappings/" + GENOME_VERSION + "_UCSC2ensembl.txt",
        tmp=directory("tmp"),
        fa="data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa",
        dict="data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.dict",
    output:
        ensVcf=temp("data/ensembl/" + SPECIES + ".orig.ensembl.vcf"),
        dictVcf="data/ensembl/" + SPECIES + ".ensembl.vcf",
    shell:
        "python scripts/convert_ucsc2ensembl.py && "
        "gatk UpdateVCFSequenceDictionary -R {input.fa} --sequence-dictionary {input.dict} -V {output.ensVcf} --output {output.dictVcf} --tmp-dir {input.tmp}"

rule index_ucsc2ensembl:
    input: "data/ensembl/" + SPECIES + ".ensembl.vcf"
    output: "data/ensembl/" + SPECIES + ".ensembl.vcf.idx"
    shell: "gatk IndexFeatureFile -F {input}"

rule filter_gff3:
    input: "data/ensembl/" + SPECIES + "." + GENEMODEL_VERSION + ".gff3"
    output: "data/ensembl/202122.gff3"
    shell: "grep \"^#\|20\|^21\|^22\" \"data/ensembl/" + SPECIES + "." + GENEMODEL_VERSION + ".gff3\" > \"data/ensembl/202122.gff3\""

rule filter_fa:
    input: "data/ensembl/" + SPECIES + "." + GENOME_VERSION + ".dna.primary_assembly.fa"
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
