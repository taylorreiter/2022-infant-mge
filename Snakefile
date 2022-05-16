import pandas as pd

m = pd.read_csv("inputs/filereport_read_run_SRR3726337.tsv", header = 0, sep = "\t")
SAMPLES = m['run_accession'].unique().tolist()

rule all:
    input:  expand("outputs/sgc_abund/{sample}.reads.gz.dom_abund.csv", sample = SAMPLES)

###################################################################
## Download reads and databases for workflow
###################################################################

rule download_metagenome_reads:
    output: 
        r1="inputs/raw_reads/{sample}_1.fastq.gz",
        r2="inputs/raw_reads/{sample}_2.fastq.gz"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 600
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastqs = row['fastq_ftp'].values[0]
        fastqs = fastqs.split(";")
        fastq_1 = fastqs[0]
        fastq_2 = fastqs[1]
        if not os.path.exists(params.tmp_base + "_1.fastq.gz"):
            shell("wget -O {output.r1} ftp://{fastq_1}")
        if not os.path.exists(params.tmp_base + "_2.fastq.gz"):
            shell("wget -O {output.r2} ftp://{fastq_2}")


rule download_human_db_for_decontam:
    output: "inputs/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 15
    shell:'''
    wget -O {output} https://osf.io/vsr4n/download
    '''

rule download_CARD_db_for_annotation:
    output: "inputs/broardstreet-v3.2.2.tar.bz2"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 15
    shell:'''
    wget -O {output} https://card.mcmaster.ca/download/0/broadstreet-v3.2.2.tar.bz2
    '''

rule download_viral_db_for_annotation:
    output: "inputs/CHVD_virus_sequences_v1.1.tar.gz"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 15
    shell:'''
    wget -O {output} https://zenodo.org/record/4498884/files/CHVD_virus_sequences_v1.1.tar.gz?download=1
    '''

rule decompress_viral_db_for_annotation:
    input: "inputs/CHVD_virus_sequences_v1.1.tar.gz"
    output: "inputs/CHVD_virus_sequences_v1.1.fasta"
    params: outdir = "inputs"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 2
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule decompress_card_db_for_annotation:
    input: "inputs/broardstreet-v3.2.2.tar.bz2"
    output: "inputs/card_db/nucleotide_fasta_protein_homolog_model.fasta"
    params: outdir = "inputs/card_db"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 2
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule combine_query_sequences_into_single_fasta:
    input: 
        "inputs/card_db/nucleotide_fasta_protein_homolog_model.fasta",
        "inputs/CHVD_virus_sequences_v1.1.fasta"
    output: "inputs/sgc_multifasta_query_hgt_db.fasta"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 2
    shell:'''
    cat {input} > {output}
    '''

rule sketch_query_sequences:
    input: "inputs/sgc_multifasta_query_hgt_db.fasta"
    output:"inputs/sgc_multifasta_query_hgt_db.sig"
    conda: "envs/sourmash.yml"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 2
    shell:'''
    sourmash sketch dna -p k=31,scaled=2000 -o {output} {input}
    '''

###################################################################
## Quality control reads
##   + fastp: remove adapaters, very low quality, and short reads
##   + bbduk: remove host (human) sequences
##   + khmer: k-mer abundance trim 
###################################################################

rule fastp:
    input:
        R1="inputs/raw_reads/{sample}_1.fastq.gz",
        R2="inputs/raw_reads/{sample}_2.fastq.gz"
    output:
        R1="outputs/fastp/{sample}_R1.trim.fq.gz",
        R2="outputs/fastp/{sample}_R2.trim.fq.gz",
        json = "outputs/fastp/{sample}.json",
        html = "outputs/fastp/{sample}.html",
    conda: "envs/fastp.yml"
    threads: 1
    resources:
        mem_mb=16000,
        time_min = 440
    shell:'''
    fastp --in1 {input.R1} \
      --in2 {input.R2} \
      --out1 {output.R1} \
      --out2 {output.R2} \
      --detect_adapter_for_pe \
      --qualified_quality_phred 4 \
      --length_required 31 --correction \
      --json {output.json} \
      --html {output.html}
    '''

rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# original database kept here:
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
# I uploaded it to OSF in case the original disappears and to allow auto download
    input: 
        r1="outputs/fastp/{sample}_R1.trim.fq.gz",
        r2="outputs/fastp/{sample}_R2.trim.fq.gz",
        human='inputs/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    output:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz',
        human_r1='outputs/bbduk/{sample}_R1.human.fq.gz',
        human_r2='outputs/bbduk/{sample}_R2.human.fq.gz'
    conda: 'envs/bbmap.yml'
    threads: 1
    resources:
        mem_mb=64000,
        time_min = 440
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'envs/khmer.yml'
    threads: 1
    resources:
        mem_mb=61000,
        time_min = 880
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

########################################################
## Build and query assembly graph with spacegraphcats
########################################################

rule make_sgc_pangenome_multifasta_conf_files:
    input:
        reads =  "outputs/abundtrim/{sample}.abundtrim.fq.gz",
        ref_genes = "inputs/sgc_multifasta_query_hgt_db.fasta",
        ref_sig = "inputs/sgc_multifasta_query_hgt_db.sig"
    output:
        conf = "outputs/sgc_conf/{sample}_r10_multifasta_conf.yml"
    resources:
        mem_mb = 500,
        time_min = 30
    threads: 1
    run:
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.sample}
input_sequences:
- {input.reads}
radius: 10
paired_reads: true
multifasta_reference:
- {input.ref_genes}
multifasta_scaled: 2000
multifasta_query_sig: {input.ref_sig}
""", file=fp)

rule spacegraphcats_build:
    """
    Build the catlas, the internal sgc graph structure that enables efficient queries
    """
    input:
        reads="outputs/abundtrim/{sample}.abundtrim.fq.gz",
        conf="outputs/sgc_conf/{sample}_r10_multifasta_conf.yml"
    output: 
        "outputs/sgc/{sample}_k31_r10/catlas.csv",
        "outputs/sgc/{sample}_k31/cdbg.gxt",
        "outputs/sgc/{sample}_k31/bcalm.unitigs.db"
    resources: 
        mem_mb = 32000,
        time_min = 1440
    benchmark: "benchmarks/spacegraphcats_build_k31_r10_{sample}.tsv"
    params: outdir = "outputs/sgc"
    conda: "envs/spacegraphcats.yml"
    shell:"""
    python -m spacegraphcats run {input.conf} build --nolock  --outdir={params.outdir} --rerun-incomplete
    """

rule spacegraphcats_multifasta_query:
    input:
        reads =  "outputs/abundtrim/{sample}.abundtrim.fq.gz",
        ref_genes = "inputs/sgc_multifasta_query_hgt_db.fasta", 
        ref_sig = "inputs/sgc_multifasta_query_hgt_db.sig",
        conf = "outputs/sgc_conf/{sample}_r10_multifasta_conf.yml",
        catlas = "outputs/sgc/{sample}_k31_r10/catlas.csv"
    output: 
        "outputs/sgc/{sample}_k31_r10_multifasta/multifasta.cdbg_annot.csv",
        "outputs/sgc/{sample}_k31_r10_multifasta/multifasta.cdbg_by_record.csv"
    params: 
        outdir = "outputs/sgc",
    conda: "envs/spacegraphcats.yml"
    benchmark: "benchmarks/spacegraphcats_multifasta_query_k31_r10_{sample}.tsv"
    resources:
        mem_mb = 32000,
        time_min = 920
    threads: 1
    shell:'''
    python -m spacegraphcats {input.conf} multifasta_query --nolock --outdir {params.outdir} --rerun-incomplete
    '''

rule spacegraphcats_estimate_abundances_of_dominating_sets:
    input:
        cdbg = "outputs/sgc/{sample}_k31/cdbg.gxt",
        catlas = "outputs/sgc/{sample}_k31_r10/catlas.csv",
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/sgc_abund/{sample}.reads.gz.dom_abund.csv"
    conda: "envs/spacegraphcats_dom.yml"
    benchmark: "benchmarks/spacegraphcats_count_dominator_abund_k31_r10_{sample}.tsv"
    resources: 
        mem_mb = 10000,
        time_min = 920
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/sgc/" + wildcards.sample + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/sgc/" + wildcards.sample + "_k31_r10", 
        out_dir = "outputs/sgc_abund/"
    shell:'''
    # TODO: CHECK WHERE THIS SCRIPT LIVES WHEN SGC IS INSTALLED
    /home/tereiter/github/spacegraphcats/scripts/count-dominator-abundance.py {params.cdbg_dir} {params.catlas_dir} --outdir {params.out_dir} {input.reads}
    '''

