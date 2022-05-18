import pandas as pd
import os

m = pd.read_csv("inputs/filereport_read_run_SRR3726337.tsv", header = 0, sep = "\t")
SAMPLES = m['run_accession'].unique().tolist()

# TODO:
# temporarily define genome accession wildcard for BLAST queries
# This is straightforward to automate gather matches -> genome download, but this is a faster bandaid
# Requires new class to be defined that will scrape in df of gather results and solve appropriate accessions per sample
GENOMES = ['GCA_000701705.1_ASM70170v1', 'GCA_001913975.1_ASM191397v1'  'GCA_009737075.1_ASM973707v1',
           'GCA_001913895.1_ASM191389v1', 'GCA_001931885.1_ASM193188v1', 'GCA_019173815.2_ASM1917381v2'] 

rule all:
    input:  
        expand("outputs/sourmash_gather/{sample}_k31_scaled2000_gtdb-rs207.csv", sample = SAMPLES),
        expand("outputs/sgc_abund/{sample}.abundtrim.fq.gz.dom_abund.csv", sample = SAMPLES), 
        expand("outputs/sgc/{sample}_k31_r10_multifasta/multifasta.cdbg_annot.csv", sample = SAMPLES),
        expand("outputs/bandage/{sample}_done_mge.txt", sample = SAMPLES),
        expand("outputs/bandage/{sample}_{genome}_done_species.txt", sample = SAMPLES, genome = GENOMES)

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
        if not os.path.exists(output.r1):
            shell("wget -O {output.r1} ftp://{fastq_1}")
        if not os.path.exists(output.r2):
            shell("wget -O {output.r2} ftp://{fastq_2}")


rule download_sourmash_db_for_taxonomy:
    """
    NOTE: ran sourmash gather using cluster-local copy of db, added rule for provenance
    """
    output: "inputs/gtdb-rs207.genomic.k31.zip"
    threads: 1
    resources:
        mem_mb = 800,
        time_min = 30
    shell:'''
    wget -O {output} whttps://osf.io/k2u8s/download
    '''

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
    sourmash sketch dna -p k=31,scaled=1000 -o {output} {input}
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
## Determine known taxonomic composition of sample
########################################################

rule sourmash_sketch_mgx:
    input: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/sourmash_sigs/{sample}.sig"
    conda: "envs/sourmash.yml"
    threads: 1
    resources: 
        mem_mb = 800,
        time_min = 30
    shell:'''
    sourmash sketch dna -p k=31,scaled=2000 -o {output} {input}
    '''
    
rule sourmash_gather_mgx:
    input:
        sig= "outputs/sourmash_sigs/{sample}.sig",
        db = "/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.genomic.k31.zip" 
    output: "outputs/sourmash_gather/{sample}_k31_scaled2000_gtdb-rs207.csv"
    conda: "envs/sourmash.yml"
    benchmark: "benchmarks/sourmash_gather_k31_scaled2000_gtdb-rs207_{sample}.tsv"
    threads: 1
    resources: 
        mem_mb = 32000,
        time_min = 240
    shell:'''
    sourmash gather -o {output} -k 31 {input.sig} {input.db}
    '''

########################################################
## Build and query assembly graph with spacegraphcats
########################################################

rule make_sgc_multifasta_conf_files:
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
multifasta_scaled: 1000
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
    python -m spacegraphcats run {input.conf} multifasta_query --nolock --outdir {params.outdir} --rerun-incomplete
    '''

rule spacegraphcats_estimate_abundances_of_dominating_sets:
    input:
        cdbg = "outputs/sgc/{sample}_k31/cdbg.gxt",
        catlas = "outputs/sgc/{sample}_k31_r10/catlas.csv",
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/sgc_abund/{sample}.abundtrim.fq.gz.dom_abund.csv"
    conda: "envs/spacegraphcats.yml"
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

#################################################################
## extract and visualize candidate MGEs
#################################################################

rule grab_names_of_candidate_prevalent_mges:
    input: cdbg_by_record = "outputs/sgc/{sample}_k31_r10_multifasta/multifasta.cdbg_by_record.csv"
    output: candidate_mges =  "outputs/candidate_mges/{sample}_candidate_mge_names.txt"
    conda: "envs/tidyverse.yml"
    resources: 
        mem_mb = 4000,
        time_min = 20
    threads: 1
    script: "scripts/snakemake_grab_names_of_candidate_prevalent_mges.R"

rule extract_candidate_mges_from_fasta:
    input: 
        candidate_mges =  "outputs/candidate_mges/{sample}_candidate_mge_names.txt",
        fasta = "inputs/sgc_multifasta_query_hgt_db.fasta"
    output: "outputs/candidate_mges/{sample}_candidate_mges.fasta"
    conda: "envs/seqtk.yml"
    resources: 
        mem_mb = 4000,
        time_min = 20
    threads: 1
    shell:'''
    seqtk subseq {input.fasta} {input.candidate_mges} > {output}
    '''

checkpoint separate_candidate_mges_into_single_fastas:
    """
    Note will expand over samples here, making queries for all discovered MGEs across all samples if there are multiple
    """
    input: expand("outputs/candidate_mges/{sample}_candidate_mges.fasta", sample = SAMPLES)
    output: directory("outputs/candidate_mge_sequences")    
    resources: 
        mem_mb = 4000,
        time_min = 20
    threads: 1
    shell:"""
    # This should probably directly encode write folder in awk instead of having an mv command at the end
    # also not the best form to rename things with sed...probably should replace this whole rule with python something or another
    mkdir -p outputs/candidate_mge_sequences
    cat {input} | awk '{{
        if (substr($0, 1, 1)==">") {{filename=(substr($0,2) ".fa")}}
        print $0 > filename
    }}'
    for file in *fa; do mv "$file" $(echo "$file" | sed -e 's/[^A-Za-z0-9._-]/_/g'); done
    mv *.fa outputs/candidate_mge_sequences
    """

def checkpoint_separate_candidate_mges_into_single_fasta_1(wildcards):
    # Expand checkpoint to get fasta names, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.separate_candidate_mges_into_single_fastas.get(**wildcards).output[0]
    file_names = expand("outputs/candidate_mge_sequences/{mge}.fa",
                        mge = glob_wildcards(os.path.join(checkpoint_output, "{mge}.fa")).mge)
    return file_names

rule make_sgc_extract_conf_files:
    input:
        reads =  "outputs/abundtrim/{sample}.abundtrim.fq.gz",
        queries = checkpoint_separate_candidate_mges_into_single_fasta_1
    output:
        conf = "outputs/sgc_conf/{sample}_r10_extract_conf.yml"
    resources:
        mem_mb = 500,
        time_min = 30
    threads: 1
    run:
        query_list = "\n- ".join(input.queries)
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.sample}
input_sequences:
- {input.reads}
ksize: 31
radius: 10
paired_reads: true
search:
- {query_list}
""", file=fp)


checkpoint spacegraphcats_extract_mge_sequences:
    input: 
        conf = "outputs/sgc_conf/{sample}_r10_extract_conf.yml",
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: directory("outputs/sgc/{sample}_k31_r10_search_oh0/")
        #"outputs/sgc/{sample}_k31_r10_search_oh0/results.csv",
        #"outputs/sgc/{sample}_k31_r10_search_oh0/{mge}.fa.cdbg_ids.reads.gz",
        #"outputs/sgc/{sample}_k31_r10_search_oh0/{mge}.fa.contigs.sig"
    benchmark: "benchmarks/spacegraphcats_extract_contigs_{sample}.tsv"
    params: outdir = "outputs/sgc"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 32000,
        time_min = 960
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

rule bcalm_mge_nbhd:
    input: "outputs/sgc/{sample}_k31_r10_search_oh0/{mge}.fa.cdbg_ids.reads.gz"
    output: "outputs/bcalm/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.fa"
    params: prefix = "outputs/bcalm/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz"
    resources:
        mem_mb = 4000,
        time_min = 120
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    bcalm -in {input} -out-dir outputs/bcalm/{wildcards.sample}_r10 -kmer-size 31 -abundance-min 1 -out {params.prefix}
    '''

rule convert_to_gfa:
    input: "outputs/bcalm/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.fa"
    output: "outputs/bcalm/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.gfa"
    resources:
        mem_mb = 2000,
        time_min = 60
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    python ./scripts/convertToGFA.py {input} {output} 31
    '''

rule bandage_plot_gfa_mge:
    input: 
        gfa="outputs/bcalm/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.gfa",
        blastquery="outputs/candidate_mge_sequences/{mge}.fa"
    output: "outputs/bandage/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.png"
    resources:
        mem_mb = 4000,
        time_min = 120
    threads: 1
    conda: "envs/bandage.yml"
    shell:'''
    XDG_RUNTIME_DIR=~/tmp # REQUIRES THAT THIS FOLDER EXISTS
    Bandage image {input.gfa} {output} --query {input.blastquery}    
    '''

    
def checkpoint_spacegraphcats_extract_mge_sequences_1(wildcards):
    # Expand checkpoint to get fasta names, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.spacegraphcats_extract_mge_sequences.get(**wildcards).output[0]    
    file_names = expand("outputs/bandage/{{sample}}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.png",
                        mge = glob_wildcards(os.path.join(checkpoint_output, "{mge}.fa.cdbg_ids.reads.gz")).mge)
    return file_names

rule dummy_solve_mge_blast:
    input: checkpoint_spacegraphcats_extract_mge_sequences_1
    output: touch("outputs/bandage/{sample}_done_mge.txt")
    
rule bandage_plot_gfa_species:
    input: 
        gfa="outputs/bcalm/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs.gfa",
        blastquery="outputs/genbank_genomes/{genome}_genomic.fna"
    output: "outputs/bandage/{sample}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs_{genome}.png"
    resources:
        mem_mb = 4000,
        time_min = 120
    threads: 1
    conda: "envs/bandage.yml"
    shell:'''
    XDG_RUNTIME_DIR=~/tmp # REQUIRES THAT THIS FOLDER EXISTS
    Bandage image {input.gfa} {output} --query {input.blastquery}    
    '''

def checkpoint_spacegraphcats_extract_mge_sequences_2(wildcards):
    # Expand checkpoint to get fasta names, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.spacegraphcats_extract_mge_sequences.get(**wildcards).output[0]    
    file_names = expand("outputs/bandage/{{sample}}_r10/{mge}.fa.cdbg_ids.reads.gz.unitigs_{{genome}}.png",
                        mge = glob_wildcards(os.path.join(checkpoint_output, "{mge}.fa.cdbg_ids.reads.gz")).mge)
    return file_names

rule dummy_solve_species_blast:
    input: checkpoint_spacegraphcats_extract_mge_sequences_2
    output: touch("outputs/bandage/{sample}_{genome}_done_species.txt")
