# Workflow to identity prevalent mobile genetic elements in a metagenome

The goal of this repository is to demonstrate a workflow for the discovery of a prevalent mobile genetic elements (MGEs) in a metagenome.

The workflow uses a metagenome assembly graph to find MGEs and their sequence context.
For a primer on metagenomics and assembly graphs, see [here](https://spacegraphcats.github.io/spacegraphcats/0a-primer/).

**Abbreviations and definiation in this document**

- cDBG: compact de Bruijn Graph; a type of assembly graph
- MGE: mobile genetic element
- k-mer: words of length *k* in nucleotide sequences.

## Getting started with this repository

This repository uses snakemake to execute the complete workflow. 

Snakemake and pandas must be installed in the run environment for the workflow to run.
We recommend using conda to create this environment. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

To get started, create and activate a conda environment based on the `environment.yml` file packaged with this repository, using the commands:
```
conda env create --name infant --file environment.yml
conda activate infant
```

Once the environment is installed and activated, the analysis and further software installations are automated with snakemake.
A snakemake workflow can be executed with many parameters, but the command below is sufficient to get it running:

```
snakemake -j 1 --use-conda --rerun-incomplete
```

Snakemake can also parallelize job submission and modulate resource usage (RAM, CPUs, run time). 
We used the command below on a slurm cluster, but other cluster engines are also supported.

```
snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 \
    --cluster "sbatch -t {resources.time_min} -J inf -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

### Expanding to additional samples

This workflow is built around metadata tables (tsv format) exported from the European Nucleotide Archive (ENA). 
To add additional samples, export the ENA metadata table, being sure to include the "run_accession" and "fastq" columns, and replace the metadata read in by the snakefile with the path to your new metadata table.

Currently, the workflow is written for paired end data only.
An example of how to adapt the workflow to accommodate single end or paired end reads can be found [here](https://github.com/greenelab/2022-microberna/blob/main/pangenome_reference.snakefile#L444).

## Overview of Approach

![](https://i.imgur.com/1bUVynT.png)

This workflow uses assembly graph queries to identify and estimate the abundance of mobile genetic elements in a metagenome.
The workflow begins by quality controlling the raw metagenome sequence reads. 
[fastp](https://github.com/OpenGene/fastp) is used for adapter removal and general quality control, [bbduk](http://seqanswers.com/forums/archive/index.php/t-42552.html) is used to remove human (host) sequences, and khmer [trim-low-abund.py](https://khmer-recipes.readthedocs.io/en/latest/007-variable-coverage-trimming/) is used to k-mer trim the reads. 
Then, [spacegraphcats](https://doi.org/10.1186/s13059-020-02066-4) is used to build the assembly graph and to enable efficient querying by determining a dominating set.
A dominating set is a subset of nodes in the assembly graph where every other node is a maximum distance of *r* away from a dominator.
The dominating set is used to carve the assembly graph into pieces. 
Because the pieces are larger than the nodes, this structure is more efficient to search and query.
For this workflow, we use a radius (*r*) of 10 to build the dominating set, generating larger pieces that can give us more of the context surrounding MGEs. 
(For a brief comparison of the performance of different radius sizes when investigating MGEs, see [here](https://spacegraphcats.github.io/spacegraphcats/02-spacegraphcats-use-cases/#identifying-the-context-of-eg-horizontally-transferred-genes).)
Using this graph structure, we then search the graph for any dominating set that overlaps with any k-mers in a query using spacegraphcats `multifasta_queries`.
Currently, the workflow uses all sequences in the [CARD database](https://card.mcmaster.ca/download) and the [CHVD phage database](https://zenodo.org/record/4498884/#.YoPcM2DMLao) as queries.
Any nucleotide database with query sequences of interest could be substituted for these databases.

After identifying antibiotic resistance elements and viruses in the metagenome, the workflow selects abundant sequences to explore further.
It selects up to 30 sequences from each database, selecting the top 10 sequences that had the highest mean abundance (mean abundance of all k-mers in the dominating set piece), total abundance (total abundance of all k-mers in the dominating set piece), and total k-mers (most diverse matches). 
We arbitrarily selected the top 10 sequences to reduce run time; in principle, all sequences could be explored at this phase.
We refer to the selected sequences as candidate MGEs.

To understand whether these MGEs are not only abundant in the metagenome but also occur in complex backgrounds (e.g. are associated with many distinct genomic backgrounds), we next explored the sequence in the neighborhoods around the candidate MGEs. 
To do this, we queried with the candidate MGE sequences and extracted the sequences in the entire dominating set piece for which there were overlaps using spacegraphcats `extract_reads`.
This is conceptually similar to the original `multifasta_query` search that we performed, but this time we actually extracted the reads that contained sequences that were in the matching dominating set pieces.
We then built mini assembly graphs with each match and assessed the complexity of the graph structure.
When a graph has more branches and loops, it indicates that there is more diversity in that sequencing neighborhood.
Higher complexity, especially when paired with high abundance, is indicative that the sequence occurs many times in the metagenome in many different contexts (e.g. genomic backgrounds).

Lastly, we assessed the taxonomic identity of the sequences that surrounded the candidate MGEs.
We determined the the taxonomic composition of the metagenome using sourmash `gather` against the GTDB rs207 database; while this database only includes bacteria and archaea, there are sourmash GenBank databases for viruses, protozoa, and fungi as well that could easily be added to this step.
(see [issue #1](https://github.com/taylorreiter/2022-infant-mge/issues/1). We ran `gather` with all databases and didn't get different results.)
Using genomes identified in the metagenome, we BLASTed each genome against the neighborhood query graph to determine whether that genome surrounded the MGE.

The strength of this approach is that is avoids assembly and binning, both of which [disproportionately exclude sequences from plasmids and genomic islands](https://doi.org/10.1099/mgen.0.000436).

## Results

![](https://i.imgur.com/u7UM6lM.png)
We applied `sourmash gather` to determine the taxonomic profile of the metagenome. 
We detected nine genomes representing seven species. 
This agrees with the original [publication of this data set](http://www.genome.org/cgi/doi/10.1101/gr.213256.116), which showed this community to be low diversity.
Sourmash likely missed viruses and plasmids that were present in the sample. 
See [notebooks/20220517_explore_taxonomic_composition_sourmash_gather.ipynb](https://github.com/taylorreiter/2022-infant-mge/blob/main/notebooks/20220517_explore_taxonomic_composition_sourmash_gather.ipynb).

We identified two interesting MGEs, both of which were antibiotic resistance genes.

The first encoded an *ErmB* genes, which is associated with the antibiotic [erythromycin](https://card.mcmaster.ca/ontology/36514).

![](https://i.imgur.com/W6HDbmk.png)

We detected this gene to be present in both *Enterococcus faecium* and *Clostridioides difficile* genomic backgrounds.

The second encoded *CfxA5*, which is associated with [beta lactamase resistance](https://card.mcmaster.ca/ontology/39649).
Unlike *Ermb*, we detected this only in *Bacteroides vulgatus*.
However, `sourmash gather` showed that the *B. vulgatus* genome we identified in this sample was 100% present, suggesting the exact MAG binned from this sample has been deposited in GenBank. 
Given that the entire contig was not a BLAST match for the *CfxA5* gene within the *B. vulgatus* bin, this suggests there is strain variation around this resistance element in the metagenome.

![](https://i.imgur.com/SDbm4b7.png)


## Limitations and future directions

This approach is primarily limited by the database used to search. 
If an MGE is not represented in the database, it cannot be identified in the metagenome.

In the long term, I would like to enable the [mobileOG-db](https://www.biorxiv.org/content/10.1101/2021.08.27.457951v1), however this database is currently only available as protein sequences. 
The code to execute protein searches with spacegraphcats [exists](https://github.com/spacegraphcats/spacegraphcats/pull/444), however I [haven't finished evaluating](https://github.com/taylorreiter/2021-cami-sgc-prot) this method for accuracy.

To save on time, I also used a scaled value of 2000 for the initial multifasta queries that identified candidate MGEs.
Spacegraphcats uses [FracMinHash sketches](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract) for some internal operations, specifically for `multifasta_queries`.
I generally use a lower scaled value such as 1000, which would increase the number of identifiable sequences.
For reference, in [the *Pseudomonas aeruginosa* genome, 2,617 of 5,892 ORFs were not represented in the FracMinHash sketch when using a scaled value of 2000](https://github.com/greenelab/2022-test-sketch-core-accessory-interactome/blob/main/notebooks/20220325_hashes_in_transcripts_and_operons.ipynb).
Decreasing the scaled value to 1000 would approximately decrease the missed sequences at least in half. 
This only impacts the initial multifasta query, and not the downstream investigation.

In the future, I hope to improve taxonomic recall of this workflow. 
Currently, the workflow identifies the genomic background of candidate MGEs using graph BLAST searches. 
Instead of this, I would like to implement a graph-based approach where shortest path traversals between a candidate MGE and single copy marker genes are used to infer the genomic context.
[Demo code to do this already exists](https://github.com/dib-lab/2020-ibd/blob/master/sandbox/test_corncob_dda/08_graph_explore.R), but still needs to be generalized and scaled.

Lastly, I would like to use graph structures themselves to discover new MGEs not represented in current databases.
This will require a better understanding of the graph structures of MGEs and other sequences.