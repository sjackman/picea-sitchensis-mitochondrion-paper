---
title: "Largest Complete Mitochondrial Genome of a Gymnosperm, Sitka Spruce (*Picea sitchensis*), Assembled using Long Reads"
author: [Shaun D. Jackman, Lauren Coombe, René L. Warren, Chen Yang, Benjamin P. Vandervalk, Richard A. Moore, Stephen Pleasance, Robin J. Coope, Joerg Bohlmann, Robert A. Holt, Steven J. M. Jones, Inanc Birol]
bibliography: psitchensismt.bib
csl: psitchensismt.csl
rangeDelim: "&ndash;"
eqnPrefix: "Equation"
figPrefix: "Fig."
tblPrefix: ["Table", "Tables"]
keywords: [mitochondrion, genome, assembly, plant, gymnosperm, conifer, spruce, organelle, Oxford, Nanopore, Illumina, 10x Genomics, Chromium, linked reads]
---

# Abstract

Plant mitochondrial genomes vary widely in size. Although many plant mitochondrial genomes have been sequenced and assembled, the vast majority are of angiosperm, and few are of gymnosperm. Most plant mitochondrial genomes are smaller than a megabase, with a few notable exceptions. We have sequenced and assembled the 5.5 Mbp mitochondrial genome of Sitka spruce (*Picea sitchensis*), the largest complete mitochondrial genome of a gymnosperm. We sequenced the whole genome using Oxford Nanopore MinION and identified contigs of mitochondrial origin assembled from these long reads. The assembly graph shows a multipartite genome structure, composed of one smaller 168 kbp circular segment of DNA, and a larger 5.4 Mbp component with a branching structure. The branching points of the assembly graph may represent active sites of recombination and give insight into an underlying complex linear branching physical genome structure.

# Introduction

Annotating plant mitochondrial genomes is made difficult by numerous features of plant mitochondria that are not typical of most organisms. RNA editing of C-to-U is pervasive and creates AUG start codons by editing ACG to AUG and can also create stop codons in a similar fashion [@Hiesel_1989]. An alternative GUG start codon is used by some genes [@Sakamoto_1997], which may be created by RNA editing of a GCG codon. The typical GU-AG splice site expected by most splice-aware alignments tools is instead GUGCG-AY (Y denotes C or T) for group II introns [@Lambowitz_2010]. Trans-spliced genes are common in mitochondrial genomes [@Kamikawa_2016], and no purpose-built software tool exists for identifying and annotating trans-spliced genes. For these reasons, annotating a plant mitochondrial genome remains a laborious and manual task.

# Methods

## Genome sequencing and assembly

Genomic DNA was extracted from young needles. Using 18 Oxford Nanopore MinION flow cells, we sequenced 98 Gbp in 9.6 million reads of whole genome sequencing, which yields 5 fold depth of the roughly 20 Gbp nuclear genome and 26 fold depth of coverage of the mitochondrial genome. Separating putative mitochondrial reads by homology to known mitochondrial sequences could discard mitochondrial sequences that are unique Sitka spruce. We instead chose to first assemble the whole genome reads and then compare contigs to known mitochondrial sequences. Assembling such a large number of Nanopore reads is not yet straight forward, and so we adopted an iterative approach to assembly. We first obtained a rough but computationally efficient assembly using Miniasm [@Li_2016], after trimming adapter sequences with Porechop [@Wick_2017_Porechop]. Miniasm produces an assembly whose sequencing error rate is comparable to that of the original reads, but no better. We polished this assembly using Racon [@Vaser_2017]. We selected contigs with homology to the white spruce (*Picea glauca* isolate PG29) mitochondrial genome [@Jackman_2015] using Bandage [@Wick_2015], retaining contigs with at least one 5 kbp alignment to the white spruce mitochondrion by BLASTN [@Altschul_1990]. This selection process would also discard any mitochondrial plasmids smaller than 5 kbp.

We select putative mitochondrial reads by aligning the Nanopore reads to this assembly using Minimap2 [@Li_2018] and retaining reads with an alignment score of 5000 or more. We assembled these reads using Unicycler [@Wick_2017_Unicycler]. This assembly yielded one circular contig and many linear contigs with no adjacent contigs, indicating that the assembly is not yet complete. We repeated the alignment of the Nanopore reads to the assembly and again retained reads with an alignment score of 5000 or more. We assembled these reads using Flye [@Kolmogorov_2018], taking the assembly graph `2-repeat/graph_final.gfa`, which identifies repeats that are longer than the read length and determines their precise boundaries. This Flye assembly was polished using Racon. Contigs with homology to the white spruce mitochondrion (alignment length at least 5000 and percent identity at least 90) were selected using Bandage. Contigs with unambiguous adjacent contigs were merged using the Bandage operation "Merge all possible nodes".

Nanopore reads have difficulty accurately determining the length of homopolymer repeats. We polished the assembly using one flow cell of Illumina HiSeq sequencing reads of the same DNA extraction, yielding 59 fold depth of coverage of the mitochondrion, to correct homopolymer errors. We used Unicycler Polish to iteratively align the reads to the assembly using Bowtie2 [@Langmead_2012] and correct the consensus sequence using Pilon [@Walker_2014]. Ten rounds of polishing yielded no further corrections on the tenth round. Unicycler Polish applies Assembly Likelihood Estimate (ALE) [@Clark_2013] to each round to verify that the assembly of the final round of polishing most resembles the reads.

## Annotation

We annotated coding genes and non-coding rRNA and tRNA genes using automated methods where possible, and manual inspection to refine these automated annotations. We used Prokka [@Seemann_2014], which uses Prodigal [@Hyatt_2010], to identify single-exon coding genes and open reading frames (ORFs). We used MAKER [@Holt_2011], which uses BLASTP and Exonerate [@Slater_2005], to identify cis-spliced coding genes. We used tRNAscan-SE [@Lowe_1997] and Aragorn [@Laslett_2004] to identify tRNA. We used RNAmmer [@Lagesen_2007] and Barrnap [@Seemann_2014] to identify rRNA. We used RNAweasel [@Lang_2007] to identify domain V of group II introns.

Following automated annotation, we reviewed coding genes for completeness, compared to their best BLASTP match and corrected the annotation, most often for aspects that are particular to plant mitochondria. We manually corrected the annotation of genes to address start codons created by RNA editing ACG to the start codon AUG and editing GCG to the alternative start codon GUG (see results for details). Three genes use atypical start codons: *rpl16* uses a GUG start codon [@Sakamoto_1997]; *rps19* uses a GUG start codon created by RNA editing GCG, seen also in *Pinus strobus* AJP33554.1; *matR* appears to use an unusual GGG start codon, seen also in *Cycas taitungensis* YP_001661429.1 [@Chaw_2008] and *Pinus strobus* AJP33535.1. Then gene *sdh4* was missed by automatic annotation, as its coding sequence was found to overlap on the same strand that of *cox3* by 73 bp.

We reviewed splice sites and adjusted their position to agree with the expected splicing motifs of group II introns when possible, ensuring not to introduce insertions or deletions into the peptide sequence compared to homologous proteins. We confirmed the presence of domain V of the group II intron upstream of the 3' splice site, identified by RNAweasel. We manually annotated trans-spliced introns by comparing alignments of homologous proteins to the genome. We determined the 5' and 3' splice sites similarly to cis-spliced introns, looking for expected group II splicing motifs, and domain V upstream of the 3' splice site.

# Results

## Genome assembly

The mitochondrial genome of Sitka spruce is 5.52 Mbp assembled in 13 segments, with 44.7% GC content. The genome assembly is composed of two components: a 168 kbp circular segment, and a larger 5.36 Mbp component composed of 12 segments, visualized by Bandage [@Wick_2015] in @fig:genome. The two smallest segments have an estimated copy number of 2 based on their depth of sequencing coverage, and all other segments have similar depth of coverage. The single-copy segments range in size from 84 kbp to 1.65 Mbp. The two repeat segments are 24 kbp and 27 kbp. No sequence variation is evident in these repeats. An absence of variation in the repeat implies that they may be involved in active recombination. Though 10% of reads are larger than 24 kbp, no reads fully span these repeats.

![The assembly graph of the mitochondrial genome of Sitka spruce. Each segment is labeled with its size and named 01 through 13 by rank of size.](figures/genome.png){#fig:genome}

## Genome annotation

The mitochondrial genome of Sitka spruce has 41 distinct protein coding genes with known function, 3 distinct rRNA genes, 26 distinct tRNA genes representing 17 distinct anticodons. The 41 protein coding genes are those seen in the gymnosperm mitochondria *Cycas taitungensis* and *Ginkgo biloba* [@Guo_2016]. 27 introns are found in 10 protein coding genes, 14 cis-spliced and 13 trans-spliced, shown in @tbl:introns.

The 5S rRNA gene *rrn5* is present in four copies. The small subunit rRNA gene *rrn18* is present in three copies, though one copy is found on the 27 kbp repeat segment with an estimated copy number of two, yielding a total copy number of four. One copy of the large subunit rRNA gene *rrn26* is present, though it is found on the 24 kbp repeat segment, which has an estimated copy number of two.

26 tRNA genes are seen, representing 17 distinct anticodons, coding for 13 distinct amino acids, DEHKLMNPQRTWY. *trnM-CAU* has seven copies, *trnD-GUC* three copies, and *trnY-GUA* two copies. All other tRNA genes are single copy. No introns are observed in tRNA. This account does not consider possible RNA editing of the anticodon.

Table: 27 introns are seen in ten protein coding genes, 14 cis-spliced and 13 trans-spliced. \*A second partial copy of *nad5* has one cis-spliced intron. {#tbl:introns}

| Gene  | Sitka introns | Cis | Trans | Cycad introns | Cis | Trans |
|-------|--------------:|----:|------:|--------------:|----:|------:|
| ccmFC |             1 |   1 |     0 |             1 |   1 |     0 |
| cox2  |             2 |   0 |     2 |             2 |   2 |     0 |
| nad1  |             4 |   1 |     3 |             4 |   2 |     2 |
| nad2  |             4 |   2 |     2 |             4 |   3 |     1 |
| nad4  |             3 |   2 |     1 |             3 |   3 |     0 |
| nad5  |             4 |   2 |     2 |             4 |   2 |     2 |
| nad5\*|             1 |   0 |     0 |               |     |       |
| nad7  |             4 |   2 |     2 |             4 |   4 |     0 |
| rpl2  |             1 |   0 |     1 |             1 |   1 |     0 |
| rps3  |             2 |   2 |     0 |             2 |   2 |     0 |
| rps10 |             1 |   1 |     0 |             1 |   1 |     0 |
| Total |            27 |  14 |    13 |            26 |  21 |     5 |

# References
