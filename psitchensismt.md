---
title: "Largest Complete Mitochondrial Genome, Assembled with 10x Genomics Chromium Linked Reads: Sitka Spruce (*Picea sitchensis*)"
author: [Shaun D. Jackman, Lauren Coombe, René L. Warren, Chen Yang, Benjamin P. Vandervalk, Richard A. Moore, Stephen Pleasance, Robin J. Coope, Joerg Bohlmann, Robert A. Holt, Steven J. M. Jones, Inanc Birol]
bibliography: psitchensismt.bib
csl: psitchensismt.csl
rangeDelim: "&ndash;"
eqnPrefix: "Equation"
figPrefix: "Fig."
tblPrefix: ["Table", "Tables"]
keywords: [genome, assembly, plant, gymnosperm, conifer, spruce, organelle, mitochondrion, 10x, chromium, ABySS]
---

# Abstract

Long-read sequencing technologies have greatly improved assembly contiguity, but at a cost roughly ten times that of short-read sequencing technology. For population studies and when sequencing large genomes, such as conifer genomes and other economically important crop species, this cost may be prohibitive. The 10x Genomics Chromium technology generates linked reads from large molecules at a cost comparable to standard short-read sequencing technologies. Whereas paired-end sequencing gives two reads from a small DNA fragment, Chromium yields roughly a hundred reads from molecules with a typical size of 10 to 100 kilobases. This technology has been applied to assemble genomes as large as the human genome.

Here we demonstrate its utility to localize assemblies of genomic loci, in particular isolating the organellar genomes from the 20 gigabase Sitka spruce (*Picea sitchensis*) nuclear genome. We have reported earlier the assembly of the plastid genome of Sitka spruce using the Gemcode technology, the precursor of Chromium. Using the Chromium technology, we have now assembled the mitochondrial genome of Sitka spruce. Our preliminary draft assembles 75% of the estimated six megabase mitochondrial genome in 15 scaffolds larger than 100 kbp.

Animal mitochondrial genomes are typically quite small, tens of kilobases, and plant mitochondria are often much larger, hundreds of kilobases. Conifer mitochondria are, in contrast, surprisingly large, nearly six megabases for white spruce (*Picea glauca*). With a mere two gymnosperm mitochondria in NCBI Genbank and a single conifer mitochondrion, few examples of these mysteriously large organellar genomes are available. Our results will thus be of interest to evolutionary biology researchers beyond the field of forestry genomics.

Chromium reads permit cost-effective assembly of large genomes with high-throughput, short-read sequencing technology, while also providing large-molecule scaffolding data. Here we assemble Chromium reads using ABySS 2.0, which uses one-tenth the memory of previous versions, and only a few hundred gigabytes of RAM for large conifer genomes.

# Introduction

+ White spruce (*Picea sitchensis*) organellar genomes [@Jackman_2015]
+ Sitka spruce (*Picea glauca*) plastid genome [@Coombe_2016]

# Methods

+ ABySS 2.0 [@Jackman_2017]
+ ARCS [@Yeo_2017]

# References
