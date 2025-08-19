# Operon Conservation Analysis Methods

## Genome Annotation

Prior to homology searches, all bacterial genome assemblies were annotated using Prokka version 1.14.6 to ensure consistent gene prediction and annotation across the dataset. We processed 8,587 *Enterococcus faecalis* genome assemblies obtained from public databases, with each genome annotated using species-specific parameters (--kingdom Bacteria --genus Enterococcus --species faecalis).

## Reference Operon Sequence Extraction

Reference sequences for homology searches were extracted from a curated GenBank file containing the complete fructoselysine/glucoselysine utilization operon from *Enterococcus faecalis* V583. This seven-gene operon encodes a complete pathway for the utilization of advanced glycation end products through a phosphoenolpyruvate:carbohydrate phosphotransferase system (PTS).

The operon consists of seven protein-coding genes arranged in the following order: *frpC* (fructoselysine-6-phosphate deglycase), *glpC* (glucoselysine-6-phosphate deglycase), *ptsD* (PTS system IID component), *ptsC* (PTS system IIC component), *ptsB* (PTS system IIB component), *ptsA* (PTS system IIA component), and *fruR* (sigma-54 dependent transcriptional regulator). Additionally, two regulatory elements were identified: the operon promoter region and the Pribnow box (-10 consensus sequence).

Sequence extraction was performed using BioPython to parse the GenBank format and extract both nucleotide and amino acid sequences. All genes in this operon are encoded on the complement strand and were appropriately reverse-complemented during extraction. The extraction process generated multiple query sequence sets: 7 protein sequences (average length 352 amino acids, range 154-919 residues), 7 nucleotide gene sequences (average length 1,058 bp, range 465-2,760 bp), and 2 regulatory nucleotide sequences for non-coding elements.

These reference sequences were formatted as FASTA files with standardized headers containing locus tags, gene names, and functional annotations. The resulting query sets provided comprehensive coverage of both coding and regulatory elements necessary for identifying complete operon homologs across diverse bacterial genomes in subsequent BLAST searches.

## BLAST-based Homology Search

We employed a multi-strategy BLAST search approach to comprehensively identify operon homologs across bacterial genomes. Reference sequences for the seven-gene operon (frpC, glpC, ptsD, ptsC, ptsB, ptsA, fruR) were extracted from the *Enterococcus faecalis* V583 genome and used as queries in four complementary search strategies.

The reference query set comprised 7 protein sequences and 9 nucleotide sequences (7 coding sequences plus 2 non-coding regulatory elements including the promoter region). BLAST searches were performed against 8,287 prokaryotic genomes using NCBI BLAST+ version 2.12.0 with default parameters except where noted.

For protein-based searches, we used tblastn to align the 7 reference protein sequences against genomic nucleotide sequences (.fna files) from Prokka-annotated genomes. This approach enabled detection of homologous genes regardless of annotation quality, with an E-value threshold of 1×10⁻⁵. To validate findings and capture nucleotide-level conservation patterns, we performed complementary blastn searches using the corresponding nucleotide sequences as queries against the same genomic targets.

To generate sequence variants suitable for phylogenetic analysis, we implemented a third strategy where Prokka-predicted coding sequences (.ffn files) from each genome were searched against a nucleotide database constructed from our reference operon genes using blastn with a minimum identity threshold of 80%. This approach captured query sequence fragments (qseq) that reflect the actual sequence diversity present in the target genomes.

Finally, we searched for non-coding regulatory elements using blastn to align promoter and regulatory sequences against genomic nucleotide sequences, enabling identification of conserved regulatory motifs upstream of operon genes.

Across all search strategies, we generated 473,867 BLAST hits, of which 55,826 (15.8%) met high-quality criteria (≥90% sequence identity and ≥80% query coverage). The mean sequence identity across all hits was 43.4% with 94.5% mean query coverage, indicating broad phylogenetic distribution of operon homologs with varying degrees of conservation.

For each genome and query gene combination, we retained only the highest-scoring hit based on bitscore to avoid redundancy. Both DNA strands were searched automatically by BLAST, with reverse complement matches identified by coordinate order (sstart > send). This comprehensive search strategy ensured robust detection of operon components across diverse bacterial lineages while capturing the sequence variation necessary for downstream phylogenetic and conservation analyses.

Multiple sequence alignments were generated using MAFFT version 7.490 with automatic algorithm selection, and temporary files were directed to high-performance local storage (/vol/tmp) to optimize computational efficiency during large-scale alignment operations.

## Core Gene Analysis

To establish a baseline for evaluating operon gene conservation, we performed comprehensive core gene analysis across all annotated *E. faecalis* genomes. Core genes were defined as those present in ≥95% of genomes, representing the conserved genetic backbone of the species that could serve as a reference for comparison with the more variable operon genes.

Core gene identification was performed by parsing all Prokka-generated GFF annotation files to extract gene names and calculate prevalence statistics across the 8,587 genome dataset. This analysis identified 5,540 unique genes across all genomes, with 1,251 genes meeting the ≥95% prevalence threshold for core gene designation. Notably, no genes were present in 100% of genomes, reflecting the inherent genomic diversity within *E. faecalis*, while 1,071 genes achieved ≥99% prevalence, representing the most highly conserved genetic elements.

For each core gene, nucleotide sequences were extracted from Prokka annotations across all genomes where the gene was present, yielding an average of 7,693 sequences per gene (range dependent on gene prevalence). The extracted sequences averaged 987 bp in length (range 171-2,340 bp), reflecting the diversity of core gene functions from small regulatory genes to large metabolic enzymes.

Multiple sequence alignments were generated for all 1,280 core genes using MAFFT with automatic algorithm selection and parallel processing for computational efficiency. The resulting alignments averaged 1,095 positions in length and contained an average of 7,709 sequences per alignment, providing robust datasets for conservation analysis.

Conservation analysis was performed using Shannon entropy-based scoring, where conservation scores range from 0 (highly variable) to 1 (perfectly conserved). The core genes exhibited high overall conservation with a mean score of 0.931, confirming their fundamental importance to *E. faecalis* biology. Conservation analysis revealed that 935 genes (73.0%) were highly conserved (score ≥0.9), 273 genes (21.3%) were moderately conserved (0.7-0.9), and only 72 genes (5.6%) showed lower conservation scores. Mean pairwise sequence identity across core genes was 0.924, with an average gap percentage of 8.2% in alignments.

This comprehensive core gene analysis established the baseline conservation patterns for essential *E. faecalis* genes, against which the conservation of the fructoselysine/glucoselysine operon genes could be compared to assess their relative evolutionary constraint and functional importance within the species.
