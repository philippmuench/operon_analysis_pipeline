# Operon Conservation Analysis Methods

## Genome Annotation

Prior to homology searches, all bacterial genome assemblies were annotated using Prokka version 1.14.6 to ensure consistent gene prediction and annotation across the dataset. We processed 8,587 *Enterococcus faecalis* genome assemblies obtained from public databases, with each genome annotated using species-specific parameters (--kingdom Bacteria --genus Enterococcus --species faecalis).

## BLAST-based Homology Search

We employed a multi-strategy BLAST search approach to comprehensively identify operon homologs across bacterial genomes. Reference sequences for the seven-gene operon (frpC, glpC, ptsD, ptsC, ptsB, ptsA, fruR) were extracted from the *Enterococcus faecalis* V583 genome and used as queries in four complementary search strategies.

The reference query set comprised 7 protein sequences and 9 nucleotide sequences (7 coding sequences plus 2 non-coding regulatory elements including the promoter region). BLAST searches were performed against 8,287 prokaryotic genomes using NCBI BLAST+ version 2.12.0 with default parameters except where noted.

For protein-based searches, we used tblastn to align the 7 reference protein sequences against genomic nucleotide sequences (.fna files) from Prokka-annotated genomes. This approach enabled detection of homologous genes regardless of annotation quality, with an E-value threshold of 1×10⁻⁵. To validate findings and capture nucleotide-level conservation patterns, we performed complementary blastn searches using the corresponding nucleotide sequences as queries against the same genomic targets.

To generate sequence variants suitable for phylogenetic analysis, we implemented a third strategy where Prokka-predicted coding sequences (.ffn files) from each genome were searched against a nucleotide database constructed from our reference operon genes using blastn with a minimum identity threshold of 80%. This approach captured query sequence fragments (qseq) that reflect the actual sequence diversity present in the target genomes.

Finally, we searched for non-coding regulatory elements using blastn to align promoter and regulatory sequences against genomic nucleotide sequences, enabling identification of conserved regulatory motifs upstream of operon genes.

Across all search strategies, we generated 473,867 BLAST hits, of which 55,826 (15.8%) met high-quality criteria (≥90% sequence identity and ≥80% query coverage). The mean sequence identity across all hits was 43.4% with 94.5% mean query coverage, indicating broad phylogenetic distribution of operon homologs with varying degrees of conservation.

For each genome and query gene combination, we retained only the highest-scoring hit based on bitscore to avoid redundancy. Both DNA strands were searched automatically by BLAST, with reverse complement matches identified by coordinate order (sstart > send). This comprehensive search strategy ensured robust detection of operon components across diverse bacterial lineages while capturing the sequence variation necessary for downstream phylogenetic and conservation analyses.

Multiple sequence alignments were generated using MAFFT version 7.490 with automatic algorithm selection, and temporary files were directed to high-performance local storage (/vol/tmp) to optimize computational efficiency during large-scale alignment operations.
