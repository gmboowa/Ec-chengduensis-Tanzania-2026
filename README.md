# Genomic characterization of *Enterobacter chengduensis* from a Tanzanian hospital

Workflow configurations, parameters, accessions and analysis scripts supporting
the comparative genomic analysis of *Enterobacter chengduensis* isolates
recovered from a patient and from a hospital environmental surface in an
orthopaedic ward in Mwanza, Tanzania, set against all publicly deposited genomes
of the species that could be processed uniformly.

Everything required to reproduce the reported results is recorded here, with
every threshold and database version stated explicitly.

## Panel

A systematic search of the NCBI Sequence Read Archive on 11 November 2026 using

```
("Enterobacter chengduensis"[Organism] OR Enterobacter chengduensis[All Fields])
  AND "library layout paired"[Properties]
```

returned 57 runs. Two are the isolates deposited by this study (BioProject
PRJNA1447346), three are metagenomic rather than isolate libraries
(SRR35223178, SRR35223179, SRR35223195), and two (SRR37266573, SRR38980048)
were released after the analysis had been run. The remaining 50 public runs
together with the two study isolates give the **52 genomes assessed**, of which
**43 were retained**.

Inclusion required all four of: GAMBIT assignment to *E. chengduensis*; average
nucleotide identity to GCF_001984825.2 of at least 95%; fewer than 300 contigs;
estimated coverage of at least 30x. `accessions/excluded_9.txt` records the
reason for each exclusion.

## Software and database versions

| Component | Version |
|---|---|
| TheiaProk_Illumina_PE_PHB | PHB v4.3.0 |
| AMRFinderPlus | 4.2.7, database 2026-08-07.1 |
| PlasmidFinder | 2.1.6, bundled database, unmodified |
| GAMBIT | database 3.0.0-20260601 |
| GTDB-Tk | 2.7.2, GTDB release R232 |
| Shovill / SPAdes | 1.1.0 / 3.14.1 |
| BUSCO | 5.7.1, enterobacterales_odb10 |
| Snippy / snippy-core | 4.6.0 |
| Gubbins | 3.2.1 (RAxML 8.2.12 as internal tree builder) |
| IQ-TREE | 2.4.0 |
| snp-dists | 0.8.2 |

## Parameters

Gene detection was run at thresholds stricter than the tool defaults. Partial
alignments across contig boundaries in fragmented draft assemblies otherwise
produce calls that the underlying sequence does not support.

| Tool | Parameter | Tool default | Used here |
|---|---|---|---|
| AMRFinderPlus 4.2.7 | minimum identity | 0.90 | **0.95** |
| AMRFinderPlus 4.2.7 | minimum coverage | 0.50 | **0.80** |
| PlasmidFinder 2.1.6 | minimum identity | 0.90 | **0.95** |
| PlasmidFinder 2.1.6 | minimum coverage | 0.60 | **0.80** |

Average nucleotide identity was computed with MUMmer inside the TheiaProk
workflow against GCF_001984825.2 (*E. chengduensis* type strain WCHECl-C4,
assembly ASM198482v2) and is reported with the proportion of bases aligned.
Shovill capped input depth for assembly at 150x.

Phylogenetic analysis mapped quality-controlled reads to the same reference with
Snippy 4.6.0. Gubbins 3.2.1 ran five iterations with RAxML 8.2.12 under GTRGAMMA
and a maximum of 30% missing data per sequence (`--filter-percentage 30`). The
filter was set above the Gubbins default of 25% because one genome, DRR390471,
has 25.26% missing positions in the core-genome alignment and is dropped at the
default. That reduced alignability reflects divergence from the single
upper-band reference rather than any defect of the assembly, and the genome sits
5 SNPs from DRR390463, so all 43 retained genomes are present in the phylogeny
and in the SNP-distance matrix. ModelFinder selected TVM+F+ASC
under the Bayesian information criterion, and branch support came from 1,000
ultrafast bootstrap replicates.

## GTDB-Tk classification

All 52 assemblies were classified with GTDB-Tk 2.7.2 (`ani_rep`, GTDB release
R232), which identifies the closest GTDB species representative with skani and
reports whether a genome falls inside that species cluster's circumscription
radius.

| Closest GTDB species | Genomes |
|---|---|
| *Enterobacter chengduensis* | 47 |
| *Enterobacter kobei* | 3 |
| *Enterobacter asburiae* | 1 |
| *Enterobacter* sp913778545 | 1 |

The 47 *E. chengduensis* genomes all have GCF_001984825.2 as their closest
representative and all satisfy the GTDB circumscription criteria, including the
22 in the lower of the two ANI bands at 95.2-95.8% skani ANI. The two bands are
therefore structure within the species as currently circumscribed, not separate
species clusters. The five genomes placed outside *E. chengduensis* are the same
five excluded by the GAMBIT and ANI criteria, so GTDB-Tk confirms those
exclusions by an independent comparison.

Per-genome results are in `results/gtdbtk/`.

## Contents

```
workflows/   GTDBTK_ANI_REP.wdl and its Terra inputs/outputs JSON
tables/      merged supplementary table; BioSample metadata; virulence tabulation
accessions/  the retained panel, the exclusions with reasons, Terra set membership
results/     GTDB-Tk per-genome classification and species counts; phylogeny
             outputs, SNP-distance matrix, selected model and tool versions
itol/        the 43-tip tree and the iTOL annotation files for Figure 1
phylogeny/   tree files, with and without the mapping reference
scripts/     BioSample parsing, supplementary table construction, the GTDB
             column merge, virulence tabulation, phylogeny drawing, iTOL
             annotation build and checker, and a static checker run over each
             WDL before import. Paths are relative to the repository root.
```

`tables/Table_S1_merged_genome_panel.xlsx` is the merged supplementary table:
run and BioSample accession, strain, origin, geography, isolation source, host
and collection date as recorded in NCBI BioSample, assembly metrics, ANI, GAMBIT
assignment, GTDB-Tk classification, resistance determinants, replicon markers, and inclusion status.
Blank metadata cells are attributes the original submitter did not populate;
they were not inferred.

`scripts/draw_tree.py` regenerates the phylogeny figure from a Newick file and
the supplementary table using only matplotlib, so it runs without a
phylogenetics package installed.

## Reuse

The workflow files and scripts are released under the MIT Licence; see
`LICENSE`. Genome accessions and their metadata belong to the groups that
deposited them; please cite the original submitters where their data are used.

## Citation

Mboowa G, Kidenya BR, Sserwadda I, Kanyerezi S, Seni J. Genomic characterization
of patient-associated and abiotic-surface *Enterobacter chengduensis* isolates
from a Tanzanian hospital in a global comparative context. Submitted.
