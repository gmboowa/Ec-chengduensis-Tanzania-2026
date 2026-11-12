iTOL annotation bundle
======================

Load Ec_43_phylo_iqtree.nwk at https://itol.embl.de/upload.cgi, then drag each
annotation file onto the rendered tree. Legends are defined inside each file and
appear automatically.

  Ec_43_phylo_iqtree.nwk   maximum-likelihood tree, 43 genomes
  itol_country.txt         country of origin (15 categories)
  itol_source.txt          isolation source, grouped (9 categories)
  itol_amr_load.txt        acquired resistance genes: low <=5, medium 6-8, high >=9
  itol_binary.txt          plasmid replicon and ybtP/ybtQ presence
  itol_styles.txt          cluster shading; red bold labels for this study


The tree
--------

Inferred with IQ-TREE 2.4.0 under TVM+F+ASC from a Snippy 4.6.0 core-genome
alignment after recombination masking with Gubbins 3.2.1. Node support is from
1,000 ultrafast bootstrap replicates. The mapping reference has been removed
from this copy: SRR5247447 is a sequencing run of the same type strain and
differs from it at zero sites, so retaining both would duplicate a tip.

Tip labels are SRA run accessions, except A55822 and A55955, the isolates
sequenced for this study.


Display setting
---------------

Branch lengths differ by roughly three orders of magnitude between and within
the two clusters. On the default Phylogram setting the within-cluster structure
collapses to a vertical line. Either:

  - switch Display mode to "Normal" (topology only), or
  - collapse the upper cluster into a triangle (select the node, Collapse)

and give the divergence in the caption. The distance-distribution figure carries
the quantitative comparison, so the tree does not need to.


Isolation-source grouping
-------------------------

BioSample free-text values were grouped for legibility; the ungrouped values are
in tables/Table_S1_merged_genome_panel.xlsx.

  blood                        9     wound, skin or soft tissue   4
  stool or rectal swab         9     other sterile site           4
  urinary tract                5     not provided                 3
  environmental                5     respiratory                  2
                                     other clinical               2

"environmental" covers wastewater, soil, a farm and the walking-frame surface.
"not provided" means the original submitter did not populate the field; nothing
was inferred.


Cluster shading
---------------

itol_styles.txt shades each cluster using a pair of tips whose common ancestor
is the whole cluster. If iTOL shades more of the tree than expected, that pair
needs changing - and it would indicate the cluster is not monophyletic in the
displayed rooting, which is worth knowing.

  lower cluster  95.67-95.80% ANI to the type strain   n=22
  upper cluster  99.03-99.99%                          n=21

The two clusters differ by 121,375-122,118 core-genome SNPs. No pair of genomes
falls between 24,557 and 121,375 SNPs.
