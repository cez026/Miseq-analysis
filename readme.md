## MiSeq Analysis for HIV Deep Sequencing Paper

The scripts contained in this repository were used in part to produce the
results found in ["Deep sequencing of protease inhibitor resistance HIV patient
isolates reveals patterns of correlated mutations in gag and protease"][1].

### Revisions

Most of these scripts were written between 2013-2014, and the manuscript
preparation and review took most of 2014, so these scripts are a bit outdated.
I will attempt to clean them up over time, but wanted to commit to originals so
there is a public review history.  The majority of the scripts need:

-   File paths cleaned up
-   Commented-out/unused/broken methods revised or removed
-   Directory structure fixed as they did not previously exist in a single
    directory

### Where is the original data?

The raw sequencing data (fastq), aligned and cleaned files (bam/sam), and
snp/amino acid substition databases are excluded from this repository due to the
private nature of patient information.  Additional information regarding IRB
documentation or data sources can be found in [published manuscript][1].  The
scripts herein have been cleaned of any analysis methods that invoke specific
patient identifiers (anonymized ids and sample numbers) so that nothing is
revealed.

[1] http://dx.doi.org/10.1371/journal.pcbi.1004249
