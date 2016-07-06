## MiSeq Analysis for HIV Deep Sequencing Paper

The scripts contained in this repository were used in part to produce some or
all of the results found in ["Deep sequencing of protease inhibitor resistance
HIV patient isolates reveals patterns of correlated mutations in gag and
protease"][1].

### Revisions

Most of these scripts were written between 2013-2014, and the manuscript
preparation and review took most of 2014, so these scripts are a bit outdated.

Also note that there is (a) no guarantee that the output of these scripts ended
up in the manuscript or (b) these are the correct/final versions of the scripts
because there are several unversioned copies flying around in different places
and I snagged these versions from a backup made around 2014/04/04.  That said,
I will attempt to clean them up over time, but wanted to commit (what may be)
the originals so there is a public review history.  The majority of the scripts
need:

-   File paths cleaned up
-   Commented-out/unused/broken methods revised or removed
-   Directory structure fixed as they did not previously exist in a single
    directory

### Where is the original data?

The raw sequencing data (fastq), aligned and cleaned files (bam/sam), and
SNP/amino acid substitution databases are excluded from this repository due to the
private nature of patient information.  Additional information regarding IRB
documentation or data sources can be found in [published manuscript][1].  The
scripts herein have been cleaned of any analysis methods that invoke specific
patient identifiers (anonymized IDs and sample numbers) so that nothing is
revealed.

[1]: http://dx.doi.org/10.1371/journal.pcbi.1004249
