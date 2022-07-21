# Battenberg

This repository contains code for the whole genome sequencing subclonal copy number caller Battenberg, as described in [Nik-Zainal, Van Loo, Wedge, et al. (2012), Cell](https://www.ncbi.nlm.nih.gov/pubmed/22608083).

## Package Update

Battenberg can now fully accommodate GRCh38 aligned data with the associated reference data generated with Beagle, including proper handling of chrX.

This repository contains the updates made to easily deploy the package for immediate use with GRCh38 aligned BAMs. Furthermore, this repository is fully incorporated into the [MGP1000](https://github.com/pblaney/mgp1000).

All original code and documents for the Battenberg R package were developed by the members of the [Wedge Lab](https://wedge-group.netlify.app/software/battenberg/) and the updates outlined in this README were incorporated by [Patrick Blaney](https://github.com/pblaney)

## Installation

### Dependencies

Battenberg depends on the following software:
* [HTSlib](https://github.com/samtools/htslib)
* [alleleCount](https://github.com/cancerit/alleleCount)
* [IMPUTE](https://mathgen.stats.ox.ac.uk/impute/impute.html)
* [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)

Additionally, the following R libraries:
* `devtools`
* `splines`
* `readr`
* `doParallel`
* `ggplot2`
* `RColorBrewer`
* `gridExtra`
* `gtools`
* `parallel`
* `igordot/copynumber`
* `GenomicRanges`
* `VanLoo-lab/ASCAT`
* `pblaney/battenberg`

While the package can be installed piecewise and run locally, the recommendation is to utilize the [Docker container](https://hub.docker.com/layers/mgp1000/patrickblaneynyu/mgp1000/battenberg-2.2.9/images/sha256-738003047651c426069f4e86bc21b383a4ba59bf83ac77f50e44b2ab3a9b766b?context=repo) created for this specific repository due to the complexity of the build, execution time, and computational resources needed.

Additionally, a [Singularity container](https://github.com/pblaney/mgp1000/tree/master/containers) is provided in the MGP1000 to deploy Battenberg within an HPC environment.

### Reference Files

This repository was developed to run exclusively with GRCh38 aligned BAMs. For users working with GRCh37 aligned BAMs, the [main Battenberg repository](https://github.com/Wedge-lab/battenberg) will work without need for adjustments.

While the Wedge Lab kindly generated the necessary reference files to be compatible with GRCh38 aligned BAMs, not all needed files can be found in a single place and their incorporation into the workflow is not straightforward. Therefore, an accessory script `battenberg_reference_downloader.sh` has been included here for users to easily download and prep the reference files for use with this repository's build.

The following code snippet will create the necessary directory structure and automate the download/prep steps.

**NOTE: The workflow expects the exact names of each directory as built in the code snippet, do not change the names. Additionally, this step can take up to a few hours due to the size of some files**
```
mkdir -p battenberg_reference
cd battenberg_reference/
mkdir -p GC_correction_hg38/
mkdir -p RT_correction_hg38/

./battenberg_reference_downloader.sh
```

## Workflow Execution

In order to simplify the execution of Battenberg, an accessory script `battenberg_executor.sh` has been included here for users. For proper execution, the script expects the `battenberg_reference/` directory compiled in the previous step to be present in the directory as this script, along with all other input files listed in the following code snippet.

```
battenberg_executor.sh \
[TUMOR_ID] \
[NORMAL_ID] \
[TUMOR_BAM] \
[NORMAL_BAM] \
[SEX_OF_SAMPLE] \
[OUTPUT_DIR] \
[NUM_OF_CPUS] \
[MIN_DEPTH_IN_NORMAL]
```

The `[SEX_OF_SAMPLE]` parameter should be expressed as either `male` or `female`.

The default value for `[MIN_DEPTH_IN_NORMAL]` is 10. However, this should be adjusted to a lower value depending on the average coverage for the input samples. It should be noted that this is an adaptation made specifically for this repository to utilize Battenberg with low coverage samples.

## Output Files

All files are deposited in the `[OUTPUT_DIR]` specified by the user.

* `[TUMOR_ID]_subclones.txt` contains the copy number data (see table below)
* `[TUMOR_ID]_rho_and_psi.txt` contains the purity estimate (make sure to use the FRAC_genome, rho field in the second row, first column)
* `[TUMOR_ID]_BattenbergProfile*png` shows the profile (the two variants show subclonal copy number in a different way)
* `[TUMOR_ID]_subclones_chr*.png` show detailed figures of the copy number calls per chromosome
* `[TUMOR_ID]_distance.png` This shows the purity and ploidy solution space and can be used to pick alternative solutions

The copy number profile saved in the `[TUMOR_ID]_subclones.txt` is a tab delimited file in text format. Within this file there is a line for each segment in the tumour genome.
Each segment will have either one or two copy number states:

* If there is one state that line represents the clonal copy number (i.e. all tumour cells have this state)
* If there are two states that line represents subclonal copy number (i.e. there are two populations of cells, each with a different state)

A copy number state consists of a major and a minor allele and their frequencies, which together add give the total copy number for that segment and an estimate fraction of tumour cells that carry each allele.

The following columns are available in the `[TUMOR_ID]_subclones.txt` output:

| Column | Description |
| ------------- | ------------- |
| chr | The chromosome of the segment |
| startpos | Start position on the chromosome |
| endpos | End position on the chromosome |
| BAF | The B-allele frequency of the segment |
| pval | P-value that is obtained when testing whether this segment should be represented by one or two states. A low p-value will result in the fitting of a second copy number state |
| LogR | The log ratio of normalised tumour coverage versus its matched normal sequencing sample |
| ntot | An internal total copy number value used to determine the priority of solutions. NOTE: This is not the total copy number of this segment! |
| nMaj1_A | The major allele copy number of state 1 from solution A |
| nMin1_A | The minor allele copy number of state 1 from solution A |
| frac1_A | Fraction of tumour cells carrying state 1 in solution A |
| nMaj2_A | The major allele copy number of state 2 from solution A. This value can be NA |
| nMin2_A | The minor allele copy number of state 2 from solution A. This value can be NA |
| frac2_A | Fraction of tumour cells carrying state 2 in solution A. This value can be NA |
| SDfrac_A | Standard deviation on the BAF of SNPs in this segment, can be used as a measure of uncertainty |
| SDfrac_A_BS | Bootstrapped standard deviation |
| frac1_A_0.025 | Associated 95% confidence interval of the bootstrap measure of uncertainty |

Followed by possible equivalent solutions B to F with the same columns as defined above for solution A (due to the way a profile is fit Battenberg can generate a series of equivalent solutions that are reported separately in the output).

### QC Plots

It also produces a number plots that show the raw data and are useful for QC.

* `[TUMOR_ID].tumour.png` and `[TUMOR_ID].germline.png` show the raw BAF and logR
* `[TUMOR_ID]_coverage.png` contains coverage divided by the mean coverage of both tumour and normal
* `[TUMOR_ID]_alleleratio.png` shows BAF * logR, a rough approximation of what the data looks like shortly before copy number calling
* `[TUMOR_ID]_chr*_heterozygousData.png` shows reconstructed haplotype blocks in the characteristic Battenberg cake pattern
* `[TUMOR_ID]_RAFseg_chr*.png` and `[TUMOR_ID]_segment_chr*.png` contains segmentation data for step 1 and step 2 respectively
* `[TUMOR_ID]_nonroundedprofile.png` shows the copy number profile without rounding to integers
* `[TUMOR_ID]_copynumberprofile.png` shows the copy number profile with (including subclonal copy number) rounding to integers
