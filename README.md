# mpra_starr_qc
QC pipeline for MPRA and STARR-seq experiments

### Scripts
- [create_STARR_counts_in_common_file_format.R](scripts/create_STARR_counts_in_common_file_format.R): Create fragment count files from STARR-seq BAM files.
- [bin_fragment.r](scripts/bin_fragment.r): Aggregate fragments in fixed size bins.
- [mpra_starr_qc.R](scripts/mpra_starr_qc.R): Main script to generate QC plots and summary statistics.
