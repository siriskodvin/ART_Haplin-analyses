README

This folder contains script files for the Haplin analyses in ART triads.
The software used are PLINK and R.

Script descriptions:
FRchr.sh is the master shell script
FRhap.R is the main analysis script
FRmpl.R merges and plots results
FRhaplo.R is the haplotype analysis script
FRmplh.R merges and plots the haplotype results

Input files required are binary genetic files (.bim, .bed, and .fam).
Also required is a key file linking the family triads/dyads.
Wherever an input file is expected in the code it is denoted by FILE.filetype, e.g. GENETICS.bim.
Paths are denoted by /PATH_DESCRIPTION/, e.g. /PATH_TO_SCRIPT_FILES/.

All code is written by Siri Nærland Skodvin (2024).
