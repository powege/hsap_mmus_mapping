#!/bin/bash

### Script that downloads ClinVar data 

# set file directory
PED_ROOT=/enter_file_directory_here/

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz \
-P "$PED_ROOT"

