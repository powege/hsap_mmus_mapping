#!/bin/bash

### Script that downloads and unzips GWAS catalog

# set file directory
PED_ROOT=/enter_file_directory_here/

# set date 
todays_date=$(date +"%Y_%m_%d")

### associations
wget https://www.ebi.ac.uk/gwas/api/search/downloads/full \
-P "$PED_ROOT"
mv "$PED_ROOT"full "$PED_ROOT"gwas_catalog_associations_"$todays_date".tsv

### ontology index
wget https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings \
-P "$PED_ROOT"
mv "$PED_ROOT"trait_mappings "$PED_ROOT"gwas_catalog_trait_mappings_"$todays_date".tsv
