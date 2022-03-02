#!/bin/bash

########
### SET VARS
########

# set working directory for files to download into
PED_ROOT=/enter_file_directory_here/

########
### download and unzip gene annotation
########

wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz \
-P "$PED_ROOT"
wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz \
-P "$PED_ROOT"

########
### download and unzip Regulatory Build annotation
########

# All cell types
wget ftp://ftp.ensembl.org/pub/release-101/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz \
-P "$PED_ROOT"
wget ftp://ftp.ensembl.org/pub/release-101/regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz \
-P "$PED_ROOT"

# heart 
wget ftp://ftp.ensembl.org/pub/release-101/regulation/homo_sapiens/RegulatoryFeatureActivity/heart/homo_sapiens.GRCh38.heart.Regulatory_Build.regulatory_activity.20190329.gff.gz \
-P "$PED_ROOT"
wget ftp://ftp.ensembl.org/pub/release-101/regulation/mus_musculus/RegulatoryFeatureActivity/heart_adult_8_weeks/mus_musculus.GRCm38.heart_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz \
-P "$PED_ROOT"

# kidney
wget ftp://ftp.ensembl.org/pub/release-101/regulation/homo_sapiens/RegulatoryFeatureActivity/kidney/homo_sapiens.GRCh38.kidney.Regulatory_Build.regulatory_activity.20190329.gff.gz \
-P "$PED_ROOT"
wget ftp://ftp.ensembl.org/pub/release-101/regulation/mus_musculus/RegulatoryFeatureActivity/kidney_adult_8_weeks/mus_musculus.GRCm38.kidney_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz \
-P "$PED_ROOT"

# spleen
wget ftp://ftp.ensembl.org/pub/release-101/regulation/homo_sapiens/RegulatoryFeatureActivity/spleen/homo_sapiens.GRCh38.spleen.Regulatory_Build.regulatory_activity.20190329.gff.gz \
-P "$PED_ROOT"
wget ftp://ftp.ensembl.org/pub/release-101/regulation/mus_musculus/RegulatoryFeatureActivity/spleen_adult_8_weeks/mus_musculus.GRCm38.spleen_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz \
-P "$PED_ROOT"

# thymus
wget ftp://ftp.ensembl.org/pub/release-101/regulation/homo_sapiens/RegulatoryFeatureActivity/thymus_1/homo_sapiens.GRCh38.thymus_1.Regulatory_Build.regulatory_activity.20190329.gff.gz \
-P "$PED_ROOT"
wget ftp://ftp.ensembl.org/pub/release-101/regulation/mus_musculus/RegulatoryFeatureActivity/thymus_adult_8_weeks/mus_musculus.GRCm38.thymus_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz \
-P "$PED_ROOT"
