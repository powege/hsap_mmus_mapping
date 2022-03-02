#!/bin/bash

########
### SET VARS
########

# set working directory for files to download into
PED_ROOT=/enter_file_path_here/
cd $PED_ROOT

########
### download and unzip pairwise alignment
########


# human mouse
wget ftp://ftp.ensembl.org/pub/release-101/maf/ensembl-compara/pairwise_alignments/hsap_grch38.v.mmus_grcm38.lastz_net.tar.gz \
-P "$PED_ROOT"

# mouse mus caroli
wget ftp://ftp.ensembl.org/pub/release-101/maf/ensembl-compara/pairwise_alignments/mmus_grcm38.v.mcar_caroli_eij_v1.1.lastz_net.tar.gz \
-P "$PED_ROOT"

# mouse mus pahari
wget ftp://ftp.ensembl.org/pub/release-101/maf/ensembl-compara/pairwise_alignments/mmus_grcm38.v.mpah_pahari_eij_v1.1.lastz_net.tar.gz \
-P "$PED_ROOT"

### UNZIP 
tar xvzf "$PED_ROOT"hsap_grch38.v.mmus_grcm38.lastz_net.tar.gz
tar xvzf "$PED_ROOT"mmus_grcm38.v.mcar_caroli_eij_v1.1.lastz_net.tar.gz
tar xvzf "$PED_ROOT"mmus_grcm38.v.mpah_pahari_eij_v1.1.lastz_net.tar.gz




