#!/bin/bash -l

#### Ideally you wouold only need to execute this code once per new R installation from current directory ####
#### The current directory should be home/user ####
#### The ouput is a file pkg_list with all the packages needing installation ####
#### The output file will be read by install_pkgs.R #### 

# Load R version -- modified as needed --
module load bioinfo-tools
module load R/3.6.0
module load R_packages/3.6.0

<<'_//comment//_' 
#### This desn't work for bianca because it doesn't have internet connection ####
date_time=$(date +'%Y-%m-%d_%T')

# Create R libraries directory -- modify as needed, version directory should be the same as the R version loaded --
args1="/castor/project/home/mararc/bin/R/4.0.0"
if [ ! -d $args1 ]; then mkdir -p $args1; fi

# Get all the libraries that are used in my R codes
grep -r -o -h --exclude-dir={routs.d,submitters.d} 'library(.*)' ./*/code.d | sort -u | grep -o  '".*"' | tr -d \" > tmp_pkgs 

# Install packages
Rscript --vanilla install_pkgs.R $args1 > install_pkgs_${date_time}.ROUT

# Save pkgs installed
mv tmp_pkgs Rpkg_list_${date_time}
-//comment_//
