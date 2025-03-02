#!/bin/bash

# Manual download from Ensemble mart
# 20240823

# Output 1: mart_export_GRCh38p14.txt
# Output 2: mart_export_GRCh37p13.txt

# http://mart.ensembl.org/biomart/martview/9590404cf4a6dbb61097944d3a06b7eb
# https://grch37.ensembl.org/biomart/martview/51d60049917fc0a40d6f416fb6fddefb
# Dataset
# Human genes (GRCh38.p14)
# Filters
# [None selected]
# Attributes
# Gene stable ID
# Gene stable ID version
# Gene start (bp)
# Gene end (bp)
# Chromosome/scaffold name
# Gene name
# MIM morbid description

cat mart_export_GRCh38p14.txt | cut -f1,2,3,4 | uniq > mart_export_GRCh38p14_short.txt
cat mart_export_GRCh37p13.txt | cut -f1,2,3,4 | uniq > mart_export_GRCh37p13_short.txt
