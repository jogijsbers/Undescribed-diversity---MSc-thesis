#!/usr/bin/env bash
# Author: Johanna Gijsbers Alejandre
# Date created: April, 2021

### Runs admixture for RADseq data

# Argument variables
DATA_DIREC=$1
PATH_TO_VCF=$2
VCF_NAME=$3
DIREC=$4
NEW_DIREC=$5
START_K=$6
END_K=$7
THREADS=$8

## Make directory to save output
mkdir -p -v ${DATA_DIREC}/

## Convert 'vcf' file to plink
echo "Converting vcf file to plink format"
plink --vcf ${PATH_TO_VCF}/${VCF_NAME}_${DIREC}.vcf --const-fid --allow-extra-chr 0 --make-bed --out ${VCF_NAME}_${NEW_DIREC}

## Run admixture for N number of K
echo "---------------------------------- Starting admixture ----------------------------------"
for ((K=${START_K}; K<=${END_K}; K++)); \
do admixture --cv=10 "${VCF_NAME}_${NEW_DIREC}.bed" $K -j${THREADS} -C=100 | tee ./${DATA_DIREC}/log${K}.out;
done

## Move output Q and P to output
mv {*.P,*.Q} ${DATA_DIREC}/

## Move unlinked snps file and plink outputs to 'admixture' directory
mv ${VCF_NAME}_${NEW_DIREC}* ${DATA_DIREC}/

## save the likelihood results of each K to a single file

grep -h CV ./${DATA_DIREC}/log*.out > ./${DATA_DIREC}/${VCF_NAME}_Kerror_${NEW_DIREC}.txt

echo "---------------------------------- Finished admixture ----------------------------------"