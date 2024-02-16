#!/bin/bash

#######################################
# Arguments:
#   $1: plink reference
#   $2: plink2 reference
#   $3: PRSice path
# Returns:
#   None
#######################################


# current_directory='/data/users/itreccani/'
# local plink=$1
# local plink2=$2
# local path_PRSice=$3
# pathRcode=$current_directory/Rcode
# pathCommonFile= $current_directory/CommonFile
# pathADNI1=$current_directory/ADNI/ADNI1
# pathADNIGO2=$current_directory/ADNI/ADNIGO_2
# pathADNIomni=$current_directory/ADNI/ADNIomni
# pathADNIGO2_2nd=$current_directory/ADNI/ADNI_GO2_2nd
# pathADNI3=$current_directory/ADNI/ADNI3
# pathADNI3_Final=$current_directory/ADNI/ADNI3_Final
# paths=($pathADNI1 $pathADNIGO2 $pathADNIomni $pathADNIGO2_2nd $pathADNI3_Final $pathADNI3)


#####################################
#                                   #                               
#           FILEs REQUESTED         #
#                                   #
#####################################

cd $current_directory
mkdir -p CommonFile
cd $pathCommonFile

# 1000Genome 
wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst
wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst
wget https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam
# Decompressing the downloaded files using PLINK2
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
mv phase3_corrected.psam all_phase3.psam
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
plink2 --pfile all_phase3 vzs --snps-only just-acgt --max-alleles 2 --rm-dup exclude-mismatch --set-missing-var-ids '@_#_$1_$2' --make-bed --out all.phase3


# GRCh37
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz


# Download from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6836675/#S8
wget https://ctg.cncr.nl/documents/p1651/AD_sumstats_Jansenetal_2019sept.txt.gz
gunzip AD_sumstats_Jansenetal_2019sept.txt.gz
awk '!seen[$6]++' AD_sumstats_Jansenetal_2019sept.txt > unique_file.txt
Rscript $pathRcode/Base_file.R $pathCommonFile $current_directory/ADNI/PRS

wget https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip

# Some usefull links 
# https://rpubs.com/maffleur/452627 
# https://imputationserver.readthedocs.io/en/latest/prepare-your-data/#additional-tools
# https://shicheng-guo.github.io/archive
# https://rdrr.io/cran/plinkQC/man/evaluate_check_ancestry.html
# https://meyer-lab-cshl.github.io/plinkQC/
# https://choishingwan.github.io/PRSice/step_by_step/
# https://github.com/isglobal-brge/imputeInversion/blob/master/postimputation.sh
# https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
# https://cran.r-project.org/web/packages/plinkQC/vignettes/AncestryCheck.pdf


#---------------------------------------------------------------------------------------------------------------------


#####################################
#                                   #                          
#             QC FILTERS            #       
#                                   #                       
#####################################


for path in "${paths[@]}"; do

    if [ $path == $pathADNI1 ]; then
        inputfile='ADNI_cluster_01_forward_757LONI'
    else if [ $path == $pathADNIGO2 ]
        inputfile='ADNI_GO_2_Forward_Bin'
    else  if [ $path == $pathADNIomni ]
        inputfile='WGS_Omni25_BIN_wo_ConsentsIssues'
    else  if [ $path == $pathADNI3 ]
        inputfile='ADNI3_PLINK_FINAL_2nd'
    else  if [ $path == $pathADNI3_Final ]
        inputfile='ADNI3_PLINK_Final'
      else [ $path == $pathADNIGO2_2nd ]
        inputfile='ADNI_GO2_GWAS_2nd_orig_BIN'
    fi

    cd $path
    # Remove the variants whitout rsID
    #Rscript $pathRcode/Check_rs.R $path $inputfile.bim
    #plink --bfile $inputfile --extract rs.txt --make-bed --out ${inputfile}-rs

    # Update case and control
    Rscript $pathRcode/upload_pheno_AD.R $path $pathCommonFile ${inputfile}-rs

    cd $path
    # Check for consistency of sex assignments and inspect the subject-level missingness rates.
    plink --bfile ${inputfile}-rs --check-sex 
    plink --bfile ${inputfile}-rs --missing
    # Check if the status is OK and check if Proportion of missing SNPs is less of 10%
    Rscript $pathRcode/sex_missing.R $path $pathCommonFile
    # Remove the wrong samples
    plink --bfile ${inputfile}-rs --keep sex.valid  --make-bed --out ${inputfile}-sex 
    plink --bfile ${inputfile}-sex --keep miss.valid --make-bed --out ${inputfile}-qc 
    rm -f plink* ${inputfile}-sex* miss.valid sex.valid ${inputfile}-rs* rs* 
      
    cd $path
    # File needs to be consistent with the HRC site list. 
    plink --bfile ${inputfile}-qc --freq
    # Download in $pathCommonFile HRC-1000G   https://rpubs.com/maffleur/452627 
    perl $pathCommonFile/HRC-1000G/HRC-1000G-check-bim.pl -b ${inputfile}-qc.bim -f plink.frq -r $pathCommonFile/HRC-1000G/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
    # Script to check plink .bim files against HRC/1000G for strand, id names, positions, alleles, ref/alt assignment
    plink --bfile ${inputfile}-qc --exclude Exclude-${inputfile}-qc-HRC.txt --write-snplist --make-bed --out TEMP1
    plink --bfile TEMP1 --update-map Chromosome-${inputfile}-qc-HRC.txt --update-chr --make-bed --out TEMP2
    plink --bfile TEMP2 --update-map Position-${inputfile}-qc-HRC.txt --make-bed --out TEMP3
    plink --bfile TEMP3 --flip Strand-Flip-${inputfile}-qc-HRC.txt --make-bed --out TEMP4
    plink --bfile TEMP4 --a2-allele Force-Allele1-${inputfile}-qc-HRC.txt --make-bed --out TEMP6
    plink --bfile TEMP6 --recode vcf --out ${inputfile}-qc-updated
    rm -f TEMP1* TEMP2* TEMP3* TEMP4* plink.* *.txt Run-* *.log *.hh
  
    cd $path
    bcftools sort ${inputfile}-qc-updated.vcf |  bgzip -c > fixref.sorted1.vcf.gz
    bcftools index fixref.sorted1.vcf.gz
    bcftools +/data/users/itreccani/bcftools/plugins/fixref.so fixref.sorted1.vcf.gz -Oz -o fixref.sorted2.vcf.gz -- -d -f $pathCommonFile/checkVCF/hs37d5.fa -i $pathCommonFile/All_20180423.vcf.gz
    bcftools sort fixref.sorted2.vcf.gz |  bgzip -c > fixref..vcf.gz
    bcftools index fixref.sorted.vcf.gz
    bcftools norm --check-ref w -f $pathCommonFile/checkVCF/hs37d5.fa fixref.sorted.vcf.gz -Ou -o null        
    rm -f  test* fixref.sorted1* fixref.sorted2* 

for path in "${paths[@]}"; do  
    cd $path
    mkdir -p CHR
    cd $path/CHR
    # Splitting the VCF file by chromosome
    for i in {1..22}
      do
        bcftools view $path/fixref.sorted.vcf.gz --regions ${i} -o $path/CHR/VCF.${i}.vcf.gz -Oz
    done

done


#---------------------------------------------------------------------------------------------------------------------


######     Michigan Imputation Server      ######   
# https://imputationserver.sph.umich.edu/index.html#!run/minimac4
## STEPS
# 1 Select: Run -> Genotype Imputation (Minimac4)
# 2 Reference Panel: 1000G Phase 1 v5
# 3 Input files: VCF.gz
# 4 Array Build: GRCh37/hg19
# 5 rsq Filter: NULL
# 6 Phasing: Eagle2
# 7 Population: EUR
# 8 Mode: Quality Control & Imputation


#---------------------------------------------------------------------------------------------------------------------


#####################################
#                                   #
#          POST IMPUTATION          #
#                                   #
#####################################

for path in "${paths[@]}"; do

    # Create the directory e move the file in the specific directory
    cd $path/CHR
    for i in {1..22}
    do
        mkdir -p chr_$i
        mv $path/CHR/chr_$i.zip $path/CHR/chr_$i/chr_$i.zip
    done

    # Unzip the files, the carefolly that it necessary of the password
    cd $path/CHR
    for i in {1..22}
    do
        cd chr_$i
        unzip chr_$i.zip 
        cd ..
    done

    cd $path/CHR
    for i in {1..22}
    do
        cd chr_$i
        input_basename1=chr$i.dose.vcf.gz
        input_basename='chr.vcf.gz'
        # Extract the INFO metric from VCF file for ALL variants
        vcftools --gzvcf ${input_basename1} --get-INFO INFO --out ${input_basename}
        # Write a list of variants with INFO metric < 0.3
        awk '$NF < 0.3' ${input_basename}.INFO | cut -f 1-2 > ${input_basename}_ExcludePositions.txt
        # Call VCFtools to exclude those variants by position
        vcftools --gzvcf ${input_basename1} --exclude-positions ${input_basename}_ExcludePositions.txt --recode --out ${input_basename}_clean
        # For non-header lines, it checks the length of the third field. If it exceeds 16000 characters, it truncates the field to 16000 characters.
        awk 'BEGIN { OFS="\t" }
        {
            if ($0 ~ /^#/) {
                print $0
            } else {
                if (length($3) > 16000) {
                    $3 = substr($3, 1, 16000)
                }
                print $0
            }
        }' ${input_basename}_clean.recode.vcf > ${input_basename}_clean2.recode.vcf

      # Sord and index the .vcf files
        bcftools sort ${input_basename}_clean2.recode.vcf |  bgzip -c > ${input_basename}_clean2.recode.vcf.gz
        bcftools index ${input_basename}_clean2.recode.vcf.gz
        cd ..
    done

    # Create a lift of file that will be merged
    cd $path/CHR
    rm -f file.list
    for i in {1..22}
    do 
        echo $path/CHR/chr_$i/${input_basename}_clean2.recode.vcf.gz >> file.list
    done
    # Merge the files 
    bcftools concat -f file.list --allow-overlaps --remove-duplicates  -Ov -o filtered_BC_merged.vcf
    plink2 --vcf filtered_BC_merged.vcf --double-id  --make-bed --out ${input_basename}
    # Performe some QC filtering
    plink2 --bfile ${input_basename} --snps-only just-acgt --max-alleles 2 --maf 0.05 --geno 0.1 --hwe 5e-7 --make-bed --out prova.merged-QC2
    # Update sex, RID and FID names, phenotype
    Rscript $pathRcode/post_imputation_fam.R $path/CHR $pathCommonFile prova.merged-QC2 $path
    # Remove samples not conforms to the study
    plink2 --bfile prova.merged-QC2 --remove-fam hispanic.txt --make-bed --out Clean1list_biallelic

done

# Merge the datasets
cd $current_directory/ADNI
plink --bfile $pathADNIomni/CHR/Clean1list_biallelic --bmerge $pathADNI1/CHR/Clean1list_biallelic  --make-bed --out $current_directory/ADNI/file_merge1
plink --bfile   $pathADNI3/CHR/Clean1list_biallelic --bmerge $pathADNIGO2/CHR/Clean1list_biallelic   --make-bed --out $current_directory/ADNI/file_merge2
plink --bfile  $pathADNIGO2_2nd/CHR/Clean1list_biallelic --bmerge $pathADNI3_Final/CHR/Clean1list_biallelic --make-bed --out $current_directory/ADNI/file_merge3
plink --bfile  $current_directory/ADNI/file_merge1 --bmerge $current_directory/ADNI/file_merge2 --make-bed --out $current_directory/ADNI/file_merge4
plink --bfile  $current_directory/ADNI/file_merge4 --bmerge $current_directory/ADNI/file_merge3 --make-bed --out $current_directory/ADNI/file_merge5
# Recode rsid
plink2 --bfile $current_directory/ADNI/file_merge5 --recover-var-ids $pathCommonFile/All_20180423.vcf.gz force partial --make-bed --out $current_directory/ADNI/file_merge6
# Perform some QC filers
plink --bfile $current_directory/ADNI/file_merge6 --maf 0.05 --geno 0.1 --make-bed --out  $current_directory/ADNI/postrecover


#---------------------------------------------------------------------------------------------------------------------

 
#####################################
#                                   #                                   
#       POPULATION ANCESTRY         # 
#                                   #                                 
#####################################


cd $current_directory/ADNI
mkdir -p PCA
cd $current_directory/ADNI/PCA
name='postrecover'
refname='all.phase3'

plink --bfile $current_directory/ADNI/postrecover --make-bed --out Target
# Filter reference and study data for non A-T or G-C SNPs
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' Target.bim > Target.ac_gt_snps
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}'  $pathCommonFile/all.phase3.bim > all.phase3.ac_gt_snps
# Create new data sets by excluding the markers identified in the previous step
plink --bfile  $pathCommonFile/all.phase3 --exclude all.phase3.ac_gt_snps --make-bed --out all.phase3.no_ac_gt_snps  
plink --bfile Target --exclude Target.ac_gt_snps --make-bed  --out Target.no_ac_gt_snps  
# Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' Target.no_ac_gt_snps.bim Target.no_ac_gt_snps.bim | sed -n '/^[XY]/!p' > all.phase3.toUpdateChr
plink --bfile all.phase3.no_ac_gt_snps --update-chr all.phase3.toUpdateChr 1 2  --make-bed --out all.phase3.updateChr
# Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} ($2 in a && a[$2] != $4) {print a[$2],$2}' Target.no_ac_gt_snps.bim all.phase3.no_ac_gt_snps.bim > all.phase3.toUpdatePos
# Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' Target.no_ac_gt_snps.bim all.phase3.no_ac_gt_snps.bim > all.phase3.toFlip
# Upate positions and flip alleles
plink --bfile all.phase3.updateChr --update-map all.phase3.toUpdatePos 1 2 --flip all.phase3.toFlip --make-bed  --out all.phase3.flipped     
# Remove mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' Target.no_ac_gt_snps.bim all.phase3.autosome.flipped.bim > all.phase3.autosome.mismatch
plink --bfile all.phase3.flipped --exclude all.phase3.autosome.mismatch --make-bed  --out all.phase3.clean
plink --bfile Target.no_ac_gt_snps --bmerge all.phase3.clean  --make-bed --out Target.merge.all.phase3
# PCA
plink --bfile Target.merge.all.phase3 --pca 5 --out Target.all.phase3

# Check the ancestry
Rscript $pathRcode/evaluate_check_ancestry.R $current_directory/ADNI/PCA Target

# Check the ancestry second modality
Rscript $pathRcode/PCA-plot.R $current_directory/ADNI/PCA $pathCommonFile PCA

# Remove the wrong samples
plink --bfile $current_directory/ADNI/postrecover --remove remove.txt --make-bed --out ADNI.Final 
plink --bfile ADNI.Final --pca --out PCA

# --indep-pairwise takes the same first two parameters as --indep. Its third parameter is a pairwise r2 threshold:
# at each step, pairs of variants in the current window with squared correlation greater than the threshold are noted,
# and variants are greedily pruned from the window until no such pairs remain.
cd $current_directory/ADNI/PCA
plink --bfile ADNI.Final --indep-pairwise 1000 50 0.1 --out ADNI.Final_snps
plink --bfile ADNI.Final --extract ADNI.Final_snps.prune.in --make-bed --out Target.Final.pruned2
plink --bfile Target.Final.pruned2 --rel-cutoff 0.1 --out Target.Final.pruned2 
plink --bfile Target.Final.pruned2 --keep Target.Final.pruned2.rel.id --make-bed --out $current_directory/ADNI/PRS/Target.Final.pruned
    

#---------------------------------------------------------------------------------------------------------------------

 
#####################################
#                                   #
#               PRS                 #
#                                   #
#####################################


cd $current_directory/ADNI
mkdir -p PRS
#cd $current_directory/PRS
cd $current_directory/ADNI/PRS

# Create Pheno and Covariate files
Rscript $pathRcode/pheno_covariate.R $current_directory/ADNI/PRS $current_directory/ADNI/PCA Target.Final.pruned

# Calculate PRS
cd $current_directory/ADNI/PRS
Rscript $path_PRSice/PRSice.R \
  --prsice $path_PRSice/PRSice_linux \
  --base $pathCommonFile/snp_ref_prscs_PRSice.txt \
  --target Target.Final.pruned \
  --clump-kb 1000 \
  --clump-r2 0.1 \
  --pheno recodedpheno.txt \
  --cov covariate.txt \
  --all-score \
  --out PRS1

# Desity plot
Rscript $pathRcode/PRSice.plot.R $current_directory/ADNI/PRS

#cd $current_directory/PRS-cs
#python $current_directory/PRS-cs/PRScs/PRScs.py --ref_dir=$current_directory/PRS-cs/ldblk_1kg_eur --bim_prefix=$current_directory/PRS/Target.Final.pruned --sst_file=$pathCommonFile/snp_ref_prscs.txt --n_gwas=1274 --out_dir=$current_directory/PRS-cs
#plink --bfile $current_directory/ADNI/PRS/Target.Final.pruned --score merge_result_prscs.txt --make-bed --out prscs
    

#---------------------------------------------------------------------------------------------------------------------


#####################################
#                                   #
#             PLINK PRS             #
#                                   #
#####################################

cd $current_directory/ADNI/PRS
plink --bfile $current_directory/ADNI/PRS/Target.Final.pruned --clump-p1 1 --clump-r2 0.1 --clump-kb 1000 --clump $pathCommonFile/snp_ref_prscs_PRSice.txt --clump-snp-field SNP --clump-field P --out clumping
      
awk 'NR!=1{print $3}' clumping.clumped >  clumping.valid.snp

# range labels in the first column, p-value lower bounds in the second column, and upper bounds in the third column
awk '{print $1,$6}'  $pathCommonFile/snp_ref_prscs_PRSice.txt > SNP.pvalue
echo "0.001 0 0.001" > range_list 
echo "0.01 0 0.01" >> range_list
echo "0.02 0 0.02" >> range_list
echo "0.03 0 0.03" >> range_list
echo "0.04 0 0.04" >> range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

# Command to calculate PRS 
plink --bfile $current_directory/ADNI/PRS/Target.Final.pruned \
      --score $pathCommonFile/snp_ref_prscs_PRSice.txt 1 2 4 header \
      --q-score-range range_list SNP.pvalue \
      --extract clumping.valid.snp \
      --out prs.plink

Rscript $pathRcode/clumping.fit.R
Rscript $pathRcode/plink.barplot.R 
Rscript $pathRcode/plink.plot.R $current_directory/ADNI/PRS
