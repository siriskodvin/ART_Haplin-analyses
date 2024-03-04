#!/bin/sh

#### -------------------------------- Script header -------------------------------- ####
# Date:         04 March 2024                                                           #
# Author:       Siri Naerland Skodvin                                                   #
# Filename:     FRchr.sh                                                                #
# Description:  Sequential runs of PLINK, to recode genetic data, and R, to analyze     #
#               data using Haplin and plot results.                                     #
#### --------------------------------------------------------------------------------####

# Input chromosome - CHANGE HERE
chrchr=21       # chromosome number in two characters
chrnum=21       # chromosome number with minimum characters
startmbp=0      # start mbp of chromosome
stopmbp=50      # stop mbp of chromosome

# Input chunks - CAN USE DEFAULT
bymbp=1         # number of mbp to include in each gen recode file

# Input pruning - CAN USE DEFAULT
window=50
step=5
r2=0.95

# Preprocess input
cpath=/PATH_TO_CODE_AND_RESULTS/
dpath=/PATH_TO_TEMPORARY_STORAGE_AND_MAIN_R_SCRIPTS/
seq=`expr $stopmbp - $startmbp`
mod=`expr $seq % $bymbp`
if [ $mod -ne 0 ]
then
	stopmbp=`expr $stopmbp - $mod + $bymbp`
	seq=`expr $stopmbp - $startmbp`
fi
chunkno=`expr $seq / $bymbp`

# Initiate log file
now=`date`
cat <<EOF > $dpath/log_chr$chrchr.txt
INPUT
Chromosome: $chrnum
Start-mbp: $startmbp
Stop-mbp: $stopmbp
By-mbp: $bymbp
Sequences: $chunkno
Prune window-step-r2: $window - $step - $r2
ANALYSIS STARTED
$now
ANALYZED SEQUENCES
EOF

# Main analysis - Plink and R in sequence
i=1
while [ $i -le $chunkno ]
do  
    for file in `ls $dpath`
	do
		filepre=`expr "$file" | cut -c1-3`
		if [ $filepre = "gen" ]
		then
			rm $dpath/$file
		fi
		if [ $filepre = "ld." ]
		then
			rm $dpath/$file
		fi
	done
    
    echo ---------------------------------------
    echo ---- PLINK SEQUENCE $i of $chunkno ----
    echo ---------------------------------------
    
    frommb=`expr $startmbp + \( $bymbp \* \( $i - 1 \) \)`
    tomb=`expr $frommb + $bymbp`
    
    /opt/plink --fam GENETICS.fam --bim GENETICS.bim --bed GENETICS.bed --maf 0.01 --chr $chrnum --from-mb $frommb --to-mb $tomb --indep-pairwise $window $step $r2 --out $dpath/ld
    
    /opt/plink --fam GENETICS.fam --bim GENETICS.bim --bed GENETICS.bed --maf 0.01 --chr $chrnum --extract $dpath/ld.prune.in --recode --out $dpath/gen
    
    echo ------------------------------------------
    echo ---- R HAPLIN SEQUENCE $i of $chunkno ----
    echo ------------------------------------------
    
    k=0
	for file in `ls $dpath`
	do
		ext=`echo "$file" | cut -d'.' -f2`
		if [ $ext = "ped" ]
		then
			k=`expr $k + 1`
		fi
	done
	
	if [ $k -eq 0 ]
	then
		echo ---- => NO VARIANTS ----
        newline1="-----------------------------------------------"
        newline2="seq $i: ${frommb} to ${tomb} mbp => no variants"
	else
		Rscript --vanilla $dpath/FRhap.R
        newline1="--------------------------------"
        newline2="seq $i: ${frommb} to ${tomb} mbp"
	fi
    
    # Update log file
    echo $newline1 >>$dpath/log_chr$chrchr.txt
    echo $newline2 >>$dpath/log_chr$chrchr.txt
    echo $newline1 >>$dpath/log_chr$chrchr.txt
    for file in `ls $cpath/Results/Chr$chrchr`
    do
        ext=`echo "$file" | cut -d'.' -f2`
        if [ $ext = "txt" ]
        then
            cat $cpath/Results/Chr$chrchr/$file >> $dpath/log_chr$chrchr.txt
            rm $cpath/Results/Chr$chrchr/$file
        fi
    done
    
    i=`expr $i + 1`
done
Rscript --vanilla $dpath/FRmpl.R

# Complete log file
lastline="ANALYSIS FINISHED"
now=`date`
vartotline="NUMBER OF VARIANTS"
vartot=`cat $cpath/Results/Chr$chrchr/numvar$chrchr.txt`
echo $lastline >>$dpath/log_chr$chrchr.txt 
echo $now >>$dpath/log_chr$chrchr.txt
echo $vartotline >>$dpath/log_chr$chrchr.txt
echo $vartot >>$dpath/log_chr$chrchr.txt
cp $dpath/log_chr$chrchr.txt $cpath/Results/Chr$chrchr
rm $cpath/Results/Chr$chrchr/numvar$chrchr.txt

echo ---- ANALYSIS FINISHED ----









