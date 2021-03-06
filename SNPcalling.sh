
#trim and align
	while read ind
	do
	TrimGalore-0.4.3/trim_galore -o GBSpops/trimmed/ GBSpops/40014_${ind}_il.fastq.gz
	Cibiv-NextGenMap-33e92fb/bin/ngm-0.5.5/ngm -t 10 -r hastate_28Sep2018_nbKXS_50KUP.masked.fasta -b -p -1 GBShybrid/trimmed/R1hyb_Rh${ind}_R1_val_1.fq.gz -2 GBShybrid/trimmed/R2hyb_Rh${ind}_R2_val_2.fq.gz -o alignments/HiC/ngm/${ind}.ngm.bam 
	done <GBSpops/inds.list

#Picard; sort + AddOrReplace
	
	while read ind
	do
	java -jar picard/build/libs/picard.jar SortSam SORT_ORDER=coordinate I=alignments/HiC/ngm/${ind}.ngm.bam O=alignments/HiC/ngm/${ind}.sort.temp.bam VALIDATION_STRINGENCY=LENIENT
	java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=HiC/ngm/${ind}.sort.temp.bam O=HiC/ngm/${ind}.add.bam RGID=${ind} RGLB=${pop} RGPL=GBS RGPU=unit1 RGSM=${ind}
	done <GBSpops/inds.list
	
#call variant sites
bcftools/bcftools mpileup -a AD,DP -Ou -f hastate_28Sep2018_nbKXS_50KUP.masked.fasta alignments/HiC/ngm/*.add.bam | bcftools/bcftools call -mO z -f gq | bcftools/bcftools view --threads 10 -O z -o alignments/HiC/ngm/GBS.vcf.gz

#### filtering#### alignments/HiC/ngm/
#all filtering
	vcftools --gzvcf alignments/HiC/ngm/GBS.vcf.gz --remove-indels --max-alleles 2 --minQ 10  --recode --out GBS.filt1 --remove bad.inds
	vcftools --vcf GBS.filt1.recode.vcf --min-meanDP 5 --minGQ 5 --recode --out GBS.filt2
	vcftools --vcf GBS.filt2.recode.vcf --recode --max-missing 0.6 --out GBS.mis60
	vcftools --vcf GBS.mis60.recode.vcf --recode --maf 0.01 --out alignments/HiC/ngm/GBS.mis60.maf01 

awk 'NR==FNR{a[$1]=1;next}{if(!a[FNR])print$0}' bad_SNPs_vcf.txt GBS.mis60.maf01.recode.vcf  > GBS.mis60.maf01.noY.recode.vcf

#format conversion
	vcftools --vcf GBS.mis60.maf01.noY.recode.vcf --012 --out GBS.mis60.maf01.noY
	vcftools --vcf GBS.filt2.recode.vcf --plink --maf 0.01 --out alignments/HiC/ngm/GBS.mis60.maf01

##convert positions from contigs to chromosomes##
	awk 'split($1,a,"_") split(a[2],b,"-") {print b[1]","$2}' GBS.mis60.maf01.012.pos >GBS.mis60.maf01.012.pos_trim
	python3 /ohta/felix.beaudry/scripts/HiC2Chrom/SNP_list_Interactive_Chromonomer_LG_position_switcher_FB.py GBS.mis60.maf01.012.pos_trim GBS.mis60.maf01.012.chrom.pos

####pop structure####
	##STRUCTURE##
	#Rscript --vanilla StructureMaker.R #012 matrix #inds|pops #outfile name
	Rscript --vanilla StructureMaker.R GBS.mis60.maf01.noY.012 GBS.pop.num GBS.mis60.maf01.noY.popd.012

	for k in {1..6}
	do
	for iter in {1..10}
	do
	structure -m mainparams -e extraparams -o GBS.mis60.maf01.noY.strc.k${k}.it${iter}.strc -K ${k}
	done
	done

	mv GBS.mis60.maf01.noY.strc.k*.it*.strc_f structure/ 
	python /ohta/felix.beaudry/scripts/structureHarvester-master/structureHarvester.py --dir=structure/ --out=structure_results/ --clumpp 

	for k in {2..4}
	do
	for class in pop ind 
	do
	/ohta/felix.beaudry/scripts/CLUMPP_Linux64.1.1.2/CLUMPP clumppParam.${class}.k${k}
	done
	done

##EEMS##
	plink-1.07-x86_64/plink --noweb --file HZEEMS --make-bed --out HZEEMS
	eems/bed2diffs/src/bed2diffs_v2 --bfile ./HZEEMS
	
	for chain in {1..6}
	do
	tr chain1 chain${chain} <params-chain.ini >params-chain${chain}.ini 
	eems/runeems_snps/src/runeems_snps --params HiC/ngm/params-chain${chain}.ini --numThinIter 9999 --numBurnIter 1000000 --numMCMCIter 2000000 --mcmcpath ./GBShyb-EEMS-nDemes800-chain${chain} --nDemes 800 --nIndiv 150 --nSites 13640
	done

##heterozygosity
plink-1.07-x86_64/plink --noweb --file HiC/ngm/GBS.mis60.maf01 --het


####FST####
for pop1 in WEIR UTIC TXROS TXMTP TXOAK TXLIV TXATH SYLV STRO SIBL SCPRO SCMAR SCBRA RUST OKRAT OKBAC OGAC NWHO NCROS NCELI NCBAT MTHO LOUI LEES LABEN KELL JONE GLEN GASTA GAGLA GABEL FOUK FLMAR FLJAS ENON ELDO CROS CALV CALH BUCK BRAN BIEN ARCA ALBRU ALBRE
do
for pop2 in UTIC WEIR TXROS TXMTP TXOAK TXLIV TXATH SYLV STRO SIBL SCPRO SCMAR SCBRA RUST OKRAT OKBAC OGAC NWHO NCROS NCELI NCBAT MTHO LOUI LEES LABEN KELL JONE GLEN GASTA GAGLA GABEL FOUK FLMAR FLJAS ENON ELDO CROS CALV CALH BUCK BRAN BIEN ARCA ALBRU ALBRE
do
vcftools --gzvcf /ohta/felix.beaudry/alignments/HiC/ngm/GBS.mis60.maf01.recode.vcf.gz --out HiC_hyb.fst.${pop1}_${pop2} --weir-fst-pop ${pop1}.inds --weir-fst-pop ${pop2}.inds 2>fst.${pop1}_${pop2}.txt
#awk -v a=${pop1} -v b=${pop2} '$4 ~ "mean" {print $7" "a" "b}' fst.${pop1}_${pop2}.txt >>fst.mean.txt
awk -v a=${pop1} -v b=${pop2} '$4 ~ "weighted" {print $7" "a" "b}' fst.${pop1}_${pop2}.txt >>fst.weighted.txt
done
done

##BGC #https://sites.google.com/site/bgcsoftware/
	for pop in XY XYY admx
	do
	vcftools --vcf GBS.mis60.maf01.noY.recode.vcf --out GBS.mis60.maf01.noY.bgc.${pop} --keep GBS.${pop}.inds --recode #--positions GBS.mis60.maf01.noY.012.pos 
	done
	
	perl vcf2bgc.pl -i GBS.mis60.maf01.noY.bgc.XY.recode.vcf -a F >XY.bgc.txt
	perl vcf2bgc.pl -i GBS.mis60.maf01.noY.bgc.XYY.recode.vcf -a F >XYY.bgc.txt
	perl vcf2bgc.pl -i GBS.mis60.maf01.noY.bgc.admx.recode.vcf -a T >admx.bgc.txt

	/ohta/apps/bgcdist/bgc -a XY.bgc.txt -b XYY.bgc.txt -h admx.bgc.txt -O 0 -x 50000 -n 10000 -p 1 -q 1 -N 1 -t 5 -F GBS.mis60.maf01.noY.bgc40k

	for stat in hi alpha beta gamma-quantile zeta-quantile
	do
	/ohta/apps/bgcdist/estpost -i GBS.mis60.maf01.noY.bgc40k.hdf5 -p ${stat} -o GBS.mis60.maf01.noY.bgc40k.${stat}.out -s 0 -c 0.95 #-w 0
	done
	
for pop in XY XYY admx
	do
	vcftools --vcf GBS.mis60.maf01.recode.vcf  --out GBS.mis60.maf01.bgc.${pop} --keep GBS.${pop}.inds --recode #--positions GBS.mis60.maf01.noY.012.pos 
	done
	
	perl vcf2bgc.pl -i GBS.mis60.maf01.bgc.XY.recode.vcf -a F >XY.nofilt.bgc.txt
	perl vcf2bgc.pl -i GBS.mis60.maf01.bgc.XYY.recode.vcf -a F >XYY.nofilt.bgc.txt
	perl vcf2bgc.pl -i GBS.mis60.maf01.bgc.admx.recode.vcf -a T >admx.nofilt.bgc.txt

	/ohta/apps/bgcdist/bgc -a XY.nofilt.bgc.txt -b XYY.nofilt.bgc.txt -h admx.nofilt.bgc.txt -O 0 -x 50000 -n 10000 -p 1 -q 1 -N 1 -t 5 -F GBS.mis60.maf01.bgc40k

	for stat in hi alpha beta gamma-quantile zeta-quantile
	do
	/ohta/apps/bgcdist/estpost -i GBS.mis60.maf01.bgc40k.hdf5 -p ${stat} -o GBS.mis60.maf01.bgc40k.${stat}.out -s 0 -c 0.95 #-w 0
	done

	vcftools --vcf GBS.mis60.maf01.recode.vcf  --out GBS.mis60.maf01.bgc --012
	awk 'split($1,a,"_") split(a[2],b,"-") {print b[1]","$2}' GBS.mis60.maf01.bgc.012.pos >GBS.mis60.maf01.bgc.012.pos_trim
	python3 /ohta/felix.beaudry/scripts/HiC2Chrom/SNP_list_Interactive_Chromonomer_LG_position_switcher_FB.py GBS.mis60.maf01.bgc.012.pos_trim GBS.mis60.maf01.bgc.012.chrom.pos


##Bayenv 
#https://gitlab.com/YDorant/Toolbox/blob/master/reshaper_baypass.py
vcftools --vcf GBS.mis60.maf01.noY.recode.vcf --out GBS.mis60.maf01.noY.bayenv.XY --keep GBS.XY.bayenv.inds --recode #--positions GBS.mis60.maf01.noY.012.pos 

python3 /ohta/felix.beaudry/scripts/reshaper_baypass.py GBS.mis60.maf01.noY.bayenv.XY.recode.vcf HiC_hyb.XY.pop.sort GBS.mis60.maf01.noY.XY.geno


#for rep in {1..3}
#do
#baypass_2.2/sources/g_baypass -gfile GBS.mis60.maf01.noY.XY.geno -nthreads 10 -outprefix GBS.mis60.maf01.noY.XY.baypass.rep${rep}.AUX
#done

/ohta/apps/baypass_2.2/sources/g_baypass -gfile GBS.mis60.maf01.noY.XY.geno -nthreads 10 -efile bayenv_XY_elev.txt -scalecov -outprefix GBS.mis60.maf01.noY.XY.baypass.STD 

####STRUCTURE GWAS####

vcftools --vcf GBS.mis60.maf01.noY.recode.vcf --plink --keep GBS.XY.inds --out GBS.mis60.maf01.noY.XY.plink


~/scripts/plink-1.07-x86_64/plink --file GBS.mis60.maf01.noY.XY.plink --make-bed --out GBS.mis60.maf01.noY.XY.plink --noweb 

cp strc_k3.fam GBS.mis60.maf01.noY.XY.plink.fam
##manual step of assigning phenotypes in .fam file
/ohta/haoran.xue/programs/gemma-0.98.1-linux-static -bfile GBS.mis60.maf01.noY.XY.plink -lm 2 -o GBS.mis60.maf01.noY.XY.k3.gemma


##LD within XX/XY
#python3 /ohta1/joanna.rifkin/Chromonomer_analyses/Rumex/VCF_converter/VCF_Interactive_Chromonomer_LG_position_switcher_recognize_dash.py
/ohta/joanna.rifkin/Chromonomer_analyses/Rumex/TX_final/Final_chromonomer_run_sorted_positions_fixed/backup_copy_TX_final_CHRR_genome.agp
/ohta1/felix.beaudry/alignments/HiC/ngm/GBS.mis60.maf01.XY.recode.vcf
/ohta1/felix.beaudry/alignments/HiC/ngm/converted_GBS.mis60.maf01.XY.recode.vcf

#/ohta1/felix.beaudry/scripts/plink-1.07-x86_64/plink --noweb --bfile GBS.A2.XY  --make-bed --out GBS.A2.XY
#/ohta1/felix.beaudry/scripts/plink_linux_x86_64_20200616/plink --bfile GBS.A2.XY  --recode --out GBS.A2.XY.plink
awk '$1 ~ "#" {print}' converted_GBS.mis60.maf01.XY.recode.vcf | awk '$1 !~ "contig" {print}' >converted_GBS.mis60.maf01.XY.reorder.vcf

awk '$1 !~ "#" {print}' converted_GBS.mis60.maf01.XY.recode.vcf | sort -k1,1 -k2,2n >>converted_GBS.mis60.maf01.XY.reorder.vcf
vcftools --vcf converted_GBS.mis60.maf01.XY.reorder.vcf --plink --out GBS.XY.v2  --keep GBS.XY.a.inds 

awk 'split($2,a,":")  {print a[1]" "$2" "$3" "$4}' GBS.XY.v2.map >GBS.XY.v2.map.adj

#/ohta1/felix.beaudry/scripts/plink_linux_x86_64_20200616/plink --file GBS.XY --r2 --ld-window-kb 100000 --ld-window-r2 0 --out GBS.XY.plink #   # inter-chr 
/ohta1/felix.beaudry/scripts/plink_linux_x86_64_20200616/plink --file GBS.XY.v2 --r2 inter-chr --out GBS.XY.v2 --ld-window-r2 0  --allow-extra-chr #--ld-window 1000  --ld-window-r2 0  #   #
awk '$1 == $4 {print}' GBS.XY.v2.ld >GBS.XY.v2.ld.within


#/ohta1/felix.beaudry/scripts/plink_linux_x86_64_20200616/plink --file GBS.mis60.maf01.XY.plink --r2 inter-chr  --out GBS.XY.plink.ld 

 awk 'split($3,a,"_") split(a[2],b,"-HR") split($6,c,"_") split(c[2],d,"-HR") {print b[1]" "$2" "d[1]" "$5" "$7}' GBS.XY.plink.ld.ld >GBS.XY.plink.ld

 python3 VCF_Interactive_Chromonomer_LG_position_switcher.py /ohta1/felix.beaudry/alignments/HiC/ngm/GBS.mis60.maf01.XY.recode.vcf 

cd /ohta/felix.beaudry/scripts/HiC2Chrom/
python3 VCF_Interactive_Chromonomer_LG_position_switcher.py /ohta1/felix.beaudry/alignments/HiC/ngm/GBS.mis60.maf01.XY.recode.vcf /ohta1/felix.beaudry/alignments/HiC/ngm/GBS.mis60.maf01.XY.recode.vc.chromd