
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
	python structureHarvester.py --dir=structure/ --out=structure_results/ --clumpp 

	for k in {2..4}
	do
		for class in pop ind 
		do
			CLUMPP_Linux64.1.1.2/CLUMPP clumppParam.${class}.k${k}
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
		vcftools --gzvcf GBS.mis60.maf01.recode.vcf.gz --out HiC_hyb.fst.${pop1}_${pop2} --weir-fst-pop ${pop1}.inds --weir-fst-pop ${pop2}.inds 2>fst.${pop1}_${pop2}.txt
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

bgcdist/bgc -a XY.bgc.txt -b XYY.bgc.txt -h admx.bgc.txt -O 0 -x 50000 -n 10000 -p 1 -q 1 -N 1 -t 5 -F GBS.mis60.maf01.noY.bgc40k

for stat in hi alpha beta gamma-quantile zeta-quantile
do
	bgcdist/estpost -i GBS.mis60.maf01.noY.bgc40k.hdf5 -p ${stat} -o GBS.mis60.maf01.noY.bgc40k.${stat}.out -s 0 -c 0.95 #-w 0
done
	

##Bayenv 
#https://gitlab.com/YDorant/Toolbox/blob/master/reshaper_baypass.py

python3 reshaper_baypass.py GBS.mis60.maf01.noY.bayenv.XY.recode.vcf HiC_hyb.XY.pop.sort GBS.mis60.maf01.noY.XY.geno

baypass_2.2/sources/g_baypass -gfile GBS.mis60.maf01.noY.XY.geno -nthreads 10 -efile bayenv_XY_elev.txt -scalecov -outprefix GBS.mis60.maf01.noY.XY.baypass.STD 

####STRUCTURE GWAS####

gemma-0.98.1-linux-static -bfile GBS.mis60.maf01.noY.XY.plink -lm 2 -o GBS.mis60.maf01.noY.XY.k3.gemma


##LD within XX/XY

plink_linux_x86_64_20200616/plink --file GBS.mis60.maf01.XY.plink --r2 inter-chr  --out GBS.XY.plink.ld 
