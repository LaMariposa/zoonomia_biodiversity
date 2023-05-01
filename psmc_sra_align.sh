#!/bin/bash

#Table ST1 from first zoonmia paper with SRA info
#add mutation rate and generation length for each species 
infile=/scratch1/zoonomia_biodiversity/FinalSupplementalTables_ST1.tsv

threads=6

genomesdir=/scratch1/zoonomia_biodiversity/broad_genomes/
samplesdir=/scratch1/zoonomia_biodiversity/psmc/

sratoolkit=/home/megan/sratoolkit.2.10.9-ubuntu64/bin/
psmc=/home/megan/psmc


#loop over samples from csv (line 5 is first entry)
for ((a=5; a<246; a++))  #*******************************************
do
	# read in entry metadata
	line=$(sed -n "${a}p" "$infile")
	#echo $line

	# get metadata
	IFS=$'\t' read -ra deets <<< "$line"
	echo -e "\n\n\nprocessing sample number: ${deets[21]}, species: ${deets[1]}"
	temp="${deets[10]%\"}"
        refID="${temp#\"}"
	temp="${deets[6]%\"}"
	sampID="${temp#\"}"
	temp="${deets[1]%\"}"
        temp="${temp#\"}"
	sp=${temp// /_}
	echo $sp
	echo $sampID
	echo $refID

	cd ${samplesdir}

	# get SRA reads
	${sratoolkit}/prefetch --output-directory ${sp} --max-size 200GB ${sampID}
	sraID=`ls ${sp}`
	echo ${sraID}
	cd ${sp}/${sraID}
	${sratoolkit}/fasterq-dump --threads ${threads} ${sraID}.sra
	bgzip --threads ${threads} ${sraID}_1.fastq
	bgzip --threads ${threads} ${sraID}_2.fastq

	# reference genome
       	ref=${genomesdir}/${sp}.fa.gz
	bwa index ${ref}
	samtools faidx ${ref}
	gzip ${genomesdir}/${sp}.fa.fai

	# align reads to reference
	bwa mem -t ${threads} -PM -R "@RG\tID:${sraID}\tLB:${sraID}\tSM:${sampID}\tPL:ILLUMINA" \
		${ref} \
         	${sraID}_1.fastq.gz ${sraID}_2.fastq.gz |
		samtools view -@ ${threads} -b |
		samtools sort -@ ${threads} -T tmpsort -o ${sraID}.bam
	samtools index -@ ${threads} ${sraID}.bam

        # extract regions
        zless ${genomesdir}/${sp}.fa.fai.gz | awk '{if ($2>50000) print $1, "0", $2}' > ${sp}.50k.bed
	samtools view -@ ${threads} -bh -L ${sp}.50k.bed -o ${sraID}.50k.bam ${sraID}.bam
	samtools index -@ ${threads} ${sraID}.50k.bam

	# stats
	samtools flagstat -@ ${threads} ${sraID}.50k.bam > ${sraID}.flagstat
	samtools depth ${sraID}.50k.bam | 
		awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${sraID}.50k.coverage.txt

	# prep for PSMC
	cov=`awk '{printf("%d"), $3}' ${sraID}.50k.coverage.txt`
	mincov=`expr $cov / 3`
	maxcov=`expr $cov \* 2`
	samtools mpileup -A -C50 -uf ${ref} ${sraID}.50k.bam | bcftools call -c - | 
		vcfutils.pl vcf2fq -d ${mincov} -D ${maxcov} | gzip > ${sraID}.con.fq.gz
	${psmc}/utils/fq2psmcfa -q20 ${sraID}.con.fq.gz > ${sraID}.psmcfa
	
	# run PSMC
	${psmc}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sraID}.psmc ${sraID}.psmcfa
  	${psmc}/utils/psmc_plot.pl ${sraID} ${sraID}.psmc

done
