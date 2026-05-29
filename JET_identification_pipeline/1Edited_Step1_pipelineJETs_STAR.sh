#!/bin/bash

#### created by Alexandre Houy and Christel Goudot
#### modified and adapted by Muhammad Faramir

#-------------------------------------------------------
##-- Configuration
#-------------------------------------------------------

##### Software paths
samtoolsBinDir=/home/faramir/repos/neondisco-jet/.pixi/envs/default/bin
starBinDir=/home/faramir/repos/neondisco-jet/.pixi/envs/default/bin

##### Run variables
day=`date +"%Y%m%d"`
time=`date +"%Hh%Mm%Ss"`
date="${day}_${time}"
readLength=100
organism="Human"
genome="GRCh38"
database="ensembl"

##### Directory paths
outputsDir=/home/faramir/data/neondisco-jet/outputData
logDir=/home/faramir/data/neondisco-jet/logDir
infoDir=/home/faramir/data/neondisco-jet/infoDir
starIndexesDir=/home/faramir/data/neondisco-jet/starIndexesDir
metadataDir=/home/faramir/data/neondisco-jet/metadataDir
ErrorDir=/home/faramir/data/neondisco-jet/ErrorDir

mkdir -p ${logDir}

##### Input files
infoFile=/home/faramir/repos/neondisco-jet/infoDir/samples.txt
fastaFile=/home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtfGeneFile=/home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf

##### Log file
logFile="${logDir}/STARrun_${date}.log"
touch ${logFile}


#-------------------------------------------------------
##-- STEP 0: STAR genome indexing (run once per genome)
#-------------------------------------------------------

threads=8
cmd="${starBinDir}/STAR \
    --runThreadN ${threads} \
    --runMode genomeGenerate \
    --genomeDir ${starIndexesDir} \
    --genomeFastaFiles ${fastaFile} \
    --sjdbGTFfile ${gtfGeneFile} \
    --sjdbOverhang ${readLength}"
echo ${cmd}
#eval "time ${cmd}"


#-------------------------------------------------------
##-- STEP 1: STAR alignment (per sample loop)
#-------------------------------------------------------

# metadata columns: R1 path, R2 path, sample name
while read fastqR1Filegz fastqR2Filegz name day;
do
    echo -e "\e[1m${name}\e[0m" >> ${logFile}

    ############################################################################
    ##### Sample-specific output directory

    outputSampleDir="${outputsDir}/${name}_${day}"
    mkdir -p ${outputSampleDir}

    tmpDir=${ErrorDir}

    echo -e "\e[1m${name}\t Creating output files and setting tmpDir\e[0m" >> ${logFile}

    ############################################################################
    ##### Sample-specific input and output file paths

    # Derive uncompressed path from gz path using bash parameter expansion
    fastqR1File="${fastqR1Filegz%.gz}"
    fastqR2File="${fastqR2Filegz%.gz}"

    echo -e "\e[1m${name}\t Reading fastq files\t${fastqR1Filegz} and ${fastqR2Filegz}\e[0m" >> ${logFile}

    prefix="${outputSampleDir}/${name}"

    samFile="${prefix}_Chimeric.out.sam"
    bamFile="${prefix}_Aligned.sortedByCoord.out.bam"
    bamChimericFile="${prefix}_Chimeric.out.bam"
    bamChimericSortFile="${prefix}_Chimeric.out.sort.bam"
    chimericFile="${prefix}_Chimeric.out.junction"
    junctionFile="${prefix}_SJ.out.tab"

    echo -e "\e[1m${name}\t Output prefix: ${prefix}\e[0m" >> ${logFile}

    ############################################################################
    ##### Check input files exist

    flag=false

    # Decompress fastq files if only the gz version exists
    if [ ! -f ${fastqR1File} -a -f ${fastqR1Filegz} ] ; then gunzip -c $fastqR1Filegz > $fastqR1File; fi
    if [ ! -f ${fastqR2File} -a -f ${fastqR2Filegz} ] ; then gunzip -c $fastqR2Filegz > $fastqR2File; fi

    # Verify at least one version of each file exists
    if [ ! -f ${fastqR1File} -a ! -f ${fastqR1Filegz} ] ; then echo "${fastqR1File} not found !" 1>&2 ; fastqR1File="NA" ; flag=true ; fi
    if [ ! -f ${fastqR2File} -a ! -f ${fastqR2Filegz} ] ; then echo "${fastqR2File} not found !" 1>&2 ; fastqR2File="NA" ; flag=true ; fi

    if ${flag} ; then continue; fi
    echo -e "${name}\t${fastqR1File}\t${fastqR2File}" >> ${logFile}

    ############################################################################
    ##### STAR alignment

    threads=16
    cmd="${starBinDir}/STAR \
            --quantMode GeneCounts \
            --twopassMode Basic \
            --runThreadN ${threads} \
            --genomeDir ${starIndexesDir} \
            --sjdbGTFfile ${gtfGeneFile} \
            --sjdbOverhang 100 \
            --readFilesIn ${fastqR1File} ${fastqR2File} \
            --outFileNamePrefix ${prefix}_ \
            --outTmpDir ${tmpDir}/STAR_${name}_${date} \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --bamRemoveDuplicatesType UniqueIdentical \
            --outFilterMismatchNoverLmax 0.04 \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 1000 \
            --winAnchorMultimapNmax 1000 \
            --chimOutType WithinBAM \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 ; \
        ${samtoolsBinDir}/samtools view -@ ${threads} -b ${samFile} > ${bamChimericFile} ; \
        ${samtoolsBinDir}/samtools sort -@ ${threads} -o ${bamChimericSortFile} -O bam ${bamChimericFile} ; \
        ${samtoolsBinDir}/samtools index ${bamChimericSortFile} ; \
        ${samtoolsBinDir}/samtools index ${bamFile}"

    echo ${cmd}
    eval "${cmd}"

    date=`date +"%Y%m%d_%Hh%Mm%Ss"`
    echo -e "${name}\tSTAR alignment:\t${cmd}\t${date}" >> ${logFile}

done < ${infoFile}
