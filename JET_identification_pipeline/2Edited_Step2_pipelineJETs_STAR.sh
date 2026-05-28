#!/bin/bash

#### created by Alexandre Houy and Christel Goudot
#### modified and adapted by Muhammad Faramir

#-------------------------------------------------------
##-- Configuration
#-------------------------------------------------------

##### Software paths
RprojectDir=/opt/R/4.5.3/bin

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

##### Path to JET R script
RscriptDir=/home/faramir/repos/neondisco-jet/JET_identification_pipeline

##### Input files
infoFile=/home/faramir/repos/neondisco-jet/infoDir/samples.txt

##### Log file
logFile="${logDir}/R_netMHCpan_run_${date}.log"
touch ${logFile}


#-------------------------------------------------------
##-- STEP 2: Identification and classification of junctions, selection of JETs
#-------------------------------------------------------

echo -e "Starting R -----" >> ${logFile}

# metadata columns: R1 path, R2 path, sample name
while read fastqR1Filegz fastqR2Filegz name day
do
    echo -e "\e[1m${name}\e[0m" >> ${logFile}

    ############################################################################
    ##### Sample-specific output directory

    outputSampleDir="${outputsDir}/${name}_${day}"
    mkdir -p ${outputSampleDir}

    echo -e "\e[1m${name}\t Creating output files\e[0m" >> ${logFile}

    ############################################################################
    ##### Calling sample-specific data files (outputs from Step 1)

    prefix="${outputSampleDir}/${name}"

    samFile="${prefix}_Chimeric.out.sam"
    bamFile="${prefix}_Aligned.sortedByCoord.out.bam"
    bamChimericFile="${prefix}_Chimeric.out.bam"
    bamChimericSortFile="${prefix}_Chimeric.out.sort.bam"
    chimericFile="${prefix}_Chimeric.out.junction"
    junctionFile="${prefix}_SJ.out.tab"

    logFinalOut="${prefix}_Log.final.out"
    libsize="$(grep "Uniquely mapped reads number" ${logFinalOut} | sed -r 's/[\\t]+//g' | cut -d '|' -f2)"

    echo -e "\e[1m${name}\t Calling and naming output files:\e[0m" >> ${logFile}
    echo -e "\e[1m${name}\t\t ${prefix}\e[0m" >> ${logFile}
    echo -e "\e[1m${name}\t\t ${libsize}\e[0m" >> ${logFile}

    size=11

    ############################################################################
    ##### Naming sample-specific output files

    fastaFile="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.fasta"
    idsFile="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.ids.txt"
    netmhcpanFile4="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.netmhcpan4.0.txt"

    ############################################################################
    ##### R analysis

    cmd="${RprojectDir}/Rscript ${RscriptDir}/Copy_JET_analysis_filtered.R \
        --chimeric ${chimericFile} \
        --junction ${junctionFile} \
        --genome ${genome} \
        --size ${size} \
        --libsize ${libsize} \
        --prefix ${prefix} \
        --verbose"

    echo ${cmd}
    eval "${cmd}"

    date=`date +"%Y%m%d_%Hh%Mm%Ss"`
    echo -e "${name}\tR_analysis:\t${cmd}\t${date}" >> ${logFile}

done < ${infoFile}

echo -e "END" >> ${logFile}
