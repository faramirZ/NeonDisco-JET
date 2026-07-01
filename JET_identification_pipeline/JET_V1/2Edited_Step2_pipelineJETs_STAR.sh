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
organism="Human"
genome="GRCh38"
database="ensembl"

##### Directory paths
outputsDir=/home/faramir/data/neondisco-jet/outputData
logDir=/home/faramir/data/neondisco-jet/logDir
ErrorDir=/home/faramir/data/neondisco-jet/ErrorDir

mkdir -p ${logDir}

##### Path to JET R script
RscriptDir=/home/faramir/repos/neondisco-jet/JET_identification_pipeline

##### Input files
infoFile=/home/faramir/repos/neondisco-jet/infoDir/samples.txt

##### Log file
logFile="${logDir}/R_JET_run_${date}.log"
touch ${logFile}


#-------------------------------------------------------
##-- STEP 2: Identification and classification of junctions, selection of JETs
#-------------------------------------------------------

echo -e "Starting Step 2 -----" >> ${logFile}
echo -e "Date: ${date}" >> ${logFile}

# metadata columns: R1 path, R2 path, sample name (3 columns — no date)
while read fastqR1Filegz fastqR2Filegz name
do
    echo -e "\e[1m${name}\e[0m" >> ${logFile}

    ############################################################################
    ##### Find sample-specific output directory from Step 1 automatically

    outputSampleDir=$(ls -td ${outputsDir}/${name}_* 2>/dev/null | head -1)

    if [ -z "${outputSampleDir}" ]; then
        echo "Error: No output directory found for sample ${name} in ${outputsDir}" | tee -a ${logFile}
        echo "       Make sure Step 1 has been run before Step 2."
        continue
    fi

    echo -e "\e[1m${name}\t Found output directory: ${outputSampleDir}\e[0m" >> ${logFile}

    ############################################################################
    ##### Calling sample-specific data files (outputs from Step 1)

    prefix="${outputSampleDir}/${name}"

    chimericFile="${prefix}_Chimeric.out.junction"
    junctionFile="${prefix}_SJ.out.tab"
    logFinalOut="${prefix}_Log.final.out"

    # Validate input files exist
    if [ ! -f "${chimericFile}" ]; then
        echo "Error: Chimeric junction file not found: ${chimericFile}" | tee -a ${logFile}
        continue
    fi

    if [ ! -f "${junctionFile}" ]; then
        echo "Error: SJ junction file not found: ${junctionFile}" | tee -a ${logFile}
        continue
    fi

    if [ ! -f "${logFinalOut}" ]; then
        echo "Error: STAR log file not found: ${logFinalOut}" | tee -a ${logFile}
        continue
    fi

    # Extract library size from STAR log
    libsize="$(grep "Uniquely mapped reads number" ${logFinalOut} | sed -r 's/[[:space:]]+//g' | cut -d '|' -f2)"

    if [ -z "${libsize}" ]; then
        echo "Error: Could not extract libsize from ${logFinalOut}" | tee -a ${logFile}
        continue
    fi

    echo -e "\e[1m${name}\t Output prefix: ${prefix}\e[0m" >> ${logFile}
    echo -e "\e[1m${name}\t Library size: ${libsize}\e[0m" >> ${logFile}

    ############################################################################
    ##### Run R analysis for peptide sizes 8, 9, 10, 11

    for size in `seq 8 11`
    do
        echo -e "\e[1m${name}\t Running R analysis for size ${size}\e[0m" >> ${logFile}

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

        if [ $? -eq 0 ]; then
            echo -e "${name}\tsize=${size}\tR_analysis: SUCCESS" >> ${logFile}
        else
            echo -e "${name}\tsize=${size}\tR_analysis: FAILED" >> ${logFile}
        fi

        date_done=`date +"%Y%m%d_%Hh%Mm%Ss"`
        echo -e "${name}\tR_analysis\tsize=${size}\t${date_done}" >> ${logFile}

    done  # end size loop

done < ${infoFile}  # end sample loop

echo -e "END" >> ${logFile}
