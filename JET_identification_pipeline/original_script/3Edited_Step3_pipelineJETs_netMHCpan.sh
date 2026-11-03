#!/bin/bash

#### created by Alexandre Houy and Christel Goudot
#### modified and adapted by Muhammad Faramir

#-------------------------------------------------------
##-- Usage and Argument Validation
#-------------------------------------------------------

usage() {
    echo ""
    echo "Usage: $(basename "$0") <SAMPLE_NAME> <HLA_ALLELES>"
    echo ""
    echo "  SAMPLE_NAME  : Sample identifier (e.g. 379T)"
    echo "  HLA_ALLELES  : Comma-separated HLA alleles in netMHCpan format"
    echo "                 (e.g. HLA-A11:01,HLA-A02:03,HLA-B13:01)"
    echo ""
    echo "Example:"
    echo "  $(basename "$0") 379T HLA-A11:01,HLA-A02:03,HLA-B13:01,HLA-B38:02,HLA-C07:02,HLA-C03:04"
    echo ""
    exit 1
}

# Check number of arguments
if [ "$#" -ne 2 ]; then
    echo "Error: Expected 2 arguments, got $#."
    usage
fi

SAMPLE_NAME=$1
HLA_ALLELES=$2

# Validate SAMPLE_NAME
if [ -z "$SAMPLE_NAME" ]; then
    echo "Error: SAMPLE_NAME cannot be empty."
    usage
fi

# Validate HLA_ALLELES is not empty
if [ -z "$HLA_ALLELES" ]; then
    echo "Error: HLA_ALLELES cannot be empty."
    usage
fi

# Validate HLA_ALLELES format — must contain at least one HLA- prefix
if ! echo "$HLA_ALLELES" | grep -qE "HLA-[ABC]"; then
    echo "Error: HLA_ALLELES does not appear to contain valid HLA alleles."
    echo "       Expected format: HLA-A11:01,HLA-B13:01,HLA-C07:02"
    echo "       Got: ${HLA_ALLELES}"
    usage
fi


#-------------------------------------------------------
##-- Configuration
#-------------------------------------------------------

##### Software paths
netMHCpan4Dir=/home/faramir/repos/neondisco-jet/dependencies/netMHCpan-4.2

##### Run variables
day=`date +"%Y%m%d"`
time=`date +"%Hh%Mm%Ss"`
date="${day}_${time}"
organism="Human"
genome="GRCh38"

##### Directory paths
outputsDir=/home/faramir/data/neondisco-jet/outputData
logDir=/home/faramir/data/neondisco-jet/logDir

mkdir -p ${logDir}

##### Log file
logFile="${logDir}/netMHCpan_run_${date}.log"
touch ${logFile}

echo "========================================" >> ${logFile}
echo "netMHCpan Step 3 — JET binding prediction" >> ${logFile}
echo "Sample     : ${SAMPLE_NAME}" >> ${logFile}
echo "HLA alleles: ${HLA_ALLELES}" >> ${logFile}
echo "Started    : $(date)" >> ${logFile}
echo "========================================" >> ${logFile}


#-------------------------------------------------------
##-- Find output directory automatically
#-------------------------------------------------------

# Check netMHCpan binary exists
if [ ! -f "${netMHCpan4Dir}/netMHCpan" ]; then
    echo "Error: netMHCpan binary not found at: ${netMHCpan4Dir}/netMHCpan"
    echo "       Please check netMHCpan4Dir path in the script."
    exit 1
fi

# Automatically find the most recent Step 2 output directory for this sample
outputSampleDir=$(ls -td ${outputsDir}/${SAMPLE_NAME}_* 2>/dev/null | head -1)

if [ -z "${outputSampleDir}" ]; then
    echo "Error: No output directory found for sample ${SAMPLE_NAME} in ${outputsDir}"
    echo "       Make sure Step 1 and Step 2 have been run before Step 3."
    exit 1
fi

echo "Output directory : ${outputSampleDir}" >> ${logFile}

# Verify the directory has FASTA files from Step 2
if ! ls ${outputSampleDir}/${SAMPLE_NAME}_Fusions*.fasta 1>/dev/null 2>&1; then
    echo "Error: No FASTA files found in ${outputSampleDir}"
    echo "       Make sure Step 2 has been run successfully before Step 3."
    exit 1
fi

prefix="${outputSampleDir}/${SAMPLE_NAME}"

echo "Prefix           : ${prefix}" >> ${logFile}


#-------------------------------------------------------
##-- STEP 3: MHC Class I Binding Prediction
#-------------------------------------------------------

echo -e "Starting netMHCpan binding prediction -----" >> ${logFile}

STARTTIME=$(date +%s)
overall_status=0

for size in `seq 8 11`
do
    echo "" >> ${logFile}
    echo "--- Peptide size: ${size} ---" >> ${logFile}

    # Input FASTA from Step 2
    fastaFile="${prefix}_Fusions_chim2.junc2.size${size}.fasta"

    # Validate FASTA exists
    if [ ! -f "${fastaFile}" ]; then
        echo "Warning: FASTA file not found for size ${size}: ${fastaFile}" | tee -a ${logFile}
        echo "         Skipping size ${size}."
        continue
    fi

    # Output files
    netmhcpanFile="${prefix}_Fusions_chim2.junc2.size${size}.netmhcpan4.txt"
    netmhcpanTmp="${netmhcpanFile}.tmp"

    echo "Input FASTA : ${fastaFile}" >> ${logFile}
    echo "Output file : ${netmhcpanFile}" >> ${logFile}

    # Run netMHCpan
    cmd="${netMHCpan4Dir}/netMHCpan \
        -a ${HLA_ALLELES} \
        -f ${fastaFile} \
        -l ${size} \
        -xls \
        -xlsfile ${netmhcpanFile} \
        -inptype 0 \
        > ${netmhcpanTmp}"

    echo "Command: ${cmd}" >> ${logFile}
    eval "${cmd}"

    if [ $? -eq 0 ]; then
        echo "netMHCpan size ${size}: SUCCESS" | tee -a ${logFile}
        date_done=`date +"%Y%m%d_%Hh%Mm%Ss"`
        echo -e "${SAMPLE_NAME}\tnetMHCpan\tsize=${size}\t${date_done}" >> ${logFile}
    else
        echo "Error: netMHCpan failed for size ${size}." | tee -a ${logFile}
        overall_status=1
    fi

done

ENDTIME=$(date +%s)
ELAPSED=$(($ENDTIME - $STARTTIME))

echo "" >> ${logFile}
echo "========================================" >> ${logFile}
echo "Finished: $(date)" >> ${logFile}
echo "Total runtime: ${ELAPSED} seconds" >> ${logFile}
echo "========================================" >> ${logFile}

if [ $overall_status -eq 0 ]; then
    echo "Step 3 completed successfully for sample ${SAMPLE_NAME}."
    echo "Output files in: ${outputSampleDir}"
else
    echo "Step 3 completed with errors. Check log: ${logFile}"
    exit 1
fi
