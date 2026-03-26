#!/bin/bash

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/IDs`
echo "SampleID: ${sampleID}"

inputpath=/path/to/CX_reports

gzip -cd ${inputpath}/${sampleID}.CX_report.txt.gz | grep -P 'CG\t' > /path/out/CG_reports/${sampleID}.CG_report.txt
