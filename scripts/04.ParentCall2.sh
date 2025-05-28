#!/bin/bash
 echo LM3 ParentCall2 started at `date`
zcat ${DATA_GL} | java -cp ~/software/Lep_Map3_0.5/bin ParentCall2 data=${PEDIGREE} posteriorFile=- removeNonInformative=1 XLimit=2 |gzip >${RESULTS_DIR_PC2}/data.ParentCall_q20_Q20.gz 2>data.ParentCall.err
