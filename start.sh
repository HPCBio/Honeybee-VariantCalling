#!/bin/bash
umask 0027
################################################################################################ 
# Program to calculate raw variants from honeybee samples of WGS short reads
# In order to run this pipeline please type at the command line
# start.sh <runfile>
################################################################################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 1 ]
then
        MSG="Parameter mismatch."
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "HoneyBee Variant Calling Workflow failure message" "$redmine"
        exit 1;
fi
        echo -e "\n\n########################################################################################"
        echo -e "#############                BEGIN VARIANT CALLING WORKFLOW              ###############"
        echo -e "########################################################################################\n\n"

        umask 0027
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found."
           exit 1;
        fi

        echo -e "\n\n########################################################################################"
        echo -e "#############                CHECKING PARAMETERS                         ###############"
        echo -e "########################################################################################\n\n"

        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )        
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        aligner_tool=$( cat $runfile | grep -w ALIGNERTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        aligner_mod=$( cat $runfile | grep -w ALIGNERMODULE | cut -d '=' -f2 )
        picard_mod=$( cat $runfile | grep -w PICARDMODULE | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samtools_mod=$( cat $runfile | grep -w SAMTOOLSMODULE | cut -d '=' -f2 )
        samblaster_mod=$( cat $runfile | grep -w SAMBLASTERMODULE | cut -d '=' -f2 )
        sorttool_mod=$( cat $runfile | grep -w SORTMODULE | cut -d '=' -f2 )
        markduplicates=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        gatk_mod=$( cat $runfile | grep -w GATKMODULE | cut -d '=' -f2 )        
        gatk_dir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        
        if [ -z $outputdir  ]
        then
 		MSG="Invalid value specified for OUTPUTDIR=$outputdir in the configuration file."
                echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
        fi

        if [ ! -d  $refdir  ]
        then
 		MSG="Invalid value specified for REFGENOMEDIR=$refdir in the configuration file."
                echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
        fi

        if [ ! -s  $refdir/$refgenome  ]
        then
 		MSG="Invalid value specified for REFGENOME=$refgenome in the configuration file."
                echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
        fi

        if [ ! -d  $picardir  ]
        then
 		MSG="Invalid value specified for PICARDIR=$picardir in the configuration file."
                echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
        fi

        if [ ! -d  $gatk_dir  ]
        then
 		MSG="Invalid value specified for GATKDIR=$gatk_dir in the configuration file."
                echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
        fi


        if [ $inputformat != "FASTQ"  ]
        then
            MSG="Incorrect value for INPUTFORMAT=$inputformat in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
        if [ $aligner_tool != "BWAMEM"  -a $aligner_tool != "BWA_MEM"  ]
        then
            MSG="Incorrect value for ALIGNERTOOL=$aligner_tool in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "Program $0 stopped at line=$LINENO.\n\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi
        if [ -z $email ]
        then
           MSG="Invalid value for parameter PBSEMAIL=$email in the configuration file"
           echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi

        if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ] 
	then
		MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
                echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	fi
        if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
        then
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
	    MSG="Invalid value for parameter INPUTTYPE=$input_type  in the configuration file."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
	if [ $analysis != "NORMAL" -a $analysis != "VARIANT" ]
	then
	    MSG="Invalid value for parameter ANALYSIS=$analysis  in the configuration file."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi

	if [ $markduplicates != "SAMBLASTER" -a $markduplicates != "PICARD" ]
	then
	    MSG="Invalid value for parameter MARKDUPLICATESTOOL=$markduplicates  in the configuration file."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
        
	if [ ! -s $sampleinfo ]
	then
	    MSG="SAMPLEINFORMATION=$sampleinfo  file not found."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
	
	if [ ! -d $sampledir ]
	then
	    MSG="INPUTDIR=$sampledir  directory not found."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
	
	numfiles=$( ls -1 $sampledir | wc -l )
	numsamples=$( wc -l $sampleinfo )	

	if [ $numfiles -lt 1 ]
	then
	    MSG="INPUTDIR=$sampledir  directory is empty."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;	
	fi

	if [ $numsamples -lt 1 ]
	then
	    MSG="SAMPLEINFORMATION=$sampleinfo  file is empty."
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;	
	fi
	
	if [ $paired -ne 1 ]
	then
	    MSG="Invalid value for parameter PAIRED=$paired in configuration file "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi


        echo -e "\n\n########################################################################################"
        echo -e "###########                      checking tools                       ##################"
        echo -e "########################################################################################\n\n"
        
        `module load $aligner_mod`
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
 	    MSG="Invalid value specified for ALIGNERMODULE=$aligner_mod in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        
        fi
	`module unload $aligner_mod`

        `module load $samblaster_mod`
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
 	    MSG="Invalid value specified for SAMBLASTERMODULE=$samtools_mod in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        
        fi
        `module unload $samblaster_mod`

        `module load $samtools_mod`
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
 	    MSG="Invalid value specified for SAMTOOLSMODULE=$samtools_mod in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        
        fi
        `module unload $samtools_mod`

        `module load $sorttool_mod`
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
 	    MSG="Invalid value specified for SORTMODULE=$sorttool_mod in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        
        fi
        `module unload $sorttool_mod`

        `module load $picard_mod`
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
 	    MSG="Invalid value specified for FASTQCMODULE=$fastqc_mod in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        
        fi
        `module unload $picard_mod`        
        `module unload java`
 
        `module load $gatk_mod`
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
 	    MSG="Invalid value specified for GATKMODULE=$gatk_mod in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        
        fi
        `module unload $gatk_mod`
        `module unload java`
        
        echo -e "\n\n########################################################################################"
        echo -e "#############       Everything seems ok. Now setup/configure output folder     #########"
        echo -e "########################################################################################\n\n"

        if [ ! -d $outputdir ]
        then
            mkdir -p $outputdir/logs
        else 
            echo "resetting directory"
	    `rm -r $outputdir/*`
            mkdir -p $outputdir/logs
        fi

	`cp $runfile $outputdir/runfile.txt`
	`cp $sampleinfo $outputdir/SAMPLENAMES.list`
	
        runfile=$outputdir/runfile.txt
        TopOutputLogs=$outputdir/logs

        echo -e "\n\n########################################################################################"
	echo -e "##### resetting the file with jobids, otherwise it will mess up the summary script ####"
        echo -e "########################################################################################\n\n"
        truncate -s 0 $TopOutputLogs/pbs.ANALYSIS

        echo -e "\n\n########################################################################################"
	echo -e "##### generate a qsub header so we would not have to repeat the same lines         ####"
        echo -e "########################################################################################\n\n"

	generic_qsub_header=$outputdir/qsubGenericHeader
	truncate -s 0 $generic_qsub_header
	echo "#!/bin/bash" > $generic_qsub_header
	echo "#PBS -V" >> $generic_qsub_header
	echo "#PBS -q $pbsqueue" >> $generic_qsub_header
	echo "#PBS -m ae" >> $generic_qsub_header
	echo "#PBS -M $email" >> $generic_qsub_header

	echo -e "##### let's check that it worked and that the file was created                     ####"
        if [ ! -s $generic_qsub_header ]
	then 
	    MSG="$generic_qsub_header is empty"
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit 1;
	fi

        echo -e "\n\n########################################################################################"
	echo -e "################### Documenting progress on redmine with this message ##################"
        echo -e "########################################################################################"
        echo -e "##### the first part of the Report also needs to be stored in Summary.Report      ######"
        echo -e "########################################################################################\n\n"

        truncate -s 0 $outputdir/logs/Summary.Report

        MSG="Variant calling workflow  started by username:$USER at: "$( echo `date` )
        LOGS="\nrunfile=$outputdir/runfile.txt \noutput directory=${outputdir} \nsample file=$sampleinfo"
        echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        echo -e "$MSG\n\nDetails:\n\n$LOGS" >> $outputdir/logs/Summary.Report


        echo -e "\n\n########################################################################################"
        echo -e "########################################################################################"
        echo -e "########################################################################################"
	echo -e "###################       MAIN LOOP STARTS HERE         ################################"
        echo -e "########################################################################################"
        echo -e "########################################################################################\n\n"
        
        while read sampleLine
        do
            if [ `expr ${#sampleLine}` -lt 1 ]
            then
                echo -e "\n\n########################################################################################"
                echo -e "##############                 skipping empty line        ##############################"
                echo -e "########################################################################################\n\n"
            else
                echo -e "\n\n########################################################################################"
                echo -e "##############         Processing next line $sampleLine   ##############################"
                echo -e "########################################################################################\n\n"
                
                sample=$( echo "$sampleLine" | cut -f 1 )  
                FQ_R1=$( echo "$sampleLine" | cut -f 2 )
                FQ_R2=$( echo "$sampleLine" | cut -f 3 )
 
                if [ `expr ${#sample}` -lt 1 ]
                then
                     MSG="unable to parse line $sampleLine"
                     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
                     exit 1
                fi

                if [ `expr ${#FQ_R1}` -lt 1 ]
                then
                     MSG="unable to parse line $sampleLine"
                     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
                     exit 1
                elif [ ! -s $sampledir/$FQ_R1 ]
                then
                     MSG="$sampledir/$FQ_R1 file not found"
                     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                                          
                     exit 1                
                fi

                if [ `expr ${#FQ_R2}` -lt 1 ]
                then
                     MSG="unable to parse line $sampleLine"
                     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
                     exit 1
                elif [ ! -s $sampledir/$FQ_R2 ]
                then
                     MSG="$sampledir/$FQ_R2 file not found"
                     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
                     exit 1                
                fi

                echo -e "\n\n########################################################################################"
                echo -e "###   Everything seems in order. Now creating folders where results will go  ###########"
                echo -e "########################################################################################\n\n"

                if [ -d $outputdir/$sample ]
                then
                     MSG="$outputdir/$sample already exists. Resetting it now."
                     rm -R $outputdir/$sample
                     mkdir -p $outputdir/$sample
                else 
                     mkdir -p $outputdir/$sample
                fi

                echo -e "\n\n########################################################################################"                
		echo -e "####   Launching Analysis script for SAMPLE $sample R1=$FQ_R1 R2=$FQ_R2       ##########"
                echo -e "########################################################################################\n\n"

		qsub1=$TopOutputLogs/qsub.align2vcall.$sample
		cat $outputdir/qsubGenericHeader > $qsub1
		echo "#PBS -N align2vcall.$sample" >> $qsub1
		echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
		echo "#PBS -o $TopOutputLogs/log.align2vcall.$sample.ou" >> $qsub1
		echo "#PBS -e $TopOutputLogs/log.align2vcall.$sample.in" >> $qsub1
		echo "$scriptdir/align_dedup_varcalling.sh $runfile $outputdir/$sample $sampledir/$FQ_R1 $sampledir/$FQ_R2 $TopOutputLogs/log.align2vcall.$sample.in $TopOutputLogs/log.align2vcall.$sample.ou $TopOutputLogs/qsub.align2vcall.$sample" >> $qsub1
		`chmod a+r $qsub1`               
		jobid=`qsub $qsub1` 
		echo $jobid >> $TopOutputLogs/pbs.ANALYSIS
		echo `date`
           fi 
        done <  $outputdir/SAMPLENAMES.list

        echo -e "\n\n########################################################################################"
        echo -e "########################################################################################"
	echo -e "#################           MAIN LOOP ENDS HERE                  #######################"
        echo -e "########################################################################################"
        echo -e "########################################################################################"
	echo -e "#################     Now, we need to generate summary           #######################"
        echo -e "########################################################################################"
        echo -e "########################################################################################\n\n"

        alljobids=$( cat $TopOutputLogs/pbs.ANALYSIS | sed "s/\..*//" | tr "\n" ":" )
        
        qsub2=$TopOutputLogs/qsub.summary
	cat $outputdir/qsubGenericHeader > $qsub2
	echo "#PBS -N Summary_vcall" >> $qsub2
	echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	echo "#PBS -o $TopOutputLogs/log.summary.ou" >> $qsub2
	echo "#PBS -e $TopOutputLogs/log.summary.in" >> $qsub2
	echo "#PBS -W depend=afterok:$alljobids " >> $qsub2
	echo "$scriptdir/summary.sh $runfile $TopOutputLogs/log.summary.in $TopOutputLogs/log.summary.ou $TopOutputLogs/qsub.summary" >> $qsub2
	`chmod a+r $qsub2`
	lastjobid=`qsub $qsub2`
        echo $lastjobid >> $TopOutputLogs/pbs.SUMMARY
	echo `date`     

        echo -e "\n\n########################################################################################"
	echo -e "##############                 EXITING NOW                            ##################"	
        echo -e "########################################################################################\n\n"
