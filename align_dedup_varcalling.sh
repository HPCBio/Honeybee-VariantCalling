#!/bin/bash
#
# align_dedup_varcall.sh <runfile> <outputdir> <read1> <read2> <log.in> <log.ou> <qsub>
# 
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 7 ]
then
        MSG="Parameter mismatch"
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
        echo -e "#####################################################################################"        
        echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
        echo -e "#####################################################################################"        

        echo -e "\n\n#####################################################################################"        
        echo -e "#############             DECLARING VARIABLES                         ###############"
        echo -e "#####################################################################################\n\n"        
        
        umask 0027
	set -x
	echo `date`
        scriptfile=$0
        runfile=$1
        outputdir=$2
        R1=$3
        R2=$4
        elog=$5
        olog=$6
        qsubfile=$7
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Honeybee Variant Calling Workflow failure message" "$redmine"
	    exit 1;
        fi
        
        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
	rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )        
        aligner_mod=$( cat $runfile | grep -w ALIGNERMODULE | cut -d '=' -f2  )
        aligner_parms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
        picard_mod=$( cat $runfile | grep -w PICARDMODULE | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samtools_mod=$( cat $runfile | grep -w SAMTOOLSMODULE | cut -d '=' -f2 )
        samblaster_mod=$( cat $runfile | grep -w SAMBLASTERMODULE | cut -d '=' -f2 )
        sorttool_mod=$( cat $runfile | grep -w SORTMODULE | cut -d '=' -f2 )
        markduplicates=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        gatk_mod=$( cat $runfile | grep -w GATKMODULE | cut -d '=' -f2 ) 
        gatk_dir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 ) 
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
        SampleName=`basename $outputdir`
        
        echo -e "\n\n##################################################################################"        
        echo -e "#############                       SANITY CHECK                   ###############"
        echo -e "##################################################################################\n\n"        
 
        if [ ! -d $rootdir ]
        then
  	    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
        
        if [ ! -d  $refdir  ]
        then
 	    MSG="Invalid value specified for REFGENOMEDIR=$refdir in the configuration file."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        elif [ ! -s $refdir/$refgenome ]
        then
 	    MSG="$refdir/$refgenome reference genome file not found."
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;       
        fi

         
        if [ `expr ${#R1}` -lt 1]
        then
            MSG="$R1 read one file not found"
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
            exit 1
        elif [ ! -s $R1 ]
        then
            MSG="$R1 read one file not found"
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
            exit 1                
        fi

        if [ `expr ${#R2}` -lt 1]
        then
            MSG="$R2 read two file not found"
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
            exit 1
        elif [ ! -s $R2 ]
        then
            MSG="$R2 read two  file not found"
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
            exit 1                
        fi

        if [ ! -d $outputdir ]
        then
 	    MSG="$outputdir outputdir not found"
            echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
        
        if [ `expr ${#SampleName}` -lt 1]
        then
            MSG="$SampleName sample undefined variable"
            echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
            exit 1     
        else
            sID=$SampleName
	    sPU=$SampleName
	    sSM=$SampleName
	fi
        if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ] 
	then
	    MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit 1;
	fi
      
	RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )
        rgheader=$( echo -n -e "@RG\t" )$( echo -e "${RGparms}"  | tr ":" "\t" | tr "=" ":" )

        if [ `expr ${#thr}` -lt 1 ]
        then
            thr=$PBS_NUM_PPN
        fi
        
        if [ `expr ${#markduplicates}` -lt 1 ]
        then
            markduplicates="PICARD"
        fi


        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"                  
        echo -e "##################################################################################"        
        echo -e "#############  DECIDING WHICH STEPS TO SKIP                        ###############"
        echo -e "##################################################################################" 
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"          
        
        copy2node_skip="NO"
        realign_skip="NO"
        copy2projectdir_skip="NO"


        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"                  
        echo -e "##################################################################################"        
        echo -e "#############  COPY FILES FOR SAMPLE $SampleName  TO  COMPUTE NODE ###############"
        echo -e "##################################################################################" 
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"          
        
	echo `date`
	
        echo -e "############## this is where results will go on compute node    ##############"
        echo -e "############## they will be transferred at the end back to project folder  ###"
        mkdir -p /state/partition1/$SampleName
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
                MSG="setup folders in COMPUTE NODE failed exitcode=$exitcode."
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                exit $exitcode;
        fi  

        mkdir -p /state/partition1/$SampleName/align
        mkdir -p /state/partition1/$SampleName/realign
        mkdir -p /state/partition1/$SampleName/variant

	if [ $copy2node_skip == "NO" ]
	then
	    echo -e "############## this is where reference genome will go in compute node  /dev/shm/${USER}/refdir ######"        
            mkdir -p /dev/shm/${USER}/refdir
            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="setup folders in COMPUTE NODE failed exitcode=$exitcode."
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                exit $exitcode;
            fi 
            
            cp $refdir/* /dev/shm/${USER}/refdir/
            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="copy $refdir to COMPUTE NODE failed exitcode=$exitcode."
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                exit $exitcode;
            elif [ "$(ls -A /dev/shm/${USER}/refdir)" ]
            then
                echo "#####  copy to compute node was successful. Create soft links  in compute node         ########"
                chmod -r a+r,w /dev/shm/${USER}/refdir/
                cd /dev/shm/${USER}/refdir/
                ln -s /dev/shm/${USER}/refdir/$refgenome genome.fa
                ref_local=/dev/shm/${USER}/refdir/genome.fa
            else
                MSG="copy $refdir to COMPUTE NODE failed exitcode=$exitcode."
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	        rm -r /dev/shm/${USER}/refdir
	        rm -r state/partition1/$SampleName/	        
                exit $exitcode;            
            fi          
        else
	    echo -e "######## copying reference genome to  compute node  was skipped. Using original folder  ###########"
            ref_local=${refdir}/genome.fa	    
        fi

	echo `date`	

        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"          	
	echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"          
	
	module load $samtools_mod
	module load $sorttool_mod
	SampleDir_local=/state/partition1/$SampleName
        AlignDir=/state/partition1/$SampleName/align
        RealignDir=/state/partition1/$SampleName/realign
        VarcallDir=/state/partition1/$SampleName/variant


        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"          
        echo -e "##################################################################################"        
        echo -e "#############   ALIGN-DEDUPPLICATION  FOR SAMPLE $SampleName       ###############"
        echo -e "##################################################################################" 
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"          

	echo `date` 
	

        cd $AlignDir
        rm  $AlignDir/*
	module load $aligner_mod

	echo -e "\n\n##################################################################################"
	echo -e "#############   select the dedup tool and then run the command        ############"
        echo -e "##################################################################################\n\n"          



        if [ $markduplicates == "SAMBLASTER" ]
        then
		echo -e "\n\n##################################################################################"
		echo -e "##    CASE1: dedup tool is $markduplciates and a single for align-deduplication ##"  
		echo -e "##################################################################################\n\n"

		module load $samblaster_mod

		echo -e "\n\n##################################################################################"
		echo -e "############# command one: alignment and deduplication                ############"	     
		echo -e "##################################################################################\n\n"
		bwa mem $aligner_parms -t $thr -R "${rgheader}" $ref_local $R1 $R2 | samblaster | samtools view -@ $thr -bSu -> ${SampleName}.wdups.bam 
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
                    MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	            
	            echo -e "#########       cleanup compute node and exit               ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit $exitcode;
		fi

		echo -e "\n\n##################################################################################"	     
		echo -e "#############  making sure that a file was produced with alignments       ############"
		echo -e "##################################################################################\n\n"
		
                if [ -s $AlignDir/${SampleName}.wdups.bam ]
                then            
                    echo -e "### sometimes we may have a BAM file with NO alignmnets, just the header ##"
                    numAlignments=$( samtools view -c $AlignDir/${SampleName}.wdups.bam ) 

                    echo `date`
                    if [ $numAlignments -eq 0 ]
                    then
                        MSG="bwa mem command did not produce alignments for $AlignDir/${SampleName}.wdups.bam. alignment failed"
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	       
	                echo -e "#########       cleanup compute node and exit                 ############"
	                cp -R /state/partition1/$SampleName $rootdir
	                #rm -r /dev/shm/${USER}/refdir
	                rm -r state/partition1/$SampleName/	        
                        exit 1;
                    else
                        echo -e "####### $AlignDir/${SampleName}.wdups.bam seems to be in order ###########"
                    fi
                else 
                    MSG="bwa mem command did not produce a file $AlignDir/${SampleName}.wdups.bam. alignment failed"
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	            echo -e "#########       cleanup compute node and exit                   ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit 1;          
                fi       
                
                
		echo -e "\n\n##################################################################################"
 		echo -e "#############  command two: sort                                      ############"
		echo -e "##################################################################################\n\n"

		novosort --index --tmpdir ./ --threads $thr -m 16g --compression 1 -o ${SampleName}.wdups.sorted.bam ${SampleName}.wdups.bam
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
                    MSG="alignment step  failed during sorting for sample $SampleName exitcode=$exitcode."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	            
	            echo -e "#########       cleanup compute node and exit               ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit $exitcode;
		fi
		module unload $samblaster_mod               

		echo -e "\n\n##################################################################################"
 		echo -e "#############      END SAMBLASTER BLOCK                               ############"
		echo -e "##################################################################################\n\n"             

        elif  [ $markduplicates == "PICARD" ]
        then
		echo -e "\n\n##################################################################################"
		echo -e "CASE2: dedup tool is $markduplciates one cmd for align and one for deduplication ##"  
		echo -e "##################################################################################\n\n"

	       
		module load $picard_mod

		echo -e "\n\n##################################################################################"	     
		echo -e "#############  command one: alignment                                 ############"
		echo -e "##################################################################################\n\n"

		bwa mem $aligner_parms -t $thr -R "${rgheader}" $ref_local $R1 $R2 | samtools view -@ $thr -bSu -> ${SampleName}.nodups.bam 
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
                    MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	            
	            echo -e "#########       cleanup compute node and exit              ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit $exitcode;
		fi   

		echo -e "\n\n##################################################################################"	     
		echo -e "#############  making sure that a file was produced with alignments       ############"
		echo -e "##################################################################################\n\n"
		
                if [ -s $AlignDir/${SampleName}.nodups.bam ]
                then            
                    echo -e "### sometimes we may have a BAM file with NO alignmnets, just the header ##"
                    numAlignments=$( samtools view -c $AlignDir/${SampleName}.nodups.bam ) 

                    echo `date`
                    if [ $numAlignments -eq 0 ]
                    then
                        MSG="bwa mem command did not produce alignments for $AlignDir/${SampleName}.nodups.bam. alignment failed"
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	       
	                echo -e "#########       cleanup compute node and exit                 ############"
	                cp -R /state/partition1/$SampleName $rootdir
	                #rm -r /dev/shm/${USER}/refdir
	                rm -r state/partition1/$SampleName/	        
                        exit 1;
                    else
                        echo -e "####### $AlignDir/${SampleName}.wdups.bam seems to be in order ###########"
                    fi
                else 
                    MSG="bwa mem command did not produce a file $AlignDir/${SampleName}.nodups.bam. alignment failed"
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	            echo -e "#########       cleanup compute node and exit                   ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit 1;          
                fi       
		echo -e "\n\n##################################################################################"	     
		echo -e "#############  command two: sort                                      ############"
		echo -e "##################################################################################\n\n"
		novosort --index --tmpdir ./ --threads $thr -m 16g --compression 1 -o ${SampleName}.nodups.sorted.bam ${SampleName}.nodups.bam
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
                    MSG="alignment step  failed during sorting for sample $SampleName exitcode=$exitcode."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 

	            echo -e "#########       cleanup compute node and exit              ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit $exitcode;
		fi 

		echo -e "\n\n##################################################################################"	     
		echo -e "############# command three: deduplication                            ############"
		echo -e "##################################################################################\n\n"     
	     
		java -Xmx8g -Xms1024m -jar $picardir/MarkDuplicates.jar \
                    INPUT=${SampleName}.nodups.sorted.bam \
                    OUTPUT=${SampleName}.wdups.sorted.bam \
                    TMP_DIR=$AlignDir ASSUME_SORTED=true MAX_RECORDS_IN_RAM=null \
                    CREATE_INDEX=true METRICS_FILE=${SampleName}.dup.metric VALIDATION_STRINGENCY=SILENT 
                    
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
                    MSG="alignment step  failed during deduplication for sample $SampleName exitcode=$exitcode."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 

	            echo -e "#########       cleanup compute node and exit              ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit $exitcode;
		fi 
		module unload $picard_mod
		module unload java         # this version of java clashes with the version needed to run gatk
		echo `date`

		echo -e "\n\n##################################################################################"	       
		echo -e "#############   END PICARD BLOCK                                      ############" 
		echo -e "##################################################################################\n\n"  
        else
                MSG="unrecognized deduplication tool $markduplicates"
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 

	        echo -e "#########       cleanup compute node and exit                   ############"
	        cp -R /state/partition1/$SampleName $rootdir
	        #rm -r /dev/shm/${USER}/refdir
	        rm -r state/partition1/$SampleName/	        
                exit 1;        

        fi

        echo -e "\n\n##################################################################################"
        echo -e "#############             STILL  INSIDE ALIGNMENT BLOCK               ############"
        echo -e "#############  CHECK THAT A DEDUPLICATED BAM FILE WAS CREATED         ############"
        echo -e "##################################################################################\n\n"         

        if [ -s $AlignDir/${SampleName}.wdups.sorted.bam ]
        then            
               echo -e "### sometimes we may have a BAM file with NO alignmnets, just the header ##"
               numAlignments=$( samtools view -c $AlignDir/${SampleName}.wdups.sorted.bam ) 
               echo `date`
               if [ $numAlignments -eq 0 ]
               then
                    MSG="bwa mem command did not produce alignments for $AlignDir/${SampleName}.wdups.sorted.bam. alignment failed"
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	       
	            echo -e "#########       cleanup compute node and exit                 ############"
	            cp -R /state/partition1/$SampleName $rootdir
	            #rm -r /dev/shm/${USER}/refdir
	            rm -r state/partition1/$SampleName/	        
                    exit 1;
                else
                    echo -e "####### $AlignDir/${SampleName}.wdups.sorted.bam seems to be in order ###########"
                fi
        else 
                MSG="bwa mem command did not produce a file $AlignDir/${SampleName}.wdups.sorted.bam. alignment failed"
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	        echo -e "#########       cleanup compute node and exit                   ############"
	        cp -R /state/partition1/$SampleName $rootdir
	        #rm -r /dev/shm/${USER}/refdir
	        rm -r state/partition1/$SampleName/	        
                exit 1;          
        fi       
         
	module unload $aligner_mod 
	echo `date`	   


	echo -e "\n\n##################################################################################"
        echo -e "#############     END ALIGNMENT-DEDUPLICATION BLOCK                   ############"
	echo -e "##################################################################################\n\n"


        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"          
	echo -e "##################################################################################"        
        echo -e "#############   GATK REALIGNMENT    FOR SAMPLE $SampleName         ###############"
        echo -e "##################################################################################"
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"
        
	echo `date`        
        if [ $realign_skip = "YES" ]
        then
                cd $RealignDir
                ln $AlignDir/${SampleName}.wdups.sorted.bam   ${SampleName}.realigned.bam        
        else
                cd $RealignDir
                rm -r $RealignDir/*           
                module load $gatk_mod

	        echo -e "\n\n##################################################################################"
                echo -e "########### command one: executing GATK RealignerTargetCreator          ##########" 
	        echo -e "##################################################################################\n\n"

                java -Xmx8g -Xms1024m -Djava.io.tmpdir=$RealignDir -jar $gatk_dir/GenomeAnalysisTK.jar \
 	         -R $ref_local \
 	         -I $AlignDir/${SampleName}.wdups.sorted.bam \
 	         -T RealignerTargetCreator \
                 -nt $thr \
 	         -o ${SampleName}.realignTargetCreator.list
 
 	        exitcode=$?
 	        echo `date`
 	        if [ $exitcode -ne 0 ]
 	        then
 	                 MSG="realignertargetcreator command failed exitcode=$exitcode. realignment for sample $SampleName stopped"
 	                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	                 echo -e "#########       cleanup compute node and exit                   ############"
	                 cp -R /state/partition1/$SampleName $rootdir
	                 #rm -r /dev/shm/${USER}/refdir
	                 rm -r state/partition1/$SampleName/	        
 	                 exit $exitcode;
                elif [ ! -s ${SampleName}.realignTargetCreator.list ]
 	        then
 	                 MSG="${SampleName}.realignTargetCreator.list file not created. realignment for sample $SampleName stopped"
 	                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"              
 	                 exit 1;
 	        fi

	        echo -e "\n\n##################################################################################" 
                echo -e "########### command two: executing GATK IndelRealigner and generating output #####"
	        echo -e "##################################################################################\n\n"

 	        java -Xmx8g -Xms1024m -Djava.io.tmpdir=$RealignDir -jar $gatk_dir/GenomeAnalysisTK.jar \
 	         -R $ref_local \
 	         -I $AlignDir/${SampleName}.wdups.sorted.bam \
 	         -T IndelRealigner \
 	         -targetIntervals ${SampleName}.realignTargetCreator.list \
 	         -o ${SampleName}.realigned.bam
  
 	        exitcode=$?
 	        echo `date`
 	        if [ $exitcode -ne 0 ]
 	        then
 	                 MSG="indelrealigner command failed exitcode=$exitcode. realignment for sample $SampleName stopped"
 	                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	                 echo -e "#########       cleanup compute node and exit                   ############"
	                 cp -R /state/partition1/$SampleName $rootdir
	                 #rm -r /dev/shm/${USER}/refdir
	                 rm -r state/partition1/$SampleName/	        
 	                 exit $exitcode;
 	        elif [ ! -s ${SampleName}.realigned.bam ]
 	        then
 	                 MSG="${SampleName}.realigned.bam  indelrealigner file not created. realignment for sample $SampleName stopped"
 	                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	                 echo -e "#########       cleanup compute node and exit                   ############"
	                 cp -R /state/partition1/$SampleName $rootdir
	                 #rm -r /dev/shm/${USER}/refdir
	                 rm -r state/partition1/$SampleName/	        
 	                 exit 1;
                fi	
                module unload $gatk_mod
	        echo `date`             

	        echo -e "\n\n##################################################################################"
                echo -e "#############           END REALIGNMENT BLOCK                         ############"
	        echo -e "##################################################################################\n\n"          
        fi

        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"	
        echo -e "##################################################################################"        
        echo -e "#############    GATK VARIANT CALLING   FOR SAMPLE $SampleName     ###############"
        echo -e "##################################################################################"
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"
        
	echo `date`        

        cd $VarcallDir
        rm -r $VarcallDir/*
        module load $gatk_mod

	echo -e "\n\n##################################################################################"            
        echo -e "########### command one: executing GATK HaplotypeCaller command         ##########" 
	echo -e "##################################################################################\n\n"

        java -Xmx6g -Xms1g -Djava.io.tmpdir=./ -jar $gatk_dir/GenomeAnalysisTK.jar \
                 -T HaplotypeCaller \
                 -R $ref_local \
                 -I $RealignDir/${SampleName}.realigned.bam \
                 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
                 -gt_mode DISCOVERY \
                 -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
                 -stand_call_conf 30 \
                 -stand_emit_conf 30 \
                 --sample_ploidy 1 \
                 -nt 1 -nct $thr \
                 -o ${SampleName}.rawVariants.vcf

  	exitcode=$?
 	echo `date`
 	if [ $exitcode -ne 0 ]
 	then
 	        MSG="haplotypecaller command failed exitcode=$exitcode. variant calling for sample $SampleName stopped"
 	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	        echo -e "#########       cleanup compute node and exit                   ############"
	        cp -R /state/partition1/$SampleName $rootdir
	        #rm -r /dev/shm/${USER}/refdir
	        rm -r state/partition1/$SampleName/	        
 	                 exit $exitcode;
 	elif [ ! -s ${SampleName}.rawVariants.vcf ]
 	then
 	        MSG="${SampleName}.rawVariants.vcf  file not created. variant calling for sample $SampleName stopped"
 	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	        echo -e "#########       cleanup compute node and exit                   ############"
	        cp -R /state/partition1/$SampleName $rootdir
	        #rm -r /dev/shm/${USER}/refdir
	        rm -r state/partition1/$SampleName/	         	       
 	        exit 1;
        fi	           
        module unload $gatk_mod
 	echo `date`           

	echo -e "\n\n##################################################################################"
        echo -e "#############       END VARIANT CALLING BLOCK                         ############"        
	echo -e "##################################################################################\n\n"
 

        echo -e "\n\n##################################################################################"  
        echo -e "##################################################################################"		
        echo -e "##################################################################################"        
        echo -e "#############    COPYING RESULTS, CLEANUP COMPUTE NODE AND WRAP-UP             ###"
        echo -e "#############    Destination of results: $outputdir/$SampleName                ###"        
        echo -e "##################################################################################"
        echo -e "##################################################################################"  
        echo -e "##################################################################################\n\n"	
	module unload $samtools_mod
	module unload $sorttol_mod
        
	echo `date`
        if [ $copy2projectdir_skip == "NO" ]
        then
                echo -e "#############    Source of results:      /state/partition1/$SampleName  #######"        
                echo -e "#############    Destination of results: $rootdir/$SampleName         #########"
            
                cp -R /state/partition1/$SampleName $rootdir
            
  	        exitcode=$?
 	        echo `date`
 	        if [ $exitcode -ne 0 ]
 	        then
 	                 MSG="failed to copy results from /state/partition1/$SampleName to $outputdir "
 	                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	                 echo -e "#########       cleanup compute node and exit                   ############"
	                 #rm -r /dev/shm/${USER}/refdir
	                 rm -r state/partition1/$SampleName/	        
 	                 exit $exitcode;
 	        else
	                 echo -e "#########  COPY OF RESULTS FINISHED OK. cleanup compute node and exit    ###"
	                 #rm -r /dev/shm/${USER}/refdir
	                 rm -r state/partition1/$SampleName/	        
 	        fi
        fi
        
	echo `date`
	echo -e "\n\n##################################################################################"
        echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
	echo -e "##################################################################################\n\n"


