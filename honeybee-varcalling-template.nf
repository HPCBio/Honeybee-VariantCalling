#!/projects/aces/grendon/apps/nextflow-0.27.3/nextflow
/*
HONEYBEE VARIANT CALLING PIPELINE
This is the generic nextflow template to process Honeybee WGS samples and generate raw variants 
using the GATK best practices protocol published
*/


/*Some hacks to deal with dates in output paths*/

import java.text.SimpleDateFormat

// version
version = 0.2


// Credit to Phil Ewels for this segment to generate the help message at the command-line

def helpMessage() {
    log.info"""
    =========================================
     HONEYBEE VARIANT CALLING PIPELINE v${version} 
    =========================================

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run -c <nxf-honeybee.conf>  <honeybee-template.nf>
    
    where nxf-honeybee.conf is the configuration file
    
    To override existing values from the command line, please type these parameters:


    I/O parameters
              --projectdir              path to project folder
              --batchid                 batchid-id. could be the sequencing lot number 
              --reads                   regular expression of files containing short reads. Include the full path. Example /full-path/*R{1,2}.fq
              --srcdir                  path to source folder
              --resultsdir              path to output folder
              --email                   email of user to receive emails with execution notifications
              --ticket                  redmine ticket to post execution notifications
              
    sample metadata
              --inputformat             format of short reads. Valid values 'FASTQ' 'FQ'
              --libtype                 Genomic library type. Valid values 'WGS' "WES'       
              --analysisoutput          Output of the analysis. Valid values 'RAWVARIANT' 'VARIANT'
              --paired                  Sequencing paired library. Valid values 'TRUE' 'YES' 'NO' 'FALSE'
              --skipTrim                Skip trimming step. Valid values 'TRUE' 'YES' 'NO' 'FALSE'
              --RGID                    RG header field ID. read group identifier. This value is dynamically assigned during execution.
              --RGLB                    RG header field LB. library identifier. This value is required.
              --RGPL                    RG header field PL. platform technology. This value is required.
              --RGPU                    RG header field PU. platform unit. This value is required.
              --RGSM                    RG header field SM. sample. This value is dynamically assigned during execution.
              --RGCN                    RG header field CN. sequencing center. This value is required.
              
    tools provenance
	      --trimTool                Name of trimming tool. Valid values 'trimmomatic' 'TRIMMOMATIC'
              --alignerTool             Name of aligner tool. Valid values 'BWA' 'bwa'
              --markdupTool             Name of deduplication tool. Valid values 'PICARD' 'picard'
              --vcallTool               Name of variant calling tool. Valid values 'GATK' 'gatk'
              --trimMod                 Trimmomatic module
              --bwaMod                  bwa module
              --picardMod               picard module
              --samtoolsMod             samtools module 
              --sortMod                 sort module
              --gatkMod                 gatk module
              --tabixMod                tabix module
              
    tool finetuning
              --trimdir                 path or environment variable for trimmomatic jar library 
              --picardir                path or environment variable for picard jar library
              --gatkdir                 path or environment variable for gatk jar library
              --bwaparams               bwa parameters
              --picardparams            picard parameters 
              --gatkparams              gatk parameters
              --JavaTrimmomatic         java parameters for executing a cmd with trimmomatic tool
              --JavaPicard              java parameters for executing a cmd with picard tool
              --JavaGATK                java parameters for executing a cmd with gatk tool
              
    paths to reference genome
              --bwarefgenome            bwa-indexed reference genome. Specify full path
              --refgenome               reference genome file in fasta format. Specify full path
              --dbsnp                   golden variants for reference genome in vcf format. Specify full path
              
     cluster resources
              --myQueue                 name of queue
              --alignTime               walltime per job for alignment step 
              --alignThreads            number of threads per job for alignment step 
              --alignMem                required memory per job for alignment step
              --dedupTime               walltime per job for deduplication step 
              --dedupThreads            number of threads per job for deduplication step
              --dedupMem                required memory per job for deduplication step
              --vcallTime               walltime per job for variant calling step 
              --vcallThreads            number of threads per job for variant calling step 
              --vcallMem                required memory per job for variant calling step 
              
    
    Other parameters that can be specified at the command line:
              --help                    To generate this message
	      -qs                       Maximum number of concurrent processes.  Valid value: a positive integer. Default value 1. Recommended value 20.  
              -with-report              Execution report generated by nextflow. Report in HTML format. Valid value:a file name
              -with-timeline            Timeline report of all processes executed during the execution. Report in HTML format. Valid value:a file name        
              -with-trace               Execution trace of each process. Report is a table and is called trace.txt

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

// Configurable variables

// To convert this script into a template, simply set the value to false for all variables that start with params. 

// I/O section
params.projectdir = false
params.batchid = false
params.reads = false
params.srcdir = false
params.resultsdir = false
params.email = false
params.ticket = false

// sample metadata
params.inputformat = false
params.libtype = false
params.analysisoutput = false
params.skipTrim = false
params.skipQC = "FALSE"
params.paired = false
params.RGID = false
params.RGLB = false
params.RGPL = false
params.RGPU = false
params.RGSM = false
params.RGCN = false

// tools provenance
params.trimTool = false
params.alignerTool = false
params.markdupTool = false
params.vcallTool = false
params.trimMod = false
params.bwaMod = false
params.picardMod = false
params.samtoolsMod = false
params.sortMod = false
params.gatkMod = false
params.tabixMod = false
params.RMod = false
params.vcftoolsMod = false

// tool finetuning
params.trimdir = false
params.picardir = false
params.gatkdir = false
params.bwaparams = false
params.picardparams = false
params.gatkparams = false
params.JavaTrimmomatic = false
params.JavaPicard = false
params.JavaGATK = false
params.adaptors = false
params.JavaPicardQC = false
trimOptions = " TOPHRED33 ILLUMINACLIP:${params.adaptors}:2:30:10:5:true LEADING:28 TRAILING:28 SLIDINGWINDOW:3:15 MINLEN:30 "


// paths to reference genome
params.bwarefgenome = false
params.refgenome= false
params.dbsnp = false

// cluster resources
params.executor = "pbs"
params.vcallQueue = "aces"
params.Queue = "secondary"
params.trimTime = "01:00:00"
params.trimThreads = '8'
params.trimMem = '10'
params.alignTime = '01:00:00'
params.alignThreads = '10'
params.alignMem = '100'
params.dedupTime = '01:00:00'
params.dedupThreads = '10'
params.dedupMem = '100'
params.vcallTime = '01:00:00'
params.vcallThreads = '10'
params.vcallMem = '100'


// Path for folders with temp and final results.
// Putting a timestamp on the name of the folder makes it hard to resume execution when it changes date
// timestamp            = new SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())
// outputdir            = "${params.resultsdir}/" + timestamp + "-batchid-${params.batchid}"
outputdir            = "${params.resultsdir}/batchid-${params.batchid}"
aligndir             = "${outputdir}/align"
dedupdir             = "${outputdir}/align"
realigndir           = "${outputdir}/realign"
vcalldir             = "${outputdir}/variant"
QCdir                = "${outputdir}/QC_reports"


// Validate inputs

if(  params.inputformat != 'FASTQ'  && params.inputformat != 'fastq' &&  params.inputformat != 'FQ'  &&  params.inputformat != 'fq' ){
    exit 1, "Invalid value for params.inputformat.  Valid options: 'FASTQ', 'fastq', 'FQ', 'fq'\n"
}

if(  params.libtype != 'WGS'  && params.libtype != 'wgs' &&  params.libtype != 'WES'  &&  params.libtype != 'wes' ){
    exit 1, "Invalid value for params.libtype.  Valid options: 'WGS', 'wgs', 'WES', 'wes'\n"
}

if(  params.analysisoutput != 'RAWVARIANT'  && params.analysisoutput != 'VARIANT' && params.analysisoutput != 'rawvariant'  && params.analysisoutput != 'variant'  ){
    exit 1, "Invalid value for params.analysisoutput.  Valid options: 'RAWVARIANT', 'rawvariant', 'VARIANT', 'variant'\n"
}

if(  params.paired != 'true'  && params.paired != 'TRUE' ){
    exit 1, "Invalid value for params.paired.  Valid options: 'TRUE', 'true'\n"
}

if(  params.skipTrim != 'true'  && params.skipTrim != 'TRUE' && params.skipTrim != 'false'  && params.skipTrim != 'FALSE' ){
    exit 1, "Invalid value for params.skipTrim.  Valid options: 'TRUE', 'true' 'FALSE' 'false'\n"
}

if(  params.skipTrim == 'false'  || params.skipTrim == 'FALSE' ){
    if( !params.trimTool || !params.trimMod || !params.trimdir ){
        exit 1, "Trimming tool was not correctly configured. Check params.trimTool params.trimMod params.trimdir\n"    
    }
    if( params.trimTool != 'trimmomatic' && params.trimTool != 'TRIMMOMATIC' ){
        exit 1, "Unknown Trimming tool\n"    
    }  
    params.skipTrim = 'FALSE'
    params.readprep      = " trimmomatic: ${trimOptions} "
} 

if(  params.skipTrim == 'true'  || params.skipTrim == 'TRUE' ){
    params.readprep      = " no trimming "
}

if( !params.RGLB || !params.RGPL || !params.RGPU || !params.RGCN ){
    exit 1, "Invalid value for RG values params.RGLB  params.RGPL  params.RGPU  params.RGCN\n"
}

if( !params.alignerTool || !params.markdupTool  || !params.vcallTool ){ 
    exit 1, "Invalid value for one of the following tools params.alignerTool params.markdupTool params.vcallTool\n"
}


if( params.alignerTool != 'BWA' &&  params.alignerTool != 'bwa' && params.alignerTool != 'BWAMEM' &&  params.alignerTool != 'bwamem' ){ 
    exit 1, "Invalid value for params.alignerTool. Valid values 'BWA' 'BWAMEM'\n"
}

if( params.markdupTool != 'PICARD' &&  params.markdupTool != 'picard' ){ 
    exit 1, "Invalid value for params.markdupTool. Valid values 'PICARD' 'picard'\n"
}

if( params.vcallTool != 'GATK' &&  params.vcallTool != 'gatk' ){ 
    exit 1, "Invalid value for params.vcallTool. Valid values 'GATK' 'gatk'\n"
}

if( !params.bwaMod || !params.picardMod  || !params.gatkMod || !params.samtoolsMod || !params.sortMod || !params.tabixMod ){ 
    exit 1, "Invalid value for one of the following modules params.bwaMod params.picardMod params.gatkMod params.samtoolsMod params.sortMod params.tabixMod \n"
}

if( !params.refgenome || !params.bwarefgenome ){
    exit 1, "Invalid value for reference genome params.refgenome OR params.bwarefgenome\n"
}

if( !params.reads || !params.resultsdir || !params.projectdir ){
    exit 1, "Invalid value for directory paths params.reads params.resultsdir params.projectdir\n"
}

// Redmine Ticket must be set for project mgmt tracking
if (params.ticket == 0) exit 1, "Must set Redmine ticket for pipeline summary to be sent"

// Create a channel for input read files. It should work with single reads and with paired reads

Channel
    .fromFilePairs( params.reads, size: -1, flat: true )
    .set { read_files_trimming  }


// Credit to Phil Ewels for this little bit of self-documentation
run_info = """
          
===============================================================================
                HONEYBEE VARIANT CALLING PIPELINE v${version}
 batchid ${params.batchid} and reads ${params.reads}
===============================================================================

Pipeline ver                : ${version}
Batchid                     : ${params.batchid}
Reads                       : ${params.reads}

Tools
Trimmomatic                 : ${params.trimMod}
BWA                         : ${params.bwaMod}
samtools                    : ${params.samtoolsMod}
picard                      : ${params.picardMod}
GATK                        : ${params.gatkMod}
tabix                       : ${params.tabixMod}
novosort                    : ${params.sortMod}
read preparation cmd        : ${params.readprep}
skip QC reporting           : ${params.skipQC}

Reference genome files
Reference genome            : ${params.refgenome}
DBSNP                       : ${params.dbsnp}

Paths and other options
Current user                : $USER
Current path                : $PWD
Project dir                 : ${params.projectdir}
Code dir                    : ${params.srcdir}
Output folder               : ${outputdir}
Ticket                      : ${params.ticket}
Cluster scheduler           : ${params.executor}

===============================================================================

"""

log.info run_info



/*
* Step 1. Quality Trimming and alignment with bwa
* Trimming is optional. The user must specify otherwise
*/

process TrimAlign {
	tag                  { id }
	executor             params.executor
	cpus                 params.alignThreads
	queue                params.Queue
	memory               "${params.alignMem} GB"
	time                 params.alignTime
	module               params.trimMod, params.bwaMod, params.samtoolsMod, params.sortMod
	errorStrategy        'finish'
	scratch              '/scratch'
	stageOutMode         'copy'
	publishDir           aligndir, mode: 'link'

	input:
	set val(id), file(read1), file(read2) from read_files_trimming
	file params.bwarefgenome

	output:
	file("*.sorted.bam") into alignChannel
	file("*.sorted.bam.bai") into alignIdxChannel    
	file "*.bam*"


	script:
	SampleName="${id}"    
	rgheader="@RG\\tID:${SampleName}\\tLB:${params.RGLB}\\tPL:${params.RGPL}\\tPU:${params.RGPU}\\tSM:${SampleName}\\tCN:${params.RGCN}"


	"""
	echo `hostname`
	echo "trim step"
	echo `date`
	
	if [  "${params.skipTrim}" = "FALSE" ]
	then


	       java -jar ${params.trimdir}/trimmomatic-0.36.jar PE \
	       -threads ${params.trimThreads}  ${read1} ${read2} \
	       ${id}_1.paired.fastq ${id}_1.unpaired.fastq \
	       ${id}_2.paired.fastq ${id}_2.unpaired.fastq $trimOptions

	else 

	       ln -s ${read1} ${id}_1.paired.fastq
	       ln -s ${read2} ${id}_2.paired.fastq

	fi

        echo "alignment step"
	echo `date`
	
	bwa mem ${params.bwaparams} -t ${params.alignThreads} -R \"${rgheader}\" ${params.bwarefgenome} ${id}_1.paired.fastq ${id}_2.paired.fastq | samtools view -@ ${params.alignThreads} -bSu -> ${SampleName}.nodups.unsorted.bam

	Qcheck1=\$( samtools quickcheck ${SampleName}.nodups.unsorted.bam )

	if [ `expr \${#Qcheck1}` -gt 0 ]
	then
	    echo "bwa produced corrupt output file ${SampleName}.nodups.unsorted.bam\n\n"
	    exit 1
	fi

	novosort --index --tmpdir /scratch --threads ${params.alignThreads} -m 16g --compression 1 -o ${SampleName}.nodups.sorted.bam ${SampleName}.nodups.unsorted.bam

	Qcheck2=\$( samtools quickcheck ${SampleName}.nodups.sorted.bam )

	if [ `expr \${#Qcheck2}` -gt 0 ]
	then
	    echo "novosort produced corrupt output file ${SampleName}.nodups.unsorted.bam \n\n"
	    exit 1
	fi

        echo "end alignment step"
	echo `date`
	
	"""
}

/*
* Step 2. MarkDuplicates
* MarkDuplicates with picard-tools
*/
process MarkDuplicates {
    tag                  { alignedBAM }    
    executor             params.executor
    cpus                 params.dedupThreads
    queue                params.Queue
    memory               "${params.dedupMem} GB"
    time                 params.dedupTime
    module               params.picardMod, params.samtoolsMod
    errorStrategy        'finish'
    scratch              '/scratch'
    stageOutMode         'copy'
    publishDir           dedupdir, mode: 'link'

    input:
    file alignedBAM from alignChannel
    file alignedBAI from alignIdxChannel
    
    output:
    file("*.wdups.sorted.bam") into dedupChannel
    file("*.wdups.sorted.bam") into QCBamChannel    
    file("*.wdups.sorted.bai") into dedupIdxChannel
    file("*.wdups.sorted.bai") into QCBamIdxChannel   
    file "*.wdups.sorted.bam*" 

    script:
    fileName="${alignedBAM.baseName}"
    BAM="${fileName.replace('.nodups.sorted','')}"
    
    """
	echo `hostname`

    
	java ${params.JavaPicard} -jar ${params.picardir}/picard.jar MarkDuplicates \
	INPUT=${alignedBAM} \
	METRICS_FILE=${BAM}.metrics \
	TMP_DIR=/scratch \
	COMPRESSION_LEVEL=0 \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=LENIENT \
	CREATE_INDEX=TRUE \
	OUTPUT=${BAM}.wdups.sorted.bam

	Qcheck3=\$( samtools quickcheck ${BAM}.wdups.sorted.bam )

	if [ `expr \${#Qcheck3}` -gt 0 ]
	then
	    echo "picard produced corrupt output file ${BAM}.wdups.sorted.bam \n\n"
	    exit 1
	fi    

    """

}

/*
* Step 3. Realignment
* Realignment with GATK
*/

process Realignment {
    tag                  { dedupBAM } 
    executor             params.executor
    cpus                 params.realignThreads
    queue                params.Queue
    memory               "${params.realignMem} GB"
    time                 params.realignTime
    module               params.gatkMod, params.samtoolsMod, params.tabixMod
    errorStrategy        'finish'
    scratch              '/scratch'
    stageOutMode         'copy'
    publishDir           realigndir, mode: 'link'

    input:
    file dedupBAM from dedupChannel
    file dedupBAI from dedupIdxChannel
    
    output:
    file("*realigned.bam") into realignChannel
    file("*realigned.bai") into realignIdxChannel    
    file "*realigned*"

    script:    
    fileName="${dedupBAM.baseName}"
    BAM="${fileName.replace('.wdups.sorted','')}"
    
    """
	echo `hostname`
        echo "RealignerTargetCreator step"
	echo `date`

	java ${params.JavaGATK} -jar ${params.gatkdir}/GenomeAnalysisTK.jar  \
	-T RealignerTargetCreator \
	-R ${params.refgenome} \
	-I ${dedupBAM} \
	-nt ${params.vcallThreads} \
	-o ${BAM}.realignTargetCreator.list

	if [ ! -s ${BAM}.realignTargetCreator.list ]
	then
	    echo "RealignerTargetCreator produced empty file ${BAM}.realignTargetCreator.list \n\n"
	    exit 1

	fi

        echo "IndelRealigner step"
	echo `date`
	
	java ${params.JavaGATK} -jar ${params.gatkdir}/GenomeAnalysisTK.jar  \
	-T IndelRealigner \
	-R ${params.refgenome} \
	-I ${dedupBAM} \
	-targetIntervals ${BAM}.realignTargetCreator.list \
	--bam_compression 0 \
	-o ${BAM}.realigned.bam

	Qcheck4=\$( samtools quickcheck ${BAM}.realigned.bam )

	if [ `expr \${#Qcheck4}` -gt 0 ]
	then
	    echo "IndelRealigner produced corrupt output file ${BAM}.realigned.bam \n\n"
	    exit 1
	fi
       
    """
}

/*
* Step 4. VarCallRaw
* Raw Variant discovery with GATK-HaplotypeCaller on a single sample
*/

process VarCallRaw {
    tag                  { realignBAM } 
    executor             params.executor
    cpus                 params.vcallThreads
    queue                params.vcallQueue
    memory               "${params.vcallMem} GB"
    time                 params.vcallTime
    module               params.gatkMod, params.tabixMod
    errorStrategy        'finish'
    scratch              '/scratch'
    stageOutMode         'copy'
    publishDir           vcalldir, mode: 'link'

    input:
    file realignBAM from realignChannel
    file realignBAI from realignIdxChannel

    
    output:
    file "${BAM}.rawVariants.g.vcf" into QCrawVcfChannel
    file "${BAM}.rawVariants.g.vcf.idx" into QCrawVcfIdxChannel
    file "*vcf*"
    
    script:    
    fileName="${realignBAM.baseName}"
    BAM="${fileName.replace('.realigned','')}"
    
    """
	echo `hostname`
    
	java ${params.JavaGATK} -jar ${params.gatkdir}/GenomeAnalysisTK.jar  \
	-T HaplotypeCaller \
	-R ${params.refgenome} \
	-I ${realignBAM} \
	-ERC GVCF  --variant_index_type LINEAR --variant_index_parameter 128000 \
	-gt_mode DISCOVERY \
	-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest \
	-A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
	-stand_call_conf 20 \
	-ploidy 2 \
	-o ${BAM}.rawVariants.g.vcf

	if [ ! -s ${BAM}.rawVariants.g.vcf ]
	then
	    echo "HaplotypeCaller produced empty file ${BAM}.rawVariants.g.vcf \n\n"
	    exit 1

	fi       
    """
}

/*
* Step 5. QCdedupBAM
* Run QC on dedup-BAM 
*/

process QCdedupBAM {
    tag                  { dedupBAM } 
    executor             params.executor
    cpus                 params.QCThreads
    queue                params.Queue
    memory               "${params.QCMem} GB"
    time                 params.QCTime
    module               params.picardMod, params.RMod, params.samtoolsMod
    errorStrategy        'ignore'
    publishDir           QCdir, mode: 'link'

    input:
    file dedupBAM from QCBamChannel
    file dedupBAI from QCBamIdxChannel
    file params.refgenome
    
    output:   
    file "${BAM}*"
    
    when:
    "${params.skipQC}" == 'FALSE'     
    
    script:    
    fileName="${dedupBAM.baseName}"
    BAM="${fileName.replace('.wdups.sorted','')}"
    
    """
    	echo `hostname`

	echo "Picard metrics for ${dedupBAM}"

	java ${params.JavaPicardQC} -jar ${params.picardir}/picard.jar \
        CollectMultipleMetrics \
        R=${params.refgenome} \
        I=${dedupBAM} \
        O=${BAM}_metrics \
        PROGRAM= null \
        PROGRAM= CollectAlignmentSummaryMetrics \
        PROGRAM= CollectInsertSizeMetrics \
        PROGRAM= CollectGcBiasMetrics

	echo "samtools metrics for ${dedupBAM}"
	samtools flagstat ${dedupBAM} > ${dedupBAM}_metrics.flagstat
	
    """
}


/*
* Step 6. QCVCF
* Run QC on rawVCF 
*/

process QCVCF {
    tag                  { rawVCF } 
    executor             params.executor
    cpus                 params.QCThreads
    queue                params.Queue
    memory               "${params.QCMem} GB"
    time                 params.QCTime
    module               params.tabixMod, params.vcftoolsMod
    errorStrategy        'ignore'
    publishDir           QCdir, mode: 'link'


    input:
    file rawVCF from QCrawVcfChannel
    file rawVCFidx from QCrawVcfIdxChannel
    
    output:    
    file "*.stats"
    
    when:
    "${params.skipQC}" == 'FALSE'    
    
    script:
    
    """
	    echo `hostname`	

	    echo "generate vcf stats for ${rawVCF}"
            bgzip -c ${rawVCF} > ${rawVCF}.gz	
            tabix -p vcf ${rawVCF}.gz
            vcf-stats ${rawVCF}.gz > ${rawVCF}.stats

	
    """
}

workflow.onComplete {
    def subject = "[Task #${params.ticket}]  HONEYBEE VARIANT CALLING PIPELINE with batchid ${params.batchid} and reads ${params.reads}"

    final_log = """

===============================================================================
                HONEYBEE VARIANT CALLING PIPELINE v${version}
 batchid ${params.batchid} and reads ${params.reads}
===============================================================================
Pipeline execution summary
===============================================================================
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
workDir      : ${workflow.workDir}
exit status  : ${workflow.exitStatus}
Error report : ${workflow.errorReport ?: '-'}
===============================================================================

    """
    
    ['mail', '-s', subject, params.email].execute() << "${final_log}\n\n${run_info}"
}

workflow.onError {
  // Display error message
  
  this.nextflowMessage()
  log.info "Workflow execution error message:"
  log.info "  " + workflow.errorMessage
  log.info "Detailed error report:"
  log.info "  " + workflow.errorReport
  
  // Send email too
  
  def subjectExit = "[Task #${params.ticket}]  HONEYBEE VARIANT CALLING PIPELINE with batchid ${params.batchid} and reads ${params.reads}"
  
  exit_log = """
  
===============================================================================
                HONEYBEE VARIANT CALLING PIPELINE v${version}
 batchid ${params.batchid} and reads ${params.reads}
===============================================================================  
Pipeline error report summary
===============================================================================  
exit status  : ${workflow.exitStatus}
Error message: ${workflow.errorMessage}
Error report : ${workflow.errorReport}
===============================================================================


   """
   
    ['mail', '-s', subjectExit, params.email].execute() << "${exit_log}\n\n${run_info}"   
  
}

