import re, glob, sys, os, subprocess, getopt, time, math

#Software
Fastp = "/home/suofang/Software/fastp-0.20.0/fastp"
BWA = "/home/suofang/Software/bwa-0.7.17/bwa"
Samtools = "/home/suofang/Software/samtools-1.9/samtools"
Bcftools = "/home/suofang/Software/bcftools-1.9/bcftools"
Picard = "java -jar /home/suofang/Software/picard_2.18.15.jar"
GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"

def subprocessRun(title,name,cmd):

    title1 = title + " start"
    print(time.asctime(time.localtime(time.time()))+" "+title1.center(30,' ').center(70,'*'))
    print(time.asctime(time.localtime(time.time()))+" "+name.center(30,' ').center(70,'*'))
    sys.stdout.flush()

    p = subprocess.Popen([cmd],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout_data, stderr_data = p.communicate()
    p.wait()
    
    print("stdout".center(20,' ').center(40,'='))
    print(stdout_data)
    print("stderr".center(20,' ').center(40,'='))
    print(stderr_data)
    sys.stdout.flush()

    title2 = title + " end"
    print(time.asctime(time.localtime(time.time()))+" "+title2.center(30,' ').center(70,'*'))
    print(time.asctime(time.localtime(time.time()))+" "+name.center(30,' ').center(70,'*'))
    print("\n")
    sys.stdout.flush()
    
def mkdir(*folders):
    for folder in folders:
        subprocessRun('Make folder',folder ,'mkdir -p {0}'.format(folder))
        
def cp(file, folder):
    subprocessRun('Copy files',file +" > "+ folder,'cp {0} {1}'.format(file,folder))
    
def info(project,reference,annotation,rawData,suffix):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK

    print("args")
    print("#"*70)
    print("Project: " + project)
    print("Reference: " + reference)
    print("Annotation: " + annotation)
    print("rawData: " + rawData)
    print("suffix: " + suffix)
    print("\n")

    print("Software used in this pipeline")
    print("#"*70)
    print(Fastp)
    print(BWA)
    print(Samtools)
    print(Bcftools)
    print(Picard)
    print(GATK)
    print("\n\n")
    print("Procedure")
    print("#"*70)
    sys.stdout.flush()
    
def pairEndFastq2ID(file):
    searchObj = re.search('.*\/([^\/\.\n\t\r\f]+)_1\..*',file)
    ID = searchObj.group(1)
    return(ID)
    
def filePath2ID(file):
    searchObj = re.search('.*\/([^\/\.\n\t\r\f]+)\..*',file)
    ID = searchObj.group(1)
    return(ID)

def prepare(project,reference,annotation,rawData,suffix):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    Preparation_0 = "{0}/0_Preparation".format(project)
    Files_1 = "{0}/1_Files".format(project)
    Summary_2 = "{0}/2_Summary".format(project)
    scripts = "{0}/scripts".format(project)
    dataFolder = "{0}/0_RawData".format(Files_1)
    mkdir(project, Preparation_0, Files_1, Summary_2, scripts, dataFolder)
    cp(reference,Preparation_0)
    cp(annotation,Preparation_0)
    cp(rawData+"/*"+suffix,dataFolder)

    reference = re.sub('/.*/',Preparation_0+"/",reference)
    annotation = re.sub('/.*/',Preparation_0+"/",annotation)

    referenceDict = re.sub('\..*','.dict',reference)
    
    #BWA Index
    cmd = "{0} index {1}".format(BWA,reference)
    subprocessRun('BWA Index',reference,cmd)
    
    #Samtools faidx
    cmd = "{0} faidx {1}".format(Samtools,reference)
    subprocessRun('Samtools faidx',reference, cmd)
    
    #Picard CreateSequenceDictionary
    cmd = "{0} CreateSequenceDictionary R={1} O={2}".format(Picard,reference,referenceDict)
    subprocessRun('Picard CreateSequenceDictionary', reference+" > "+referenceDict, cmd)
    return(reference,annotation,dataFolder)
    
def fastp(in1,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    in2 = re.sub('_1\.fastq','_2.fastq',in1)
    out1 = re.sub('(.*\/)',outputFolder+'/',in1)
    out1 = re.sub('_1\.[^\.\/\n\t\r\f]+','_1.fastp.fastq',out1)
    out2 = re.sub('_1\.fastp\.fastq','_2.fastp.fastq',out1)
    html = re.sub('_1\.[^\/\n\t\r\f]+','.html',out1)
    json = re.sub('_1\.[^\/\n\t\r\f]+','.json',out1)
    
    cmd = "{0} --length_required 70 \
               --in1 {1} \
               --in2 {2} \
               --out1 {3} \
               --out2 {4} \
               --html {5} \
               --json {6}".format(Fastp,in1,in2,out1,out2,html,json)
    subprocessRun("Fastp",in1,cmd)
    return(out1)
    
def BWA_MEM(ref,in1,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    ID = pairEndFastq2ID(in1)
    in2 = re.sub('_1\.fastp\.fastq','_2.fastp.fastq',in1)
    out = re.sub('_1\.[^\/\n\t\r\f]+','.bwa.sam',in1)
    out = re.sub('(.*\/)',outputFolder+'/',out)
    STR = re.sub('\.[^\/\n\t\r\f]+','',out)

    cmd = "{0} mem \
               -R \'@RG\\tID:{1}\\tLB:lib1\\tPL:illumina\\tSM:{1}\\tPU:unit1\' \
               {2} \
               {3} \
               {4} \
               1> {5} \
               2> {6}".format(BWA,ID,ref,in1,in2,out,STR)
    subprocessRun("BWA MEM",in1,cmd)
    return(out)

def samtoolsSort(inputFile,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.bwa.srt.bam',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    
    cmd = "{0} sort -O bam \
                    -o {1} \
                    {2}".format(Samtools,outputFile,inputFile)
    subprocessRun("samtools sort",inputFile,cmd)
    return(outputFile)
    
def markduplicate(inputFile,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.bwa.srt.mark_dup.bam',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    mark_dup = re.sub('\.[^\/\n\t\r\f]+','.marked_dup_metrics.txt',outputFile)
    
    cmd = "{0} MarkDuplicates I={1} \
                              O={2} \
                              M={3} \
                              CREATE_INDEX=true".format(Picard,inputFile,outputFile,mark_dup)
    subprocessRun("markduplicate",inputFile,cmd)
    return(outputFile)

def bcftoolsCall(ref,inputFile,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.Samtools.vcf',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    cmd = "{0} mpileup -a FORMAT/DP4 -Ou -f {1} {2} |\
            {0} call -f GQ -mv -Ov -o {3}".format(Bcftools,ref,inputFile,outputFile)
    subprocessRun("bcftools Call",inputFile,cmd)
    return(outputFile)

def HaplotypeCaller(ref,inputFile,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.GATK.vcf',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    cmd = "{0} HaplotypeCaller \
          -R {1} \
          -I {2} \
          -O {3} \
          --dont-use-soft-clipped-bases true \
          -A AS_BaseQualityRankSumTest \
          -A AS_FisherStrand \
          -A AS_MappingQualityRankSumTest \
          -A AS_QualByDepth \
          -A AS_ReadPosRankSumTest \
          -A AS_RMSMappingQuality \
          -A AS_StrandOddsRatio \
          -A Coverage \
          -A DepthPerAlleleBySample \
          -A DepthPerSampleHC \
          -A FisherStrand \
          -A MappingQuality \
          -A MappingQualityRankSumTest \
          -A QualByDepth \
          -A ReadPosRankSumTest \
          -A RMSMappingQuality \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,inputFile,outputFile)
    
    subprocessRun("HaplotypeCaller",inputFile,cmd)
    return(outputFile)
  
def HaplotypeCaller_GVCF(ref,inputFile,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.GATK_GVCF.g.vcf',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    cmd = "{0} HaplotypeCaller \
          -R {1} \
          -I {2} \
          -O {3} \
          --dont-use-soft-clipped-bases true \
          -ERC GVCF \
          -A AS_BaseQualityRankSumTest \
          -A AS_FisherStrand \
          -A AS_MappingQualityRankSumTest \
          -A AS_QualByDepth \
          -A AS_ReadPosRankSumTest \
          -A AS_RMSMappingQuality \
          -A AS_StrandOddsRatio \
          -A Coverage \
          -A DepthPerAlleleBySample \
          -A DepthPerSampleHC \
          -A FisherStrand \
          -A MappingQuality \
          -A MappingQualityRankSumTest \
          -A QualByDepth \
          -A ReadPosRankSumTest \
          -A RMSMappingQuality \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,inputFile,outputFile)
    subprocessRun("HaplotypeCaller GVCF",inputFile,cmd)
    return(outputFile)

def Deepvariant(ref,inputFile,outputFolder):

    file = re.sub('(.*\/)',"",inputFile)
    inputFolder = re.sub(file,"",inputFile)
    
    #Deepvariant needs the reference and the reference.fai file
    cp(ref,inputFolder)
    cp(ref+".fai",inputFolder)
    
    ref = re.sub("(.*\/)",'',ref)
    ID = filePath2ID(inputFile)
    outputFile = outputFolder + "/" + ID + ".Deepvariant.vcf"

    BIN_VERSION = "0.10.0"

    cmd = 'docker run \
    -v "{0}":"/input" \
    -v "{1}:/output" \
    google/deepvariant:"{2}" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/{3} \
    --reads=/input/{4} \
    --output_vcf=/output/{5}.Deepvariant.vcf \
    --output_gvcf=/output/{5}.Deepvariant.g.vcf'.format(inputFolder,outputFolder,BIN_VERSION,ref,file,ID)
    subprocessRun("Deepvariant",inputFile,cmd)
    return(outputFile)

def IndexFeatureFile(inputFile):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    cmd = "{0} IndexFeatureFile \
           -I {1}".format(GATK,inputFile)
    subprocessRun("IndexFeatureFile",inputFile,cmd)
    
def BaseRecalibrator(ref,inputFile,knownSite,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.BaseRecalibrated.table',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    
    cmd = "{0} BaseRecalibrator \
           -I {1} \
           -R {2} \
           --known-sites {3} \
           -O {4}".format(GATK,inputFile,ref,knownSite,outputFile)
    subprocessRun("BaseRecalibrator",inputFile,cmd)
    return(outputFile)
    
def ApplyBQSR(ref,inputFile,outputFolder,Base_Recal_table):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
        
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.BaseRecalibrated.bam',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    
    cmd = "{0} ApplyBQSR \
           -I {1} \
           -R {2} \
           -bqsr {3} \
           --create-output-bam-index true \
           -O {4}".format(GATK,inputFile,ref,Base_Recal_table,outputFile)
    subprocessRun("ApplyBQSR",inputFile,cmd)
    return(outputFile)
    
def ReorderSam(Ref,inputFile,outputFolder):

    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.BaseRecalibrated.Reordered.bam',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    
    cmd = "{0} ReorderSam \
          INPUT={1} \
          OUTPUT={2} \
          REFERENCE={3} \
          CREATE_INDEX=true".format(Picard,inputFile,outputFile,Ref)
    subprocessRun("ReorderSam",inputFile,cmd) 
    return(outputFile)
    
    
def CombineGVCFs(ref,inputFiles,outputFolder):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = outputFolder + "/cohort.g.vcf"
    
    cmd = "{0} CombineGVCFs \
          -R {1} \
          -O {2} \
          -A AS_BaseQualityRankSumTest \
          -A AS_FisherStrand \
          -A AS_MappingQualityRankSumTest \
          -A AS_QualByDepth \
          -A AS_ReadPosRankSumTest \
          -A AS_RMSMappingQuality \
          -A AS_StrandOddsRatio \
          -A Coverage \
          -A DepthPerAlleleBySample \
          -A DepthPerSampleHC \
          -A FisherStrand \
          -A MappingQuality \
          -A MappingQualityRankSumTest \
          -A QualByDepth \
          -A ReadPosRankSumTest \
          -A RMSMappingQuality \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,outputFile)
    for i in inputFiles:
        cmd += " --variant {0}".format(i)
    subprocessRun("CombineGVCFs",outputFile,cmd)
    return(outputFile)

def GenotypeGVCFs(inputFile,outputFolder,ref):
    
    global Fastp,BWA,Samtools,Bcftools,Picard,GATK
    
    outputFile = re.sub('\.[^\/\n\t\r\f]+','.GATK_GVCF.vcf',inputFile)
    outputFile = re.sub('(.*\/)',outputFolder+'/',outputFile)
    
    cmd = "{0} GenotypeGVCFs \
          -R {1} \
          -O {2} \
          -V {3} \
          -A AS_BaseQualityRankSumTest \
          -A AS_FisherStrand \
          -A AS_MappingQualityRankSumTest \
          -A AS_QualByDepth \
          -A AS_ReadPosRankSumTest \
          -A AS_RMSMappingQuality \
          -A AS_StrandOddsRatio \
          -A Coverage \
          -A DepthPerAlleleBySample \
          -A DepthPerSampleHC \
          -A FisherStrand \
          -A MappingQuality \
          -A MappingQualityRankSumTest \
          -A QualByDepth \
          -A ReadPosRankSumTest \
          -A RMSMappingQuality \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,outputFile,inputFile)
    subprocessRun("GenotypeGVCFs",inputFile,cmd)
    return(outputFile)
    

def variantCalling(project,reference,annotation,rawData,suffix,callers):

    samples = glob.glob('{0}/*{1}'.format(rawData,suffix))

    for i in samples:

        # Step1 - Quality Control of Sequencing
        # Software: fastp

        inputFolder = rawData
        outputFolder = "{0}/1_Files/1_Fastp".format(project)
        mkdir(outputFolder)
        in1 = i
        output = fastp(in1,outputFolder)

        # Step2 - Mapping
        # Software: BWA-MEM
        inputFolder = outputFolder
        outputFolder = "{0}/1_Files/2_BWA".format(project)
        mkdir(outputFolder)
        in1 = output
        output = BWA_MEM(reference,in1,outputFolder)

        # Step3 - Process SAM
        # Software: samtools + Picard
        # Sort
        inputFolder = outputFolder
        outputFolder = "{0}/1_Files/3_ProcessedBAM".format(project)
        mkdir(outputFolder)
        inputFile = output
        output = samtoolsSort(inputFile,outputFolder)
        # Mark duplicate
        inputFile = output
        output = markduplicate(inputFile,outputFolder)

        # Step4 - Call Variants
        # Software: bcftools/GATK/Deepvariant

        inputFolder = outputFolder
        inputFile = output

        for j in callers:

            if j == "Samtools":
                # bcftools
                outputFolder = "{0}/1_Files/4_Variants_Samtools".format(project)
                mkdir(outputFolder)
                output = bcftoolsCall(reference,inputFile,outputFolder)

            elif j == "GATK":
                #GATK
                outputFolder = "{0}/1_Files/4_Variants_GATK".format(project)
                mkdir(outputFolder)
                output = HaplotypeCaller(reference,inputFile,outputFolder)

            elif j == "GATK_GVCF":
                #GATK_GVCF
                outputFolder = "{0}/1_Files/4_Variants_GATK_GVCF".format(project)
                mkdir(outputFolder)
                output = HaplotypeCaller_GVCF(reference,inputFile,outputFolder)

            elif j == "Deepvariant":
                #Deepvariant
                outputFolder = "{0}/1_Files/4_Variants_Deepvariant".format(project)
                mkdir(outputFolder)
                output = Deepvariant(reference,inputFile,outputFolder)