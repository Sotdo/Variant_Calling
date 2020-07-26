# import python modules
import re, glob, sys, os, subprocess, getopt, time

def mkdir(*folders):
    for folder in folders:
        subprocessRun('Make folder',folder ,'mkdir -p {0}'.format(folder))

def cp(file, folder):
    subprocessRun('Copy files',file +" > "+ folder,'cp {0} {1}'.format(file,folder))

def filePath2ID(file):
    searchObj = re.search('.*\/([^\n\t\r\f]+)_1\..*',file)
    ID = searchObj.group(1)
    return(ID)

def subprocessRun(title,name,cmd):

    title1 = title + " start"
    print(time.asctime(time.localtime(time.time()))+" "+title1.center(30,' ').center(70,'*'))
    print(time.asctime(time.localtime(time.time()))+" "+name.center(30,' ').center(70,'*'))
    sys.stdout.flush()

    p = subprocess.Popen([cmd],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()

    print("stdout".center(20,' ').center(40,'='))
    sys.stdout.flush()
    for i in p.stdout.readlines():
        print(str(i).strip('\n'))
        sys.stdout.flush()
    print("stderr".center(20,' ').center(40,'='))
    sys.stdout.flush()
    for i in p.stderr.readlines():
        print(str(i).strip('\n'))
        sys.stdout.flush()

    title2 = title + " end"
    print(time.asctime(time.localtime(time.time()))+" "+title2.center(30,' ').center(70,'*'))
    print(time.asctime(time.localtime(time.time()))+" "+name.center(30,' ').center(70,'*'))
    print("\n")
    sys.stdout.flush()

def fastp(in1,out1):
    in2 = re.sub('_1\.fastq','_2\.fastq',in1)
    out2 = re.sub('_1\.fastq','_2\.fastq',out1)
    html = re.sub('\.[^\/\n\t\r\f]+','.html',out1)
    json = re.sub('\.[^\/\n\t\r\f]+','.json',out1)

    fastpPath = "/home/suofang/Software/fastp-0.20.0/fastp"
    cmd = "{0} --length_required 70 \
               --in1 {1} \
               --in2 {2} \
               --out1 {3} \
               --out2 {4} \
               --html {5} \
               --json {6}".format(fastpPath,in1,in2,out1,out2,html,json)

    subprocessRun("Fastp",in1,cmd)

def BWA_MEM(ID,ref,in1,out):
    in2 = re.sub('_1\.fastq','_2\.fastq',in1)
    STR = re.sub('\.[^\/\n\t\r\f]+','',out)
    BWA = "/home/suofang/Software/bwa-0.7.17/bwa"
    cmd = "{0} mem \
               -R \'@RG\\tID:{1}\\tLB:lib1\\tPL:illumina\\tSM:{1}\\tPU:unit1\' \
               {2} \
               {3} \
               {4} \
               1> {5} \
               2> {6}".format(BWA,ID,ref,in1,in2,out,STR)
    subprocessRun("BWA MEM",in1,cmd)

def samtoolsSort(inputFile,outputFile):
    samtools = "/home/suofang/Software/samtools-1.9/samtools"
    cmd = "{0} sort -O bam \
                    -o {1} \
                    {2}".format(samtools,outputFile,inputFile)
    subprocessRun("samtools sort",inputFile,cmd)

def markduplicate(inputFile,outputFile):
    mark_dup = re.sub('\.[^\/\n\t\r\f]+','\.marked_dup_metrics.txt',outputFile)
    Picard = "java -jar /home/suofang/Software/picard_2.18.15.jar"
    cmd = "{0} MarkDuplicates I={1} \
                              O={2} \
                              M={3} \
                              CREATE_INDEX=true".format(Picard,inputFile,outputFile,mark_dup)
    subprocessRun("markduplicate",inputFile,cmd)
    
def IndexFeatureFile(inputFile):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} IndexFeatureFile \
           -I {1}".format(GATK,inputFile)
    subprocessRun("IndexFeatureFile",inputFile,cmd)
    
    
def BaseRecalibrator(inputFile,outputFile,ref,knownSite):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} BaseRecalibrator \
           -I {1} \
           -R {2} \
           --known-sites {3} \
           -O {4}".format(GATK,inputFile,ref,knownSite,outputFile)
    subprocessRun("BaseRecalibrator",inputFile,cmd)
    
def ApplyBQSR(inputFile,outputFile,ref,Base_Recal_table):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} ApplyBQSR \
           -I {1} \
           -R {2} \
           -bqsr {3} \
           --create-output-bam-index true \
           -O {4}".format(GATK,inputFile,ref,Base_Recal_table,outputFile)
    subprocessRun("ApplyBQSR",inputFile,cmd)
    
def ReorderSam(inputFile,outputFile,Ref):
#     DICT = re.sub('\.[^\.]+','\.dict',Ref)
    Picard = "java -jar /home/suofang/Software/picard_2.18.15.jar"
    cmd = "{0} ReorderSam \
          INPUT={1} \
          OUTPUT={2} \
          REFERENCE={3} \
          CREATE_INDEX=true".format(Picard,inputFile,outputFile,Ref)
    subprocessRun("ReorderSam",inputFile,cmd)

    
def bcftoolsCall(inputFile,outputFile,ref):
    bcftools = "/home/suofang/Software/bcftools-1.9/bcftools"
    cmd = "{0} mpileup -a FORMAT/DP4 -Ou -f {1} {2} |\
            {0} call -f GQ -mv -Ov -o {3}".format(bcftools,ref,inputFile,outputFile)
    subprocessRun("bcftools Call",inputFile,cmd)

def HaplotypeCaller(inputFile,outputFile,ref):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} HaplotypeCaller \
          -R {1} \
          -I {2} \
          -O {3} \
          --dont-use-soft-clipped-bases true \
          -A BaseQuality \
          -A FisherStrand \
          -A MappingQuality \
          -A QualByDepth \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,inputFile,outputFile)
    subprocessRun("HaplotypeCaller",inputFile,cmd)

def HaplotypeCaller_GVCF(inputFile,outputFile,ref):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} HaplotypeCaller \
          -R {1} \
          -I {2} \
          -O {3} \
          --dont-use-soft-clipped-bases true \
          -ERC GVCF \
          -A BaseQuality \
          -A FisherStrand \
          -A MappingQuality \
          -A QualByDepth \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,inputFile,outputFile)
    subprocessRun("HaplotypeCaller GVCF",inputFile,cmd)

def Deepvariant(inputFile,outputFolder,reference):

    file = re.sub('/.*/',"",inputFile)
    inputFolder = re.sub(file,"",inputFile)
    cp(reference,inputFolder)
    cp(reference+".fai",inputFolder)
    ref = re.sub("/.*/",'',reference)
    ID = re.sub("\..*","",file)

    BIN_VERSION = "0.10.0"

    cmd = 'docker run \
    -v "{0}":"/input" \
    -v "{1}:/output" \
    google/deepvariant:"{2}" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/{3} \
    --reads=/input/{4} \
    --output_vcf=/output/{5}.bwa.Deepvariant.vcf \
    --output_gvcf=/output/{5}.bwa.Deepvariant.g.vcf'.format(inputFolder,outputFolder,BIN_VERSION,ref,file,ID)
    subprocessRun("Deepvariant",inputFile,cmd)

def variantCalling(project,reference,annotation,rawData,suffix,callers,BQSR):

    for i in samples:
        print(i)

        sys.stdout.flush()

        # Step1 - Quality Control of Sequencing
        # Software: fastp
        BQSRFile = re.sub('\.[^\/\n\t\r\f]+','.Known.vcf',i)
        BQSRFile = re.sub(rawData,BQSR,BQSRFile)
        IndexFeatureFile(BQSRFile)

        # Step4 - Call Variants
        # Software: bcftools/GATK/Deepvariant

        inputFolder = rawData
        inputFile = i
        outputFolder = "{0}/1_Files/6_BaseRecalibration".format(project)
        mkdir(outputFolder)
        outputFile = re.sub('\.[^\/\n\t\r\f]+','.BaseRecalibrated.table',inputFile)
        outputFile = re.sub(inputFolder,outputFolder,outputFile)
        BaseRecalibrator(inputFile,outputFile,reference,BQSRFile)
        
        inputFile = outputFile
        inputFolder = outputFolder
        outputFolder = "{0}/1_Files/6_BaseRecalibration".format(project)
        mkdir(outputFolder)
        outputFile = re.sub('\.[^\/\n\t\r\f]+','.BaseRecalibrated.bam',inputFile)
        outputFile = re.sub(inputFolder,outputFolder,outputFile)
        ApplyBQSR(i,outputFile,reference,inputFile)
        
        inputFile = outputFile
        inputFolder = outputFolder
        outputFolder = "{0}/1_Files/6_BaseRecalibration".format(project)
        mkdir(outputFolder)
        outputFile = re.sub('\.[^\/\n\t\r\f]+','.BaseRecalibrated.Reordered.bam',inputFile)
        outputFile = re.sub(inputFolder,outputFolder,outputFile)
        ReorderSam(inputFile,outputFile,reference)
        

        for j in callers:

            if j == "Samtools":
                # bcftools
                outputFolder = "{0}/1_Files/7_Variants_Samtools".format(project)
                mkdir(outputFolder)
                outputFile = re.sub('\.[^\/\n\t\r\f]+','.bwa.Samtools.vcf',inputFile)
                outputFile = re.sub(inputFolder,outputFolder,outputFile)
                bcftoolsCall(inputFile,outputFile,reference)

            elif j == "GATK":
                #GATK
                outputFolder = "{0}/1_Files/7_Variants_GATK".format(project)
                mkdir(outputFolder)
                outputFile = re.sub('\.[^\/\n\t\r\f]+','.bwa.GATK.vcf',inputFile)
                outputFile = re.sub(inputFolder,outputFolder,outputFile)
                HaplotypeCaller(inputFile,outputFile,reference)

            elif j == "GATK_GVCF":
                #GATK_GVCF
                outputFolder = "{0}/1_Files/7_Variants_GATK_GVCF".format(project)
                mkdir(outputFolder)
                outputFile = re.sub('\.[^\/\n\t\r\f]+','.bwa.GATK.g.vcf',inputFile)
                outputFile = re.sub(inputFolder,outputFolder,outputFile)
                HaplotypeCaller_GVCF(inputFile,outputFile,reference)
                

            elif j == "Deepvariant":
                #Deepvariant
                outputFolder = "{0}/1_Files/7_Variants_Deepvariant".format(project)
                mkdir(outputFolder)
                Deepvariant(inputFile,outputFolder,reference)


def main(argv):

    usage = '''
    variantCalling.py

    Usage: python variantCalling.py [options]

    options:
        -h/--help
        -P/--project path to the project, like /data/yangyusheng/projectName
        -R/--reference path to the reference genome
        -A/--annotation path to the annotation
        -D/--data raw data, folder stores the raw data, like /data/yangyusheng/
        -S/--suffix raw data suffix, the suffix of *.fastq is .fastq
        -C/--caller Samtools, GATK, GATK_GVCF, Deepvariant
    '''

    try:
        opts, args = getopt.getopt(argv,"hP:R:A:D:S:C:B:",["help","project=","reference=","annotation=","rawData=","suffix=","caller=","BQSR="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)

    project = ""
    reference = ""
    annotation = ""
    rawData = ""
    suffix = ""
    BQSR = ""
    callers = []
    for opt, arg in opts:
        if opt == "-h" or opt == "--help":
            print(usage)
        elif opt == "-P" or opt == "--project":
            project = arg
        elif opt == "-R" or opt == "--reference":
            reference = arg
        elif opt == "-A" or opt == "--annotation":
            annotation = arg
        elif opt == "-D" or opt == "--rawData":
            rawData = arg
        elif opt == "-S" or opt == "--suffix":
            suffix = arg
        elif opt == "-C" or opt == "--caller":
            callers.append(arg)
        elif opt == "-B" or opt == "--BQSR":
            BQSR = arg
    variantCalling(project,reference,annotation,rawData,suffix,callers,BQSR)

if __name__ == "__main__":
    main(sys.argv[1:])

