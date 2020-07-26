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
    stdout_data, stderr_data = p.communicate()
    p.wait()
    
    print("stdout".center(20,' ').center(40,'='))
    sys.stdout.flush()
    print(stdout_data)
    print("stderr".center(20,' ').center(40,'='))
    sys.stdout.flush()
    print(stderr_data)

    title2 = title + " end"
    print(time.asctime(time.localtime(time.time()))+" "+title2.center(30,' ').center(70,'*'))
    print(time.asctime(time.localtime(time.time()))+" "+name.center(30,' ').center(70,'*'))
    print("\n")
    sys.stdout.flush()

def GenomicsDBImport(inputFiles,outputFolder):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} GenomicsDBImport \
          --genomicsdb-workspace-path {1}".format(GATK,outputFolder)
    for i in inputFiles:
        cmd += " -V {0}".format(i)
    subprocessRun("GenomicsDBImport",outputFolder,cmd)

def CombineGVCFs(inputFiles,outputFile,ref):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} CombineGVCFs \
          -R {1} \
          -O {2} \
          -A BaseQuality \
          -A FisherStrand \
          -A MappingQuality \
          -A QualByDepth \
          -A StrandBiasBySample \
          -A StrandOddsRatio".format(GATK,ref,outputFile)
    for i in inputFiles:
        cmd += " --variant {0}".format(i)
    subprocessRun("CombineGVCFs",outputFile,cmd)

def GenotypeGVCFs(inputFile,outputFile,ref):
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
    cmd = "{0} GenotypeGVCFs \
          -R {1} \
          -O {2} \
          -V {3} \
          -A StrandBiasBySample".format(GATK,ref,outputFile,inputFile)
    subprocessRun("GenotypeGVCFs",inputFile,cmd)

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
        opts, args = getopt.getopt(argv,"hP:R:A:D:S:C:",["help","project=","reference=","annotation=","rawData=","suffix=","caller="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)

    project = ""
    reference = ""
    annotation = ""
    rawData = ""
    suffix = ""
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
    # jointCalling(project,reference,annotation,rawData,suffix,callers)

    inputFiles = glob.glob("{0}/*{1}".format(rawData,suffix))
    outputFile = "{0}/1_Files/7_Variants_GATK_GVCF/cohort.g.vcf".format(project)
    CombineGVCFs(inputFiles,outputFile,reference)

    inputFile = outputFile
    outputFile = "{0}/1_Files/7_Variants_GATK_GVCF/calling.vcf".format(project)
    GenotypeGVCFs(inputFile,outputFile,reference)

if __name__ == "__main__":
    main(sys.argv[1:])
