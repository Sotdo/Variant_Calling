import re, glob, sys, os, subprocess, getopt, time

def mkdir(*folders):
    for folder in folders:
        subprocessRun('Make folder',folder ,'mkdir -p {0}'.format(folder))

def cp(file, folder):
    subprocessRun('Copy files',file +" > "+ folder,'cp {0} {1}'.format(file,folder))

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


def prepare(project,reference,annotation,rawData,suffix):

    print("args")
    print("#"*70)
    print("Project: " + project)
    print("Reference: " + reference)
    print("Annotation: " + annotation)
    print("rawData: " + rawData)
    print("suffix: " + suffix)
    print("\n")

    fastp = "/home/suofang/Software/fastp-0.20.0/fastp"
    BWA = "/home/suofang/Software/bwa-0.7.17/bwa"
    samtools = "/home/suofang/Software/samtools-1.9/samtools"
    bcftools = "/home/suofang/Software/bcftools-1.9/bcftools"
    Picard = "java -jar /home/suofang/Software/picard_2.18.15.jar"
    GATK = "java -jar /home/suofang/Software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"

    print("Software used in this pipeline")
    print("#"*70)
    print(fastp)
    print(BWA)
    print(samtools)
    print(bcftools)
    print(Picard)
    print(GATK)
    print("\n\n")
    print("Procedure")
    print("#"*70)
    sys.stdout.flush()

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
    subprocessRun('BWA Index',reference,"{0} index {1}".format(BWA,reference))
    subprocessRun('Samtools faidx',reference,"{0} faidx {1}".format(samtools,reference))
    subprocessRun('Picard CreateSequenceDictionary', reference+" > "+referenceDict, "{0} CreateSequenceDictionary R={1} O={2}".format(Picard,reference,referenceDict))



def main(argv):

    usage = '''
    preparation.py

    Usage: python preparation.py [options]

    options:
        -h/--help
        -P/--project path to the project, like /data/yangyusheng/projectName
        -R/--reference path to the reference genome
        -A/--annotation path to the annotation
        -D/--data raw data, folder stores the raw data, like /data/yangyusheng/
        -S/--suffix raw data suffix, the suffix of *.fastq is .fastq
    '''

    try:
        opts, args = getopt.getopt(argv,"hP:R:A:D:S:",["help","project=","reference=","annotation=","rawData=","suffix="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)

    print(opts)
    project = ""
    reference = ""
    annotation = ""
    rawData = ""
    suffix = ""
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

    prepare(project,reference,annotation,rawData,suffix)

if __name__ == "__main__":
    main(sys.argv[1:])
