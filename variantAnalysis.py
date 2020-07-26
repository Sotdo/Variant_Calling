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

def vcfallelicprimitives(inputVCF,outputVCF):
    cmd = "/home/suofang/Software/vcflib/bin/vcfallelicprimitives {0} > {1}".format(inputVCF,outputVCF)
    subprocessRun("vcfallelicprimitives",inputVCF,cmd)

def removeMultiAllelicSite_RefCall(inputVCF,outputVCF):
    file1 = open(inputVCF,'r')
    lines = list(file1.readlines())
    file1.close()
    file2 = open(outputVCF,'w')
    for i in lines:
        if i[0] == "#":
            file2.write(i)
        else:
            info = re.split('\t',i)
            if info[6] != "RefCall" and ',' not in info[4]:
                file2.write(i)
    file2.close()

def removeRefError(project,inputVCF,outputVCF):

    refErrorTxt = open('{0}/0_Preparation/refError.txt'.format(project),'r')
    refErrorLines = list(refErrorTxt.readlines())
    refErrorTxt.close()
    refError = []
    for i in refErrorLines:
        refError.append(i.strip('\n'))

    file1 = open(inputVCF,'r')
    lines = list(file1.readlines())
    file1.close()
    file2 = open(outputVCF,'w')
    for i in lines:
        if i[0] == "#":
            file2.write(i)
        else:
            info = re.split('\t',i)
            variant = info[0]+"\t"+info[1]+"\t"+info[3]+"\t"+info[4]
            if variant not in refError:
                file2.write(i)
    file2.close()


def variantNormalization(project,reference,annotation,rawData,suffix,callers):

    for i in callers:

        samples = glob.glob("{0}_{1}/*.{1}{2}".format(rawData,i,suffix))

        for j in samples:

            tmp1 = re.sub('\.[^\/\n\t\r\f]+','.Simple.vcf',j)
            removeMultiAllelicSite_RefCall(j,tmp1)

            tmp2 = re.sub('\.[^\/\n\t\r\f]+','.Normalized.vcf',tmp1)
            vcfallelicprimitives(tmp1,tmp2)

            outputVCF = re.sub('\.[^\/\n\t\r\f]+','.Normalized.NoRefError.vcf',tmp2)
            removeRefError(project,tmp2,outputVCF)






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
    variantNormalization(project,reference,annotation,rawData,suffix,callers)

if __name__ == "__main__":
    main(sys.argv[1:])
