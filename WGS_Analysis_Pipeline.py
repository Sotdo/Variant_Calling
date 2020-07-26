import re, glob, sys, os, subprocess, getopt, time, math
import workFlowFunctions

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
        -B/--BQSR 
    '''

    try:
        opts, args = getopt.getopt(argv,"hP:R:A:D:S:C:B",["help","project=","reference=","annotation=","rawData=","suffix=","caller=","BQSR"])
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
    reference, annotation, rawData  = workFlowFunctions.prepare(project,reference,annotation,rawData,suffix)
    
    workFlowFunctions.variantCalling(project,reference,annotation,rawData,suffix,callers)

if __name__ == "__main__":
    main(sys.argv[1:])