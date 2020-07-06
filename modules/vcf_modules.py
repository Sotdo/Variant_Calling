import os, re, glob, subprocess

class vcf:

    file = ''
    header = []
    title = ''
    sites = []
    variants = []
    annotations = {}
    ID = ''

    def __init__(self,filePath):

        self.file = ''
        self.header = []
        self.title = ''
        self.sites = []
        self.variants = []
        self.annotations = {}
        self.ID = ''

        self.file = filePath

        inputs = open(self.file,'r')
        lines = list(inputs.readlines())
        inputs.close()

        searchObj = re.search('.*\/([^\.\n\t\r\f]+)\..*',self.file)
        self.ID = searchObj.group(1)

        for i in lines:

            i = i.strip('\n')
            if i[0] == "#" and i[1] == "#":

                self.header.append(i)

            elif i[0] == "#":

                self.title = i

            else:

                cells = re.split('\t',i)

                Chr = cells[0]
                Posi = cells[1]
                Ref = cells[3]
                Alt = cells[4]
                Qual = cells[5]
                Type = ''

                if len(Ref) == len(Alt) and len(Ref) == 1:
                    Type = 'SNP'
                elif len(Ref) > len(Alt):
                    Type = "Deletion"
                elif len(Ref) < len(Alt):
                    Type = "Insertion"

                site = Chr + "\t" + Posi
                self.sites.append(site)

                variant = Chr + "\t" + Posi + "\t" + Ref + "\t" + Alt
                self.variants.append(variant)
                self.annotations[variant] = {}
                self.annotations[variant]['Quality'] = Qual
                self.annotations[variant]['Type'] = Type

                info = re.split(';',cells[7])

                for j in info:

                    KeyValue = re.search('(\S+)=(\S+)',j)
                    if KeyValue != None:
                        Key = KeyValue.group(1)
                        Value = KeyValue.group(2)
                        self.annotations[variant][Key] = Value

                Formats = re.split(':',cells[8])
                FormatsValue = re.split(':',cells[9])

                for j in range(len(Formats)):
                    
                    Key = Formats[j]
                    Value = FormatsValue[j]

                    if Key == "AD":
                        AD = re.split(',',Value)
                        Value = (float(AD[0])+1)/(float(AD[1])+1)

                    self.annotations[variant][Key] = Value
