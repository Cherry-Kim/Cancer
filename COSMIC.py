import requests
import string,sys,os,glob

GENE=['BRACA1','BCL11A_ENST00000356842','IKZF1_ENST00000331340','BRACA2']

def DOWNLOAD(omics):
	print omics
	email    = "k@snu.ac"
	password = "!"
	url      = "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v90/"

	r = requests.get(url+omics, auth=(email, password))
	download_url = r.json()["url"]
	r = requests.get(download_url)

	with open(omics, "wb") as f:
	    f.write(r.content)


def Point_mutation(Pname):
	PGENE=['BRACA1_ENST00000307564','BCL11A_ENST00000356842','BRACA2_ENST00000331340']
	fp=open(Pname,'r')
	hd=fp.readline()
	fpout=open('Point_mutaion.txt','w')
	hd_temp=['Gene name','Transcript','Sample Name','SampleID','AA Mutation','CDS Mutation','Mutation Description','Primary Tissue','Tissue Subtype1','Tissue Subtype2','Histology','Histology Subtype1','Histology Subtype2','PubmedID','Somatic Status','Sample Type','Zygosity','GRCh','Mutation Genome Position','Age']
	fpout.write('\t'.join(hd_temp)+'\n')

	for line in fp:
       		line_temp=line.split('\t')
        	if not ('NS' == line_temp[7]): #Primary Tissue
               		if not ('Unknown' == line_temp[21]):    #Mutation
               	        	for i in PGENE:
                               		trans=i.split('_')[1]
                               		if trans in line_temp[1]:
	                                        temp=[line_temp[0],line_temp[1],line_temp[4],line_temp[5],line_temp[20],line_temp[19],line_temp[21],line_temp[7],line_temp[8],line_temp[9],line_temp[11],line_temp[12],line_temp[13],line_temp[32],line_temp[31],line_temp[34],line_temp[22],line_temp[24],line_temp[25],line_temp[36]]
        	                                fpout.write('\t'.join(temp)+'\n')
	fp.close()
	fpout.close()

def Methylation(Mname):
	fp=open(Mname,'r')
	hd=fp.readline()
	fpout=open('Methylation.txt','w')
	hd_temp=['Gene name','Sample name','SampleID','ProbeID','GRCh','Chromosome','Position','Methylation Type','Level(Beta-value)','Normal Average','Primary Tissue','Tissue Subtype1','Histology','Histology Subtype1','Study','P_VALUE']
	fpout.write('\t'.join(hd_temp)+'\n')

	for line in fp:
	        line_temp=line.split('\t')
	##      if ('1'==line_temp[16]) and breast'==line_temp[4]: #positive
	        for i in GENE:
        	        if i == line_temp[17]:
                	        temp=[line_temp[17],line_temp[2],line_temp[1],line_temp[12],line_temp[13],line_temp[14],line_temp[15],line_temp[18],line_temp[20],line_temp[19],line_temp[4],line_temp[5],line_temp[8],line_temp[9],line_temp[0],line_temp[21]]
                       		fpout.write('\t'.join(temp))
	fp.close()
	fpout.close()


def CNV(Cname):
	fp=open(Cname,'r')
	hd=fp.readline()
	fpout=open('CosmicCNV.txt','w')
	hd_temp=['SampleID','Gene','Primary site','CN Type','Minor Allele','Copy Number','GRCh','CNV segment Pos.','Study']
	fpout.write('\t'.join(hd_temp)+'\n')

	for line in fp:
	        line_temp=line[:-1].split('\t')
	#       if 'breast' in line_temp[5]:
	        for i in GENE:
	                if i == line_temp[2]:
	                        temp=[line_temp[3],line_temp[2],line_temp[5],line_temp[16],line_temp[15],line_temp[14],line_temp[18],line_temp[19],line_temp[17]]
	                        fpout.write('\t'.join(temp)+'\n')
	fp.close()
	fpout.close()

def Expression(Ename):
	fp=open(Ename,'r')
	fp.readline()
	fpout=open('expression.txt','w')
	for line in fp:
	        line_temp=line.split('\t') 
       		for i in GENE:
               		if i==line_temp[2]:
                       		fpout.write(line)
	fp.close()
	fpout.close()

def Main():
	Files = ['CosmicMutantExport.tsv.gz','CosmicCompleteDifferentialMethylation.tsv.gz','CosmicCompleteCNA.tsv.gz','CosmicCompleteGeneExpression.tsv.gz']
	for omics in Files:
		DOWNLOAD(omics)
	Pname = Files[0].split('.gz')[0]
	Point_mutation(Pname)
	Mname = Files[1].split('.gz')[0]
	Methylation(Mname)
	Cname = Files[2].split('.gz')[0]
	CNV(Cname)
	Ename = Files[3].split('.gz')[0]
	Expression(Ename)
Main()

