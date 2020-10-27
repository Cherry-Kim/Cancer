import string,sys,glob

CMS={}
fp=open("res_prediction.txt","r")
fp.readline()
for line in fp:
	line_temp=line.split('\t')
	if not 'NA' == line_temp[1]:
		if not line_temp[1] in CMS:
			CMS[line_temp[1]]=[]
		CMS[line_temp[1]].append(line_temp[0])

for key,value in sorted(CMS.items()):
#	print key,":",len(value)	#CMS4 [sample1, sample2]
	fp1=open('Normalized.t.txt','r')
	hd=fp1.readline()
	g=hd[:-1].split('\t')
	fpout=[open(g[i+1]+'_cms.txt','a')for i in range(0,20)]
	for line in fp1:
		line_temp1=line[:-1].split('\t')	#sample1.genes.results
		for j in value:
			if j in line_temp1[0]:
				for i in range(20):
					if ('CMS1' == key) and ('sample1'==j):
						fpout[i].write("CMS"+'\t'+"TPM"+'\t'+"GENE"+'\n')
						fpout[i].write(str(key)+'\t'+line_temp1[i+1]+'\t'+g[i+1]+'\n')
					else:
						fpout[i].write(str(key)+'\t'+line_temp1[i+1]+'\t'+g[i+1]+'\n')
