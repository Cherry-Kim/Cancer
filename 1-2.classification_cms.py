import string,sys,glob

CMS={}
#GENE=['ENST00000307564_AKNA','ENST00000356842_BCL11A','ENST00000290551_BTG2','ENST00000376670_GATA1','ENST00000339463_GFI1B','ENST00000595883_SPIB','ENST00000331340_IKZF1','ENST00000357234_IRF5','ENST00000268638_IRF8','ENST00000342854_TCF7','ENST00000302754_JUNB','ENST00000294339_TAL1','ENST00000526816_POU2F2','ENST00000346872_IKZF3','ENST00000248071_KLF2','ENST00000341911_MYB','ENST00000273668_EAF2','ENST00000264834_KLF1','ENST00000435572_NFE2','ENST00000378538_SPI1']

fp=open("res_prediction.txt","r")
fp.readline()
for line in fp:
	line_temp=line.split('\t')
	if not 'NA' == line_temp[1]:
		if not line_temp[1] in CMS:
			CMS[line_temp[1]]=[]
		CMS[line_temp[1]].append(line_temp[0])

for key,value in sorted(CMS.items()):
#	print key,":",len(value)	#CMS4 [PM.AU.0035.T, PM.AU.0036.T]
	fp1=open('Normalized.t.txt','r')
	hd=fp1.readline()
	g=hd[:-1].split('\t')
	fpout=[open(g[i+1]+'_cms3.txt','a')for i in range(0,20)]
	for line in fp1:
		line_temp1=line[:-1].split('\t')	#PM.AU.0035.T.genes.results
		for j in value:
			if j in line_temp1[0]:
				for i in range(20):
					if ('CMS1' == key) and ('PM.AS.0035.T'==j):
						fpout[i].write("CMS"+'\t'+"TPM"+'\t'+"GENE"+'\n')
						fpout[i].write(str(key)+'\t'+line_temp1[i+1]+'\t'+g[i+1]+'\n')
					else:
						fpout[i].write(str(key)+'\t'+line_temp1[i+1]+'\t'+g[i+1]+'\n')
