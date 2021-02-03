#!/usr/bin/env Rscript

Survival_rna <- function (
	GENE
	) {
	library(cgdsr) 
	mycgds = CGDS("http://www.cbioportal.org/")
	cancerstudy <- getCancerStudies(mycgds)
	#cancerstudy$name	#[72] "Colorectal Adenocarcinoma (TCGA, PanCancer Atlas)"
	mycancerstudy = cancerstudy[72,1]

	getCaseLists(mycgds,mycancerstudy)[,1]	#[8] "coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
	mycaselist = getCaseLists(mycgds,mycancerstudy)[8,1]

	
	getGeneticProfiles(mycgds,mycancerstudy)[,1]	#[4] "coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
	myrnaprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]
	AST_mirna = getProfileData(mycgds, genes=GENE, myrnaprofile, mycaselist)
	BCL11A_mirna_cases = rownames(AST_mirna[which(AST_mirna$BCL11A>1.5), ,drop=F])
	##### 4. Clinical data integration
	myclinicaldata = getClinicalData(mycgds,mycaselist)
	myclinicaldata$OS_STATUS[myclinicaldata$OS_STATUS == ""] <- NA	#"" -> NA  #"0:LIVING"   "1:DECEASED"
	total_sample <- rownames(myclinicaldata) 
	type <- rep('Wild', length(total_sample))
	names(type) <- total_sample
	type[BCL11A_mirna_cases] <- "BCL11A_mirna"
	type <- type[names(type) %in% rownames(myclinicaldata)]
	type <- factor(type, levels= c("Wild","BCL11A_mirna"))
	# Surival curves were fitted using the Kaplan-Meier formular in the R package 'survival'.
	# and visualized using the ggsurvplot function of the R package 'survminer'
	library(survival)
	library(survminer)
	BCL11A_fit <- survfit(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED")~type, data=myclinicaldata)
        # fit < survfit(Surv(time, status) ~ sex, data=data)
	png("b.png")
	ggsurvplot(BCL11A_fit, data=myclinicaldata, pval=T, risk.table=T, conf.int=F, break.time.by = 50, legend.title = "Patient types", risk.table.fontsize = 2.5,surv.plot.height = 0.45, legend="right")
	dev.off()
}
Survival_rna( GENE <- c('BCL11A'))
