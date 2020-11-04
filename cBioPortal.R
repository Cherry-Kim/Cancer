##########   CASE I  ##########
#install.packages('cgdsr') 
library(cgdsr) 

##### 1. Get list of cancer studies at server
mycgds = CGDS("http://www.cbioportal.org/")
cancerstudy <- getCancerStudies(mycgds)
head(cancerstudy)
cancerstudy$name
mycancerstudy = cancerstudy[24,1]

##### 2. Extract samples and features
getCaseLists(mycgds,mycancerstudy)[,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[2,1]

getGeneticProfiles(mycgds,mycancerstudy)
mymutationprofile = getGeneticProfiles(mycgds,mycancerstudy)[9,1] 
mymrnaprofile=getGeneticProfiles(mycgds,mycancerstudy)[8,1] 
#mymethylationprofile=getGeneticProfiles(mycgds,mycancerstudy)[5,1] 

##### 3-1. Get mutation profiles for genes
AST_mutation = getMutationData(mycgds,mycaselist,mymutationprofile,c('BRACA1','BRACA2') ) 
head(AST_mutation)
table(AST_mutation$gene_symbol) 

BRACA1_mutated_cases2 <- AST_mutation[which(AST_mutation$gene_symbol=='BRACA1'),3]	#3: case_id
BRACA1_mutated_cases <- paste0("X",AKNA_mutated_cases2)
head(BRACA1_mutated_cases)
#brca2_mutated_cases <- brca_mutation[which(brca_mutation$gene_symbol=='BRCA2'),3]

##### 3-2. Extract samples with BRCA1, BRCA2 methylation
AST_mirna = getProfileData(mycgds, c('BRACA1','BRACA2'),myrnaprofile, mycaselist)
#brca_methylation = getProfileData(mycgds, c('BRCA1','BRCA2'),mymethylationprofile, mycaselist)
head(AST_mirna)
AKNA_mirna_cases = rownames(AST_mirna[which(AST_mirna$AKNA>1.0 | AST_mirna$AKNA < -1.0),])
#brca1_methylation_cases = rownames(brca_methylation[which(brca_methylation$BRCA1>0.8),])
#brca2_methylation_cases = rownames(brca_methylation[which(brca_methylation$BRCA2>0.8),])

##### 4. Clinical data integration
myclinicaldata = getClinicalData(mycgds,mycaselist)
head(myclinicaldata)

myclinicaldata$OS_STATUS[myclinicaldata$OS_STATUS == ""] <- NA	#"" -> NA

total_sample <- rownames(myclinicaldata) 
type <- rep('Wild', length(total_sample))
names(type) <- total_sample
head(type)

AKNA_mutated_cases = gsub("-",".",AKNA_mutated_cases) #- -> .
AKNA_mirna_cases = gsub("-",".",AKNA_mutated_cases)

type[AKNA_mutated_cases] <- "BRACA1_mutation"
type[AKNA_mirna_cases] <- "BRACA1_mirna"

type <- type[names(type) %in% rownames(myclinicaldata)]
type <- factor(type, levels= c("Wild","BRACA1_mutation","BRACA1_mirna") )

##### 6. Survival Analysis
#install.packages('survival')
library(survival)

out <- survfit(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED")~type, data=myclinicaldata)

survdiff(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED") ~ type,data=myclinicaldata)
coxph(Surv(OS_MONTHS, OS_STATUS=="1:DECEASED") ~ type, data=myclinicaldata)

color <- c("Blue", "Red") 
#color <- c("black","Skyblue","Blue", "Red") 
plot(out, col=color , main="Association of BRCA1 Mutations with Survival", xlab="Time,days", ylab="Proportion", lty=1:2, lwd=2)
legend("topright", levels(type), col=color, lty=1:2, lwd=3) 
