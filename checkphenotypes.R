#Input: output from phenoscanner (www.phenoscanner.medschl.cam.ac.uk/)
#Output: three files with the number of times the direction of effect of all phenotypes available and a subset of the phenotypes (those of higher interest), per SNP

#GOAL
#Obtain a landscape of phenotypes associated with a list of predetermined SNPs, by counting the number of concordant and discordant effects from published GWAS available in Phenoscanner

#USAGE:
#R CMD BATCH checkphenotypes.R

#necessary package: plyr - to install in R: install.packages("plyr")
library(plyr)

#input file obtained from Phenoscanner
database<-read.table("phenoscanner_example.txt", sep="\t", header=T)

#get the list of SNPs from the Phenoscanner input
esnipes<-unique(database$SNP)
#get the list of phenotypes that show some association according to Phenoscanner
phenotypes<-unique(database$Trait)
#create and initialize the output matrix
plusminus<-matrix(0,nrow=((length(phenotypes)*2)+1),ncol=length(esnipes)+1)
princpheno<-0
for (i in 1:length(phenotypes)){
plusminus[i*2,1] = toString(phenotypes[i]);
}

#loop to select the lowest p-value of the phenotype in the LD range predetermined at Phenoscanner
#only considers phenotypes with beta available, in order to consider direction of effect
for (tala in 1:length(esnipes)){
	plusminus[1,tala+1]=toString(esnipes[tala])
	subesnipes<-subset(database,database$SNP==esnipes[tala])
	subphenotypes<-unique(subesnipes$Trait)
	for (doenca in 1:length(subphenotypes)){
		subna<-subset(subesnipes, !is.na(subesnipes$Beta))
		subsubesnipes<-subset(subna,subna$Trait==subphenotypes[doenca])
		l<-match(subphenotypes[doenca],plusminus[,1],incomparables = NULL)
		minim<-which.min(subsubesnipes$P)
		if (length(minim)>0){
				if(subsubesnipes$Beta[minim]<0){
					plusminus[l+1,tala+1]=as.numeric(plusminus[l+1,tala+1])+1
				}else if(subsubesnipes$Beta[minim]>0){
					plusminus[l,tala+1]=as.numeric(plusminus[l,tala+1])+1
			}
		}
	}
}

#list of phenotypes of particular interest - can be altered to fir other purposes
tier1 <-c("Coronary artery disease", "LDL cholesterol", "LDL", "Total cholesterol", "Triglycerides", "Systolic blood pressure SBP", "Diastolic blood pressure DBP", "Hypertension", "Type II diabetes", "Obesity class 2", "Obesity class 1", "Obesity class 3", "BMI", "Body mass index BMI", "Obesity body mass index BMI")
important<-matrix("",nrow=((length(tier1)*2)+1),ncol=length(esnipes)+1)
important[1,]<-plusminus[1,]
#loop to count directions of effect only of the phenotypes of interest
for (i in 1:length(tier1)){
			l<-match(tier1[i],plusminus[,1],incomparables = NULL)
			important[2*i,]<-plusminus[l,]
			important[2*i+1,]<-plusminus[l+1,]
			if (i==1){
				princpheno<-l
			}
}

#count directions of effect for all subset of important phenotype available in the input
vitorias<-matrix(0,nrow=(length(esnipes)+1),ncol=4)

for (l in 1:(nrow(vitorias)-1)){
	v<-match(esnipes[l],important[1,],incomparables = NULL)
	vitorias[l+1,1]<-as.character(esnipes[l])
	g<-match(esnipes[l],database$SNP,incomparables = NULL)
	vitorias[l+1,4]<-as.character(database$gene[g])
	for (i in seq(4,(2*length(tier1)-1),2)){
			if (important[2,v]>important[3,v]){
				if (important[i,v]>important[i+1,v]){
					vitorias[l+1,2]=as.integer(vitorias[l+1,2])+1
				} else if (important[i,v]<important[i+1,v]){
					vitorias[l+1,3]=as.integer(vitorias[l+1,3])+1
				}
			} else if (important[2,v]<important[3,v]){
				if (important[i,v]>important[i+1,v]){
					vitorias[l+1,3]=as.integer(vitorias[l+1,3])+1
				} else if (important[i,v]<important[i+1,v]){
					vitorias[l+1,2]=as.integer(vitorias[l+1,2])+1
				}
			}
	}
}

#count directions of effect for all phenotypes available in the input
totvitorias<-matrix(0,nrow=(length(esnipes)+1),ncol=4)

for (l in 1:(nrow(totvitorias)-1)){
	v<-match(esnipes[l],plusminus[1,],incomparables = NULL)
	totvitorias[l+1,1]<-as.character(esnipes[l])
	g<-match(esnipes[l],database$SNP,incomparables = NULL)
	#print (g)
	totvitorias[l+1,4]<-as.character(database$gene[g])
	for (i in seq(4,(nrow(plusminus)-1),2)){
			if (plusminus[princpheno,v]>plusminus[princpheno+1,v]){
				if (plusminus[i,v]>plusminus[i+1,v]){
					totvitorias[l+1,2]=as.integer(totvitorias[l+1,2])+1
				} else if (plusminus[i,v]<plusminus[i+1,v]){
					totvitorias[l+1,3]=as.integer(totvitorias[l+1,3])+1
				}
			} else if (plusminus[princpheno,v]<plusminus[princpheno+1,v]){
				if (plusminus[i,v]>plusminus[i+1,v]){
					totvitorias[l+1,3]=as.integer(totvitorias[l+1,3])+1
				} else if (plusminus[i,v]<plusminus[i+1,v]){
					totvitorias[l+1,2]=as.integer(totvitorias[l+1,2])+1
				}
			}
	}
}


#writeoutputs - one with the directions of effect per phenotype, one with the number of selected phenotypes concordant/discordant, and one with all phenotypes concordant/discordant
write.table(important,sep=",","gwastableplusminus.txt",col.names = F, row.names = F, append=F)
write.table(vitorias,sep=",","gwastablevitorias.txt",col.names = F, row.names = F, append=F)
write.table(totvitorias,sep=",","gwastabletotvitorias.txt",col.names = F, row.names = F, append=F)