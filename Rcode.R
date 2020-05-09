
#Descriptive Epi Diversity infant-mother
#DIFFERENCES IN SUM OF RELATIVE ABUNDANCE BETWEEN MOTHER AND INFANTS
library("ggplot2")
library("ggpubr")
library("vegan")
library("dplyr")
library(ggsci)
library("VennDiagram")
library(gridExtra)
library(grid)

######to set the path of the new working directory
setwd("/Users/andreasosamoreno/Desktop/LixinData/PlosManuscript/Files")

###### Load abundance data
abundance <- read.csv ("Abundance.csv", header=TRUE)
abundance <-abundance [,c("ARCHgutm","class","pair","X26","X201","X202","X206","X336","X338","X340","X342","X359","X366","X370","X372","X376","X1522","X1528","X1546","X1547","X1548","X1551","X1552","X1553","X1556","X1557","X1558","X1559","X1565","X1567","X6","X14","X49","X104","X174","X402","X404","X406","X410","X412","X417","X425","X428","X429","X431","X432","X438","X1503","X1540","X1541","X1545","X9","X11","X64","X81","X89","X234","X245","X246","X298","X331","X355","X1300","X1302","X1303","X1305","X1504","X1509","X1536","X1549","X1572","X1573","X42","X46","X106","X107","X108","X121","X153","X162","X236","X362","X1108","X1118","X1123","X1505","X1512","X1544","X54","X180","X181","X185","X191","X196","X200","X294","X506","X507","X1507","X1513","X1539","X91","X137","X138","X209","X227","X229","X283","X285","X801","X804","X806","X809","X812","X815","X817","X819","X1511","X1519","X133","X177","X208","X280","X363","X211","X213","X214","X215","X306","X309","X316","X318","X328","X1200","X1201","X1577")]

abundance$class<-as.factor(ifelse(abundance$class=="Child","Infancy","Pregnancy"))
abundance [is.na(abundance)]<- 0
abundance_mom<-abundance[which(abundance$class=="Pregnancy"),]
abundance_child<-abundance[which(abundance$class=="Infancy"),]

###### Load presence/absence data
presence<-read.csv("Binary.csv", header=T)
presence<-presence[,c("ARCHgutm","class","pair","X26","X201","X202","X206","X336","X338","X340","X342","X359","X366","X370","X372","X376","X1522","X1528","X1546","X1547","X1548","X1551","X1552","X1553","X1556","X1557","X1558","X1559","X1565","X1567","X6","X14","X49","X104","X174","X402","X404","X406","X410","X412","X417","X425","X428","X429","X431","X432","X438","X1503","X1540","X1541","X1545","X9","X11","X64","X81","X89","X234","X245","X246","X298","X331","X355","X1300","X1302","X1303","X1305","X1504","X1509","X1536","X1549","X1572","X1573","X42","X46","X106","X107","X108","X121","X153","X162","X236","X362","X1108","X1118","X1123","X1505","X1512","X1544","X54","X180","X181","X185","X191","X196","X200","X294","X506","X507","X1507","X1513","X1539","X91","X137","X138","X209","X227","X229","X283","X285","X801","X804","X806","X809","X812","X815","X817","X819","X1511","X1519","X133","X177","X208","X280","X363","X211","X213","X214","X215","X306","X309","X316","X318","X328","X1200","X1201","X1577")]

presence$class<-as.factor(ifelse(presence$class=="Child","Infancy","Pregnancy"))
presence_mom<-presence[which(presence$class=="Pregnancy"),]
presence_child<-presence[which(presence$class=="Infancy"),]

###### Demographics
###Load maternal demographics
demomom<-read.csv("demographics_pregnancy.csv")

###Load infancy demographics
demochild<-read.csv("demographics_infancy.csv")

###### Richness and abundance of ARG and MGE genes among fecal samples

###Number of genes OVERALL
overallsum<- colSums(presence[,31:136])
momsum<- colSums(presence_mom[,31:136])
childsum<- colSums(presence_child[,31:136])
which(childsum==0)


###number of genes MOTHERS same code for infants
#Number of ARGs
momsum<- rowSums(presence_mom[,31:136])
max(momsum)
min(momsum)

#Number of MGEs
momsum<- rowSums(presence_mom[,4:30])
which(overallsum==0)
max(momsum)
min(momsum)


###median relative abundance
#### To sum all RA of al ARG and MGE
abundance$totalRA_MGE <- rowSums(abundance[,4:30],na.rm=T)
abundance$totalRA_ARG <- rowSums(abundance[,31:136],na.rm=T)
medianARG<- median(abundance$totalRA_ARG)
medianMGE<- median(abundance$totalRA_MGE)
max(abundance$totalRA_ARG)
min(abundance$totalRA_ARG)
max(abundance$totalRA_MGE)
min(abundance$totalRA_MGE)

###Overall per class
indices <- abundance[,c("class","pair")]
rownames(indices)<-rownames(abundance)
mydata_RA<-abundance[4:136]

indices$SumMGE<-rowSums(mydata_RA[,1:27])
indices$Sumamino <- rowSums(mydata_RA[,28:48])
indices$SumMDR<-rowSums(mydata_RA[,49:69])
indices$Sumbeta<-rowSums(mydata_RA[,70:85])
indices$Sumtetra<-rowSums(mydata_RA[,86:98])
indices$SumMLSB<-rowSums(mydata_RA[,99:116])
indices$Sumsulfo<-rowSums(mydata_RA[,117:121])
indices$Sumvanco<-rowSums(mydata_RA[,122:129])
indices$Sumfluoro<-rowSums(mydata_RA[,130:133])
indices$SumARG<-rowSums(mydata_RA[,28:133])


#Test for differences between mother and child
wilcox.test(SumARG~class, data=indices) #0.1003
wilcox.test(SumMGE~class, data=indices) #0.0006

#######Shared ARG and MGE patterns between pregnant women and infants
multiplepie<-read.csv("SharedGenesPie.csv")
#### Median of shared genes between mothers and infants, 
multiplepie2<-multiplepie[which(multiplepie$Excel=="Common"),]
median(multiplepie2$Percentage)
max(multiplepie2$Percentage)
min(multiplepie2$Percentage)

###Figure 2
###Venn Diagram
#1st row
bg12<-draw.pairwise.venn(100, 92, 34, fill=c("red","aquamarine"), cex=1)
bg9<-draw.pairwise.venn(120, 122, 58, fill=c("red","aquamarine"), cex=1, rotation.degree=180)
b2016<-draw.pairwise.venn(103, 106, 41, fill=c("red","aquamarine"), cex=1, rotation.degree=180)
b3032<-draw.pairwise.venn(87, 100, 30, fill=c("red","aquamarine"), cex=1, rotation.degree=180)
bg1<-draw.pairwise.venn(97, 42, 18, fill=c("red","aquamarine"), cex=1)
bg18<-draw.pairwise.venn(153, 151, 73, fill=c("red","aquamarine"), cex=1)


venn<-grid.arrange(grobTree(bg12),grobTree(bg9),grobTree(b2016), grobTree(b3032),grobTree(bg1),grobTree(bg18),ncol=6, heights = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), nrow = 6)
ggsave("multivenn.pdf", venn)

#######Characteristics of the pregnancy and the infancy resistome
#Sum of all antibiotic resistance genes found per class
mydata<-presence[,4:136]
indices$RichMGE<-rowSums(mydata[,1:27])
indices$RichARG<-rowSums(mydata[,28:133])

#Test for differences between mother and child
wilcox.test(RichARG~class, data=indices) #0.78
wilcox.test(RichMGE~class, data=indices) #5.88x10-6

###Fig 3
barplot<-c(indices$Sumamino,indices$SumMDR, indices$Sumbeta, indices$Sumtetra, indices$SumMLSB, indices$Sumsulfo, indices$Sumvanco, indices$Sumfluoro,indices$SumMGE)
type<-rep(c("aminoglycoside*","MDR*","betalactamase","tetracycline","MLSB*","sulfonamide*","vancomycin*","fluoroquinolone*","MGE*"),each=91)
class<-rep(abundance$class,9)
dataframe<-data.frame(barplot, type, class)
dataframe2<-aggregate(barplot~type+class, dataframe,median)
mother<-ggplot(data = dataframe2, aes(x = class, y=barplot)) +theme_grey(base_size = 10)+ geom_bar(stat="identity", aes(fill=type))+theme_classic() + scale_fill_brewer(palette="Paired")+ labs(title = "", y = "Median (Relative abundance)", x = "") + labs(fill = "Antibiotic class")
ggsave(mother, filename = "Descriptive/Fig3.tiff", units="in", height=4, width=4,dpi=300, compression="lzw")
dev.off()

######Infancy resistome is more diverse than pregnancy resistome: Alpha diversity
indices$ShannonMGERA<-diversity(mydata_RA[,1:27])
indices$ShannonARGRA<-diversity(mydata_RA[,28:133])
a<- ggboxplot(data=indices, x="class", y="RichMGE", color="black", fill="class", palette="d3", ylab="Richness: Number of genes",xlab="", bxp.errorbar=T, outlier.colour = "gray", outlier.shape = 1, main="MGE")+ theme_classic(base_size=9) + guides(fill=FALSE) + coord_cartesian(ylim = c(0, 22))+ stat_compare_means(label.x=0.8, label.y=22, size=3)
b<- ggboxplot(data=indices, x="class", y="RichARG", color="black", fill="class", palette="d3", ylab="Richness: Number of genes",xlab="", bxp.errorbar=T, outlier.colour = "gray", outlier.shape = 1, main="ARG")+ theme_classic(base_size=9) + guides(fill=FALSE) + coord_cartesian(ylim = c(0, 80))+ stat_compare_means(label.x=0.8, label.y=75, size=3)
c<- ggboxplot(data=indices, x="class", y="ShannonMGERA", color="black", fill="class", palette="d3", ylab="Shannon diversity index",xlab="", bxp.errorbar=T, outlier.colour = "gray", outlier.shape = 1, main="MGE")+ theme_classic(base_size=9) + guides(fill=FALSE) + coord_cartesian(ylim = c(0, 2.2))+ stat_compare_means(label.x=0.8, label.y=2.2, size=3)
d<- ggboxplot(data=indices, x="class", y="ShannonARGRA", color="black", fill="class", palette="d3", ylab="Shannon diversity index",xlab="", bxp.errorbar=T, outlier.colour = "gray", outlier.shape = 1, main="ARG")+ theme_classic(base_size=9) + guides(fill=FALSE) + coord_cartesian(ylim = c(0, 3))+ stat_compare_means(label.x=0.8, label.y=3, size=3)
dev.off()
figure <- ggarrange(a, b, c, d, labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))
ggsave(figure, filename = "Diversity/Fig4.tiff", units="in", height=5, width=5,dpi=300, compression="lzw")

######Pregnancy resistome differs from infancy resistome: Beta diversity
###PCOA beta diversity and maternal covarites EXAMPLE
###Vegan vegdist command
spe.bray <- vegdist(mydata_RA[,28:133], method="bray")
spe.brayMGE <- vegdist(mydata_RA[,1:27], method="bray")
spe.b.pcoa <-cmdscale(spe.bray,k=(nrow(mydata[,28:133])-1),eig=TRUE)
spe.b.pcoaMGE <-cmdscale(spe.brayMGE,k=(nrow(mydata[,1:27])-1),eig=TRUE)

###PCOA plot for ARG Bray Curtis Index Figure 5
#For MGE plots use spe.b.pcoaMGE instead of spe.b.pcoa
#For other Sorensen plots use the dataset presence (presence/absence data)
eig.Data.Sorensen.PCoA<-spe.b.pcoa$eig
  eig.Data.Sorensen.PCoA.sum<-sum(eig.Data.Sorensen.PCoA)
  a<-(eig.Data.Sorensen.PCoA/eig.Data.Sorensen.PCoA.sum)*100
plot(scores(spe.b.pcoa, choices=c(1,2)), type="n", cex=0.7, ann=FALSE, cex.axis=0.6, )
title(xlab=paste("PC1","(",round(a[1],1),"%",")",sep=""), ylab=paste("PC2","(",round(a[2],1),"%",")",sep=""),line=2, cex.lab=0.6)
ordiellipse(spe.b.pcoa, abundance$class, draw="lines", col=c("red","blue"))
points(spe.b.pcoa$points[which(abundance$class=="Infancy"),], pch=21, col="red", bg="red", cex=0.6)
points(spe.b.pcoa$points[which(abundance$class=="Pregnancy"),], pch=21, col="blue", bg="blue", cex=0.6)
legend(x="topleft", legend=c("Infancy","Pregnancy"), col=c("red","blue"), pt.bg=c("red","blue"), pch=21, cex=0.6) 
text(x = 0.15, y = .41, labels = "p value < 0.001", xpd = NA, cex=0.6)
text(x = -0.245, y = .55, labels = "A     ARG Bray Curtis Index", xpd = NA, cex=0.8, font=2)

### To get p values from the permanova test
#PERMANOVA ANALYSIS p=0.0001
adonis<-adonis(spe.bray~class, abundance, permutations=9999, method="bray", p.adjust="fdr" ) 

######Demographics and perinatal influencing resistome patterns Fig 6 and 7
###Maternal covariates
###Merge demographics and resistome data
dataset_ab<-merge(demomom, abundance_mom, by.y="ARCHgutm",by.x="Study.ID")
dataset_pa<-merge(demomom, presence_mom, by.y="ARCHgutm",by.x="Study.ID")

###Creating new variables
dataset_pa$race<- factor(ifelse(dataset_pa$BC_RACEMOM=="White","White","Non-White"))
dataset_ab$race<- factor(ifelse(dataset_ab$BC_RACEMOM=="White","White","Non-White"))
nrow(dataset_pa)
dataset_ab$BMI <- factor(ifelse(dataset_ab$mom_Prepreg_BMI<25,"Normal or underweight",ifelse(dataset_ab$mom_Prepreg_BMI>=25& dataset_ab$mom_Prepreg_BMI<30,"Overweight","Obese")))

dataset_ab$parity <- factor(ifelse(dataset_ab$E_children_given_birth_to_inc_current<3,"1-2 children",">3 children"))

###Prepare data for PCOA
y<-dataset_ab$race
dataset3<-dataset_ab[which(y!="NA"),] #Final dataset for abundance
dataset4<-dataset_pa[which(y!="NA"),] #Final dataset for Binary data
x<-dataset3$race

#Relative
amrdata_RAARG<- dataset3[,42:147]
amrdata_RAMGE<- dataset3[,15:41]
spe.bray <- vegdist(amrdata_RAARG, method="bray")
spe.brayMGE <- vegdist(amrdata_RAMGE, method="bray")
spe.b.pcoa <-cmdscale(spe.bray,k=(nrow(amrdata_RAARG)-1),eig=TRUE)
spe.b.pcoaMGE <-cmdscale(spe.brayMGE,k=(nrow(amrdata_RAMGE)-1),eig=TRUE)
eig.Data.Sorensen.PCoA<-spe.b.pcoa$eig
  eig.Data.Sorensen.PCoA.sum<-sum(eig.Data.Sorensen.PCoA)
  a<-(eig.Data.Sorensen.PCoA/eig.Data.Sorensen.PCoA.sum)*100

### PCOA ARG maternal race
#For other ARG maternal variables use x=dataset_ab$BMI (pre-pregancy weight), dataset_ab$E_ever_smoked (smoking status) and dataset_ab$parity (parity)
ordiplot(scores(spe.b.pcoa, choices=c(1,2)), type="n", xlab=paste("PC1","(",round(a[1],1),"%",")",sep=""), ylab=paste("PC2","(",round(a[2],1),"%",")",sep="") )
ordiellipse(spe.b.pcoa, x, draw="lines", col=c("red","blue"))
points(spe.b.pcoa$points[which(x=="Non-White"),], pch=21, col="red", bg="red")
points(spe.b.pcoa$points[which(x=="White"),], pch=21, col="blue", bg="blue")
legend(x="topleft", legend=levels(x), col=c("red","blue"), pt.bg=c("red","blue"), pch=21, cex=0.8, box.lty=0) 
text(x = -0.45, y = -0.27, labels = "p value = 0.011", xpd = NA, cex=0.8)
text(x = -0.67, y = 0.56, labels = "A     Maternal Race", xpd = NA,cex=1.3, font=2)

adonis(spe.bray~x, dataset3, permutations=9999, method="bray", p.adjust="frd")

###Infant covariates
###Merge demographics and resistome data
dataset_ab<-merge(demochild, abundance_child, by.y="ARCHgutm",by.x="Study.ID")
dataset_pa<-merge(demochild, presence_child, by.y="ARCHgutm",by.x="Study.ID")

###Creating new variables
dataset_ab$breastfeeding <- factor(ifelse(dataset_ab$I6mosSIF_breastmilk_in_diet=="0%"|dataset_ab$I6mosSIF_breastmilk_in_diet=="20%", "<50%",">=50%"))
dataset_pa$breastfeeding <- factor(ifelse(dataset_pa$I6mosSIF_breastmilk_in_diet=="0%"|dataset_pa$I6mosSIF_breastmilk_in_diet=="20%", "<50%",">=50%"))

dataset_ab$solid <- factor(ifelse(dataset_ab$I6M_Infant_eating_solids=="0","non-solids","solids"))
dataset_pa$solid <- factor(ifelse(dataset_pa$I6M_Infant_eating_solids=="0","non-solids","solids"))

###Prepare data for PCOA
y<-dataset_ab$BC_sex
dataset3<-dataset_ab[which(y!="NA"),] #Final dataset for abundance
dataset4<-dataset_pa[which(y!="NA"),] #Final dataset for Binary data
x<-dataset3$BC_sex
str(amrdata_RAARG)
#Relative
amrdata_RAARG<- dataset4[,44:149]
amrdata_RAMGE<- dataset4[,17:43]
spe.bray <- vegdist(amrdata_RAARG, method="bray")
spe.brayMGE <- vegdist(amrdata_RAMGE, method="bray")
spe.b.pcoa <-cmdscale(spe.bray,k=(nrow(amrdata_RAARG)-1),eig=TRUE)
spe.b.pcoaMGE <-cmdscale(spe.brayMGE,k=(nrow(amrdata_RAMGE)-1),eig=TRUE)
eig.Data.Sorensen.PCoA<-spe.b.pcoa$eig
  eig.Data.Sorensen.PCoA.sum<-sum(eig.Data.Sorensen.PCoA)
  a<-(eig.Data.Sorensen.PCoA/eig.Data.Sorensen.PCoA.sum)*100

### PCOA ARG infant sex
#For other ARG maternal variables use x=dataset_ab$BMI (pre-pregancy weight), dataset_ab$E_ever_smoked (smoking status) and dataset_ab$parity (parity)
ordiplot(scores(spe.b.pcoa, choices=c(1,2)), type="n", xlab=paste("PC1","(",round(a[1],1),"%",")",sep=""), ylab=paste("PC2","(",round(a[2],1),"%",")",sep="") )
ordiellipse(spe.b.pcoa, x, draw="lines", col=c("red","blue"))
points(spe.b.pcoa$points[which(x=="Female"),], pch=21, col="red", bg="red")
points(spe.b.pcoa$points[which(x=="Male"),], pch=21, col="blue", bg="blue")
legend(x="topleft", legend=levels(x), col=c("red","blue"), pt.bg=c("red","blue"), pch=21, cex=0.8, box.lty=0) 
text(x = 0.15, y = 0.25, labels = "p value = 0.183", xpd = NA, cex=0.8)
text(x = -0.69, y = 0.56, labels = "A     Sex", xpd = NA,cex=1.3, font=2)

adonis(spe.bray~x, dataset4, permutations=9999, method="bray", p.adjust="frd")

###Multivariable analysis for Relative abudance data
#adjust for Maternal Age, Cohort, Shipping Time, BMI Category
##Change to MGE in spe.brayMGe
##Change dataset_ab to dataset_pa to get data on composition
drops<-c("E_ever_smoked","E_times_taken_abx_past_year","E_abx_names","ASIF_taken_Abx_past_mos")
multivariable<-dataset_ab%>%
  select(-one_of(drops))

multivariable<-na.omit(multivariable)
Cohort<-ifelse(grepl("B_", multivariable$Study.ID),1,2)
nrow(multivariable)
amrdata_RAARG<- multivariable[,38:143]
amrdata_RAMGE<- multivariable[,11:37]

#Bray Curtis Matrix for maternal samples
spe.bray <- vegdist(amrdata_RAARG, method="bray")
spe.brayMGE <- vegdist(amrdata_RAMGE, method="bray")
spe.b.pcoa <-cmdscale(spe.bray,k=(nrow(amrdata_RAARG)-1),eig=TRUE)
spe.b.pcoaMGE <-cmdscale(spe.brayMGE,k=(nrow(amrdata_RAMGE)-1),eig=TRUE)
str(multivariable$Sample_Transit_Time_days)
str(multivariable$mom_Prepreg_BMI)
str(Cohort)
str(multivariable$Mom_age_Years)
str(multivariable$race)


adonis2(spe.bray~as.numeric(multivariable$Sample_Transit_Time_days)+as.numeric(multivariable$mom_Prepreg_BMI)+as.numeric(multivariable$Mom_age_Years)+as.factor(Cohort)+as.factor(multivariable$race),by=NULL,permutations=9999)

adonis2(spe.brayMGE~as.numeric(multivariable$Sample_Transit_Time_days)+as.numeric(multivariable$mom_Prepreg_BMI)+as.numeric(multivariable$Mom_age_Years)+as.factor(Cohort)+as.factor(multivariable$race),by="margin",permutations=9999)

#adjust for breastfeeding, sex, Shipping Time, cohort, mode of delivery
##Change to MGE in spe.brayMGe
##Change dataset_ab to dataset_pa to get data on composition
drops<-c("infant_birth_weight_grams_BC","I6mosSIF_past_month_Medicine","I1WSIF_since_birth_Medicine","I6M_Infant_eating_solids","DueDate_MMDDYEAR")
multivariable<-dataset_pa%>%
  select(-one_of(drops))

multivariable<-na.omit(multivariable)
Cohort<-ifelse(grepl("B_", multivariable$Study.ID),1,2)
str(multivariable)
amrdata_RAARG<- multivariable[,39:144]
amrdata_RAMGE<- multivariable[,12:38]

#Bray Curtis Matrix for maternal samples
spe.bray <- vegdist(amrdata_RAARG, method="bray")
spe.brayMGE <- vegdist(amrdata_RAMGE, method="bray")
spe.b.pcoa <-cmdscale(spe.bray,k=(nrow(amrdata_RAARG)-1),eig=TRUE)
spe.b.pcoaMGE <-cmdscale(spe.brayMGE,k=(nrow(amrdata_RAMGE)-1),eig=TRUE)
str(multivariable$Sample_Transit_Time_days)
str(multivariable$breastfeeding)
str(Cohort)
str(multivariable$BC_sex)
str(multivariable$ROUT)


adonis2(spe.bray~as.numeric(multivariable$Sample_Transit_Time_days)+as.factor(multivariable$breastfeeding)+as.factor(multivariable$BC_sex)+as.factor(Cohort)+as.factor(multivariable$ROUT),by=NULL,permutations=9999)

adonis2(spe.brayMGE~as.numeric(multivariable$Sample_Transit_Time_days)+as.factor(multivariable$solid)+as.factor(multivariable$BC_sex)+as.factor(Cohort)+as.factor(multivariable$ROUT),by="margin",permutations=9999)








