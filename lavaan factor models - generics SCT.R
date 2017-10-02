 ##############################################################################################
#
# Phenotype generator - takes input from REDcap and then generates factor scores based on data.
#
###############################################################################################
install.packages(c("data.table","scales","lavaan"))
library(devtools) # to run script from gist and bring in data
library(lavaan)
library(semPlot)
library(knitr)
library(ggplot2)
library(reshape2)
library(psych)


# 18-09-2017

require(tidyverse)

#P.thompson - adjusted for our genetic data.
#P.thompson - updated for data fixes in raw data - 02-10-2017 

##Reference: https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/ # code originally from this site.
mydir<-"C:/Users/pthompson/Dropbox/project SCT analysis/Extrapolated_norms/"
sct.data<-read.csv(paste0(mydir,"SCTData_DATA_2017-10-02_1038.csv"))
twin.data<-read.csv(paste0(mydir,"TwinsData_DATA_2017-10-02_1037.csv"))

short.sct<-select(sct.data,ï..record_id,age_at_test,wasi_matrices_ss,wasi_block_design_ss,wasi_vocab_ss,
                  wdck_jhsn_ss,sent_rep_ss,nonword_rep_ss,oromotor_ss,nara_acc_ss,nara_comp_ss,nara_rate_ss,
                  towre_words_ss,towre_nonwords_ss,phab_pic_ss,phab_digit_ss)
short.twin<-select(twin.data,record_id,age_at_test,wasi_matrices_ss,wasi_block_design_ss,wasi_vocab_ss,
                  wdck_jhsn_ss,sent_rep_ss,nonword_rep_ss,oromotor_ss,nara_acc_ss,nara_comp_ss,nara_rate_ss,
                  towre_words_ss,towre_nonwords_ss,phab_pic_ss,phab_digit_ss)
short.sct$source<-'sct'
short.twin$source<-'twin'

names(short.sct)[1]<-"record_id"

all.data<-rbind(short.sct,short.twin)
all.data$source<-as.factor(all.data$source)


all.data$record_id<-as.factor(all.data$record_id)

for(i in 1:length(all.data))
{
  all.data[,i]<-car::recode(all.data[,i],"996=NA;997=NA;998=NA;999=NA")
}

###############################################################################

all.data$missing<-vector(mode="numeric",length=531)

for(i in 1:length(all.data[,1]))
{
  all.data$missing[i]<-ifelse(any(is.na(all.data[i,4:15]))==TRUE,1,0)
}


#amount of missing data in all cells, not just rows

#sct now 4.8% missing
sum(is.na(all.data[389:531,5:16]))/(12*143)

#twin now 5.3% missing
sum(is.na(all.data[1:388,5:16]))/(12*388)

###############################################################################

# log transform 


log.gen <- all.data[, 3:16]
gen.set <- all.data[, 17]

################################################################################
################################################################################
# SEM program. 
################################################################################


library(reshape2)
#Which data for model

data.f1<-all.data[,3:16]

#plot variable distributions
mydata.melt <- melt(data.f1)
ggplot(mydata.melt, aes(x = value)) + 
  facet_wrap(~variable,scales = "free_x") +
  geom_density() 

#Single factor model

names(data.f1)<-c("matrices_ss", "blockD_ss","Vocab_ss", "WJ_ss", "Sent_Rep_ss",
 "Nword_Rep_ss", "Oro_ss", "Nara_acc", "Nara_comp", "Nara_rate", "Words_ss", 
 "Nwords_ss", "Pics_ss", "Digit_ss") 

data.f1$Oro_ss_f<-as.factor(data.f1$Oro_ss)


data.f2<-data.f1[-c(24,67),]

##########
model.f1 <- ' f1 =~ Vocab_ss + WJ_ss + Sent_Rep_ss + Oro_ss_f + Nword_Rep_ss
              f2 =~ Nara_acc + Nara_comp + Nara_rate + Words_ss + Nwords_ss + Pics_ss + Digit_ss
              f3 =~ matrices_ss + blockD_ss'
fit.mod.A <- sem(model.f1, data = data.f2,ordered=c("Oro_ss_f"),missing="pairwise")
summary(fit.mod.A, fit.measures = TRUE, standardized=TRUE, rsq=TRUE)

semPaths(fit.mod.A, "std", title = FALSE, 
         nCharNodes=0, edge.label.cex=0.6,esize=0.5)

fs1<-data.frame(record_id=all.data$record_id[-c(24,67)],lavPredict(fit.mod.A))

write.xlsx(fs1,file="factor_scores_SCT.xlsx",sheetName="Single_factor",row.name=F)

##########


model.f2 <- ' f1 =~ Vocab_ss + WJ_ss + Sent_Rep_ss + Oro_ss_f + Nword_Rep_ss
              f2 =~ Nara_acc + Nara_comp + Nara_rate + Words_ss + Nwords_ss + Pics_ss + Digit_ss
              f3 =~ matrices_ss + blockD_ss
              m =~ f1 + f2 + f3'

fit.mod.B <- sem(model.f2, data = data.f2,ordered=c("Oro_ss_f"),missing="pairwise")

#,ordered=c("Oro_ss_f"))
summary(fit.mod.B, fit.measures = TRUE, standardized=TRUE, rsq=TRUE)

semPaths(fit.mod.B, "std", title = FALSE, 
         nCharNodes=0, edge.label.cex=0.6,esize=0.5)

fs2<-data.frame(record_id=all.data$record_id[-c(24,67)],lavPredict(fit.mod.B))

write.xlsx(fs2,file="factor_scores_SCT.xlsx",sheetName="Second_order CFA",row.names = F,append = TRUE)
#######

model.f3 <- ' f1 =~ matrices_ss + blockD_ss + Vocab_ss + WJ_ss + Sent_Rep_ss + Oro_ss_f + Nword_Rep_ss
              f2 =~ Nara_acc + Nara_comp + Nara_rate + Words_ss + Nwords_ss + Pics_ss + Digit_ss'
fit.mod.C <- sem(model.f3, data = data.f2,ordered=c("Oro_ss_f"),missing="pairwise")
summary(fit.mod.C, fit.measures = TRUE, standardized=TRUE, rsq=TRUE)

semPaths(fit.mod.C, "std", title = FALSE, 
         nCharNodes=0, edge.label.cex=0.6,esize=0.5)

fs3<-data.frame(record_id=all.data$record_id[-c(24,67)],lavPredict(fit.mod.C))

write.xlsx(fs3,file="factor_scores_SCT.xlsx",sheetName="Two_factor corr CFA",row.names = F,append = TRUE)

##########

#######

library(nFactors)

ev <- eigen(cor(data.f1[,1:14],use= "na.or.complete")) # get eigenvalues
ap <- parallel(subject=nrow(data.f1[,1:14]),var=ncol(data.f1[,1:14]),
  rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
