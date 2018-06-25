# Copyright (C) 2015  Marco Manfrini, PhD

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------------------------------------------------------
# LIBRARY
------------------------------------------------------------------------------------------------------------------------

library(glmnet)
library(corrplot)
library(survival)

source("/home/marco/Dati/codice/script-R/Functions.R")

------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
------------------------------------------------------------------------------------------------------------------------

plot_desc<-function(dset){
  
  # PLOT SET AFTER PREPROCESSING
  
  jpeg("Plot1.jpeg", width = 210, height = 297, units = 'mm', res = 300)
  pdf("Plot1.pdf", width = 8.3, height = 11.7)
  par(mfrow=c(6,5))
  
  hist(dset$ETA , ylim = c(0, 200), main = "ETA' distribution", xlab = "ETA'", ylab = "COUNTS", border = "white", col = "orange")
  barplot(table(dset$FIGO.CODE.2), ylim=c(0, 500), main="FIGO code 2", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$FIGO.CODE), ylim=c(0, 700), main="FIGO code", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$MIOMETRIO), ylim=c(0, 500), main="MIOMETRIO", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$CANALE), ylim=c(0, 700), main="CANALE", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$WASHING), ylim=c(0, 800), main="WASHING", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$VAGINA), ylim=c(0, 900), main="VAGINA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$ANNESSI), ylim=c(0, 800), main="ANNESSI", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$SIEROSA), ylim=c(0, 800), main="SIEROSA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$PARAMETRI), ylim=c(0, 800), main="PARAMETRI", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$CHIR.SIGLA), ylim=c(0, 800), main="CHIR SIGLA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$LINFADENECTOMIA), ylim=c(0, 800), main="LINFADENECTOMIA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$LINFONODI), ylim=c(0, 800), main="LINFONODI", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$CARCINOSI), ylim=c(0, 800), main="CARCINOSI", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  hist(dset$BMI, ylim=c(0, 300), main="BMI", xlab="BMI", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$MENOPAUSA), ylim=c(0, 800), main="MENOPAUSA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  hist(dset$ETAMENOPA, ylim=c(0, 400), main="ETAMENOPA", xlab="ETAMENOPA", ylab="COUNTS", border="white", col="orange")
  hist(dset$ANNIMENO, ylim=c(0, 200), main="ANNIMENO", xlab="ANNIMENO", ylab="COUNTS", border="white", col="orange")
  hist(dset$MENARCA, ylim=c(0, 300), main="MENARCA", xlab="MENARCA", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$PARITA), ylim=c(0, 300), main="PARITA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$GRAVIDANZE), ylim=c(0, 400), main="GRAVIDANZE", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$NEOADIUV), ylim=c(0, 1000), main="NEOADIUV", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$ADIUV), ylim=c(0, 500), main="ADIUV", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$CHEMIO.CODE), ylim=c(0, 800), main="CHEMIO", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$ORMONO), ylim=c(0, 800), main="ORMONO", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$RXTERAP), ylim=c(0, 600), main="RXTERAP", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$ENDO), ylim=c(0, 700), main="ENDO", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$SINT), ylim=c(0, 500), main="SINTOMO", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  hist(dset$DURATA, ylim=c(0, 600), main="DURATA", xlab="DURATA", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$SEDE.SIGLA1), ylim=c(0, 800), main="SEDE.SIGLA1", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  
  dev.off()
  
  jpeg("Plot2.jpeg", width = 210, height = 297, units = 'mm', res = 300)
  pdf("Plot2.pdf", width = 8.3, height = 11.7)
  par(mfrow=c(6,5))
  
  barplot(table(dset$ISTOTIPO.CODE), ylim=c(0, 800), main="ISTOTIPO.CODE", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$GRADING), ylim=c(0, 400), main="GRADING", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  barplot(table(dset$TIPO.CA), ylim=c(0, 700), main="TIPO.CA", xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
  
  dev.off()
  return(TRUE)
}

plot_set<-function(dset){
  
  # PLOT SET AFTER PREPROCESSING

  nvars<-dim(dset)[2]
  
  pdf("dataset_pp_v2.pdf", width = 8.3, height = 11.7)
    
  par(mfrow=c(6,5))
  
  for(i in 1:nvars){
  
    barplot(table(dset[i]), main=names(dset[i]), xlab="LEVELS", ylab="COUNTS", border="white", col="orange")
    
  }

  dev.off()
  
}

plot_surv<-function(dset){
  
  nvars<-dim(dset)[2]
  c<-4
  r<-ceiling(nvars/c)
  inc<-2*c
  inr<-2*r
  #pdf("surv.pdf", width = 8.3, height = 11.7)
  pdf("surv.pdf", width = inc, height = inr)
  
  par(mfrow=c(r,c))
  
  for(i in 1:(nvars-2)){
    
    surv<-survfit(Surv(FUM, STATUS)~dset[,i], data=dset)
    plot(surv, mark.time=F, col=c("red", "blue", "green", "orange", "black"), main=names(dset[i]))  
    legend(x=60, y=1, legend=levels(as.factor(dset[,i])), lty=c(1,1), col=c("red", "blue", "green", "orange", "black"), bty="n", cex=0.3)
  }
  
  dev.off()
  
}

stats_categorical<-function(x){

  nr<-length(x)
  m<-t(as.matrix(round(table(x)/nr*100, digits=2)))
  row.names(m)<-"%"
  return(m)
    
}

calc_surv<-function(x, s){
  
  n<-levels(x$s)
  print(n)
  surv<-survfit(Surv(FUM, STATUS)~ strata(s), data=x)
  plot(surv, mark.time=F, col=c("red", "blue", "green", "orange"))
  return(surv)
  
}

#------------------------------------------------------------------------------------------------------------------------
# LOAD DATA
#------------------------------------------------------------------------------------------------------------------------

setwd("/home/marco/Dati/Borghi/patterndirecidivecarcinomaendometriale")

dset <- read.delim("ENDOMETRIO ANONIMO PER ELABORAZIONE.csv" , sep = "\t" , dec = "," , stringsAsFactors = F)

#------------------------------------------------------------------------------------------------------------------------
# DESCRIPTIVE
#------------------------------------------------------------------------------------------------------------------------

#plot_desc(dset)

# REMOVE VARAIBLES 95% UNIMODAL

outVector<-NULL
for(i in 1:dim(dset)[2]){
  
  m<-as.data.frame(stats_categorical(dset[,i]))
  sel<-which(m>95)
  if(length(sel>0)){outVector<-c(outVector, names(dset[i]))
#  write.table(m, file=paste0(getwd(), "/", names(dset[i]), ".95.tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)
    }
#  write.table(m, file=paste0(getwd(), "/", names(dset[i]), ".tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)

}

dset<-dset[, !(colnames(dset) %in% outVector)]

#write.table(outVector, file=paste0(getwd(), "/", "removed.tsv"), sep="\t", col.names=F, row.names=F, quote=F)

#------------------------------------------------------------------------------------------------------------------------
# RECODE VARIABLES 70% UNIMODAL or <10%(5%) IN AT LAST ONE MODALITY
#------------------------------------------------------------------------------------------------------------------------
  
outVector<-NULL
for(i in 1:dim(dset)[2]){
  
  m<-as.data.frame(stats_categorical(dset[,i]))
  sel<-which(m>70 | m<5)
  if(length(sel>0)){outVector<-c(outVector, names(dset[i]))
#    write.table(m, file=paste0(getwd(), "/", names(dset[i]), ".70.5.tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)
    }
}

#write.table(outVector, file=paste0(getwd(), "/", "5.recoded.tsv"), sep="\t", col.names=F, row.names=F, quote=F)

# FIGO.CODE

#sel<-which(dset$FIGO.CODE==4)
#dset[sel, 'FIGO.CODE']=3
stats_categorical(dset$FIGO.CODE)


# FIGO.CODE.2

# sel<-which(dset$FIGO.CODE.2==6)
# dset[sel, 'FIGO.CODE.2']=5
# sel<-which(dset$FIGO.CODE.2==7)
# dset[sel, 'FIGO.CODE.2']=5
# sel<-which(dset$FIGO.CODE.2==8)
# dset[sel, 'FIGO.CODE.2']=5
stats_categorical(dset$FIGO.CODE.2)

# MIOMETRIO

# REMOVE MIOMETRIO=3(NOT EVALUATED)

sel<-which(dset$MIOMETRIO==3)
dset<-dset[-sel,]

# PARAMETRI

# REMOVED ENTIRE VARIABLE DUE TO 88% NA

dset<-dset[,!colnames(dset)=="PARAMETRI"]

# ISTOTIPO.CODE

# (0,1,2,3)->0
# 15 REMOVED
# 8,9,10,11->3

sel<-which(dset$ISTOTIPO.CODE==0 | dset$ISTOTIPO.CODE==1 | dset$ISTOTIPO.CODE==2 | dset$ISTOTIPO.CODE==3)
dset[sel, 'ISTOTIPO.CODE']<-0
sel<-which(dset$ISTOTIPO.CODE==4)
dset[sel, 'ISTOTIPO.CODE']<-1
sel<-which(dset$ISTOTIPO.CODE==5)
dset[sel, 'ISTOTIPO.CODE']<-1
sel<-which(dset$ISTOTIPO.CODE==6)
dset[sel, 'ISTOTIPO.CODE']<-1
sel<-which(dset$ISTOTIPO.CODE==7)
dset[sel, 'ISTOTIPO.CODE']<-2
sel<-which(dset$ISTOTIPO.CODE==8 | dset$ISTOTIPO.CODE==9 | dset$ISTOTIPO.CODE==10 | dset$ISTOTIPO.CODE==11)
dset[sel, 'ISTOTIPO.CODE']<-3
sel<-which(dset$ISTOTIPO.CODE==15)
dset<-dset[-sel, ]

# LINFADENECTOMIA
# 3->1

sel<-which(dset$LINFADENECTOMIA==3)
dset[sel, 'LINFADENECTOMIA']<-0

# BMI
# 85 cases are NA
# 0. <25
# 1. >=25 & <30
# 2. >=30 & <35
# 3. >=35 & <40
# 4. >=40

sel<-which(dset$BMI<18.5)
dset[sel, 'BMI']=0
sel<-which(dset$BMI>=18.5 & dset$BMI<25)
dset[sel, 'BMI']=0
sel<-which(dset$BMI>=25 & dset$BMI<30)
dset[sel, 'BMI']=1
sel<-which(dset$BMI>=30 & dset$BMI<35)
dset[sel, 'BMI']=2
sel<-which(dset$BMI>=35 & dset$BMI<40)
dset[sel, 'BMI']=3
sel<-which(dset$BMI>=40)
dset[sel, 'BMI']=4

# ETAMENOPA
# 0. <50
# 1. >=50 & <55
# 2. >=55

sel<-which(dset$ETAMENOPA<50)
dset[sel, 'ETAMENOPA']=0
sel<-which(dset$ETAMENOPA>=50 & dset$ETAMENOPA<55)
dset[sel, 'ETAMENOPA']=1
sel<-which(dset$ETAMENOPA>=55)
dset[sel, 'ETAMENOPA']=2

# ANNIMENO
# 0. <5
# 1. >=5 & <10
# 2. >10

sel<-which(dset$ANNIMENO<5)
dset[sel, 'ANNIMENO']=0
sel<-which(dset$ANNIMENO>=5 & dset$ANNIMENO<10)
dset[sel, 'ANNIMENO']=1
sel<-which(dset$ANNIMENO>=10)
dset[sel, 'ANNIMENO']=2

# MENARCA
# 0. <12
# 1. >=12 & <15
# 2. >=15

sel<-which(dset$MENARCA<12)
dset[sel, 'MENARCA']=0
sel<-which(dset$MENARCA>=12 & dset$MENARCA<15)
dset[sel, 'MENARCA']=1
sel<-which(dset$MENARCA>=15)
dset[sel, 'MENARCA']=2

# PARITA
# 0. 0
# 1. 1
# 2. 2
# 3. 3
# 4. >3

sel<-which(dset$PARITA==0)
dset[sel, 'PARITA']=0
sel<-which(dset$PARITA==1)
dset[sel, 'PARITA']=1
sel<-which(dset$PARITA==2)
dset[sel, 'PARITA']=2
sel<-which(dset$PARITA==3)
dset[sel, 'PARITA']=3
sel<-which(dset$PARITA>3)
dset[sel, 'PARITA']=4

# GRAVIDANZE
# 0. 0
# 1. 1
# 2. 2
# 3. 3
# 4. >3

sel<-which(dset$GRAVIDANZE==0)
dset[sel, 'GRAVIDANZE']=0
sel<-which(dset$GRAVIDANZE==1)
dset[sel, 'GRAVIDANZE']=1
sel<-which(dset$GRAVIDANZE==2)
dset[sel, 'GRAVIDANZE']=2
sel<-which(dset$GRAVIDANZE==3)
dset[sel, 'GRAVIDANZE']=3
sel<-which(dset$GRAVIDANZE>3)
dset[sel, 'GRAVIDANZE']=4

# CHEMIO.CODE
# 0. 0
# 1. >0

sel<-which(dset$CHEMIO.CODE==0)
dset[sel, 'CHEMIO.CODE']=0
sel<-which(dset$CHEMIO.CODE>0)
dset[sel, 'CHEMIO.CODE']=1

# RXTERAP
# 0. 0
# 1. >0

sel<-which(dset$RXTERAP==0)
dset[sel, 'RXTERAP']=0
sel<-which(dset$RXTERAP>0)
dset[sel, 'RXTERAP']=1

# ENDO
# 0. 0
# 1. >0

sel<-which(dset$ENDO==0)
dset[sel, 'ENDO']=0
sel<-which(dset$ENDO>0)
dset[sel, 'ENDO']=1

# TERAPIA RECIDIVA 1
# 0. 0
# 1. 1 & 2 & 3
# 2. 4 & 7
# 3. 5 & 8 & 9
# 4. 6

sel<-which(dset$TERAPIA.RECIDIVA.1==1 | dset$TERAPIA.RECIDIVA.1==2 | dset$TERAPIA.RECIDIVA.1==3)
dset[sel, 'TERAPIA.RECIDIVA.1']=1
sel<-which(dset$TERAPIA.RECIDIVA.1==4 | dset$TERAPIA.RECIDIVA.1==7)
dset[sel, 'TERAPIA.RECIDIVA.1']=2
sel<-which(dset$TERAPIA.RECIDIVA.1==5 | dset$TERAPIA.RECIDIVA.1==8 | dset$TERAPIA.RECIDIVA.1==9)
dset[sel, 'TERAPIA.RECIDIVA.1']=3
sel<-which(dset$TERAPIA.RECIDIVA.1==6)
dset[sel, 'TERAPIA.RECIDIVA.1']=4

# ETA
# 0. <50
# 1. >=50 & <60
# 2. >=60 & <70
# 3. >=70

sel<-which(dset$ETA<50)
dset[sel, 'ETA']=0
sel<-which(dset$ETA>=50 & dset$ETA<60)
dset[sel, 'ETA']=1
sel<-which(dset$ETA>=60 & dset$ETA<70)
dset[sel, 'ETA']=2
sel<-which(dset$ETA>=70)
dset[sel, 'ETA']=3

# CHIR.SIGLA
# 0. DG/SALV
# 1. PIVER1
# 2. PIVER2
# 3. PIVER3

sel<-which(dset$CHIR.SIGLA=="PIVER1 ")
dset[sel, 'CHIR.SIGLA']="PIVER1"
sel<-which(dset$CHIR.SIGLA=="PIVER2 ")
dset[sel, 'CHIR.SIGLA']="PIVER2"

sel<-which(dset$CHIR.SIGLA=="DG/SALV")
dset[sel, 'CHIR.SIGLA']=0
sel<-which(dset$CHIR.SIGLA=="PIVER1")
dset[sel, 'CHIR.SIGLA']=1
sel<-which(dset$CHIR.SIGLA=="PIVER2")
dset[sel, 'CHIR.SIGLA']=2
sel<-which(dset$CHIR.SIGLA=="PIVER3")
dset[sel, 'CHIR.SIGLA']=3
dset$CHIR.SIGLA<-as.numeric(dset$CHIR.SIGLA)

# DATA RECID
# Removed 0 entries

sel<-which(dset$DATA.RECID1==0)
dset[sel, 'DATA.RECID1']=NA

#------------------------------------------------------------------------------------------------------------------------
# REMOVING VARS NOT INCLUDED IN VARS LIST
#------------------------------------------------------------------------------------------------------------------------

varNames<-read.delim(paste0(getwd(), "/vars.tsv"), header = F, sep="\t")

sel<-which(colnames(dset) %in% varNames[,1])
cset<-dset[, sel] # clinical vars
tset<-dset[, c(11, 53, 5, 6)] # time vars

aset<-cbind(tset, cset) # analysis set

#------------------------------------------------------------------------------------------------------------------------
# SURVIVAL ANALYSIS
#------------------------------------------------------------------------------------------------------------------------

  
# OVERALL SURVIVAL DIAGNOSIS

oset<-aset

STATUS<-rep(0, dim(oset)[1])
STATUS[which(oset$DATAMORTE!="")]<-1

sel<-which(oset$DATACONTR=="")
oset[sel, 'DATACONTR']<-oset[sel, 'DATAMORTE']

FUM<-calc_follow_up_monthes(oset$DATA.DIAGNOSI, oset$DATACONTR, NA)

oset<-cbind(oset, FUM, STATUS)
c<-dim(oset)[2]
plot_surv(oset[,5:c])

# EVENT FREE SURVIVAL DIAGNOSIS

oset<-aset
sel<-which(oset$DATAMORTE=="")
oset<-oset[sel,]

STATUS<-rep(0, dim(oset)[1])
STATUS[which(!is.na(oset$DATA.RECID1))]<-1

sel<-which(oset$DATA.RECID1!="")
oset[sel, 'DATACONTR']<-oset[sel, 'DATA.RECID1']

FUM<-calc_follow_up_monthes(oset$DATA.DIAGNOSI, oset$DATACONTR, NA)

oset<-cbind(oset, FUM, STATUS)
c<-dim(oset)[2]
plot_surv(oset[,5:c])

# OVERALL SURVIVAL RELAPSE TREATED 

oset<-aset
sel<-which(oset$DATA.RECID1!="")
oset<-oset[sel,]

sel<-which(oset$DATA.RECID1=="07/03/1996")
oset<-oset[-sel,]
sel<-which(oset$DATA.RECID1=="03/03/2010")
oset<-oset[-sel,]

STATUS<-rep(0, dim(oset)[1])
STATUS[which(!is.na(oset$DATAMORTE))]<-1

sel<-which(oset$DATACONTR=="")
oset[sel, 'DATACONTR']<-oset[sel, 'DATAMORTE']

FUM<-calc_follow_up_monthes(oset$DATA.RECID1, oset$DATACONTR, NA)

oset<-cbind(oset, FUM, STATUS)
c<-dim(oset)[2]
plot_surv(oset[,5:c])

#-------------------------------------------------------------------------------------------------------------------
# MULTICOLLINEARITY
#-------------------------------------------------------------------------------------------------------------------

jpeg("corplot_pp_v2.jpeg", width = 100, height = 100, units = 'mm', res = 600)
m <- cor(cset, method = "spearman")
corrplot(m, method="color", na.label=" ", tl.cex=0.4, tl.offset=0.1, tl.col="grey20",
         cl.pos="r", cl.cex=0.5, 
         addgrid.col="white", title="CORRELATION", mar=c(0,1,1,1)) #type="lower"
dev.off()

#-------------------------------------------------------------------------------------------------------------------
# COX REGRESSION
#-------------------------------------------------------------------------------------------------------------------

# OVERALL SURVIVAL

oset<-aset

STATUS<-rep(0, dim(oset)[1])
STATUS[which(oset$DATAMORTE!="")]<-1

sel<-which(oset$DATACONTR=="")
oset[sel, 'DATACONTR']<-oset[sel, 'DATAMORTE']

FUM<-calc_follow_up_monthes(oset$DATA.DIAGNOSI, oset$DATACONTR, NA)

oset<-cbind(oset, FUM, STATUS)
c<-dim(oset)[2]

coxmodel<-data.frame()

for(i in 5:(c-2)){
   
  coxmodel<-rbind(coxmodel, summary(coxph(Surv(oset$FUM, oset$STATUS)~oset[,i]))$coefficients)
}

row.names(coxmodel)<-names(oset[,5:(c-2)])

write.table(coxmodel, file=paste0(getwd(), "/", "cox_os_v2.tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)

# coxmodel<-coxph(Surv(FUM, STATUS)~ETA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~FIGO.CODE, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~FIGO.CODE.2, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~ANNESSI, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~TIPO.CA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~ISTOTIPO.CODE, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~GRADING, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~CHIR.SIGLA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~LINFADENECTOMIA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~LINFONODI, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~CARCINOSI, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~BMI, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~MENOPAUSA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~ANNIMENO, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~MENARCA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~PARITA, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~GRAVIDANZE, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~ADIUV, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~CHEMIO.CODE, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~RXTERAP, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~ENDO, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~SEDE.SIGLA1, data=oset)
# coxmodel<-coxph(Surv(FUM, STATUS)~SINT, data=oset)

# EVENT FREE SURVIVAL DIAGNOSIS

oset<-aset
sel<-which(oset$DATAMORTE=="")
oset<-oset[sel,]

STATUS<-rep(0, dim(oset)[1])
STATUS[which(!is.na(oset$DATA.RECID1))]<-1

sel<-which(oset$DATA.RECID1!="")
oset[sel, 'DATACONTR']<-oset[sel, 'DATA.RECID1']

FUM<-calc_follow_up_monthes(oset$DATA.DIAGNOSI, oset$DATACONTR, NA)

oset<-cbind(oset, FUM, STATUS)
c<-dim(oset)[2]

coxmodel<-data.frame()

for(i in 5:(c-2)){
  
  coxmodel<-rbind(coxmodel, summary(coxph(Surv(oset$FUM, oset$STATUS)~oset[,i]))$coefficients)
}

row.names(coxmodel)<-names(oset[,5:(c-2)])

write.table(coxmodel, file=paste0(getwd(), "/", "cox_efs_v2.tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)

# SURVIVAL AFTER RECURRENCE 

oset<-aset
sel<-which(oset$DATA.RECID1!="")
oset<-oset[sel,]

sel<-which(oset$DATA.RECID1=="07/03/1996")
oset<-oset[-sel,]
sel<-which(oset$DATA.RECID1=="03/03/2010")
oset<-oset[-sel,]

STATUS<-rep(0, dim(oset)[1])
STATUS[which(!is.na(oset$DATAMORTE))]<-1

sel<-which(oset$DATACONTR=="")
oset[sel, 'DATACONTR']<-oset[sel, 'DATAMORTE']

FUM<-calc_follow_up_monthes(oset$DATA.RECID1, oset$DATACONTR, NA)

oset<-cbind(oset, FUM, STATUS)
c<-dim(oset)[2]

coxmodel<-data.frame()

for(i in 5:(c-2)){
  
  coxmodel<-rbind(coxmodel, summary(coxph(Surv(oset$FUM, oset$STATUS)~oset[,i]))$coefficients)
}

row.names(coxmodel)<-names(oset[,5:(c-2)])

write.table(coxmodel, file=paste0(getwd(), "/", "cox_sar_v2.tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)


#-------------------------------------------------------------------------------------------------------------------
# LOGISTIC REGRESSION (RECURRENCE)
#-------------------------------------------------------------------------------------------------------------------

oset<-aset

STATUS<-rep(0, dim(oset)[1])
STATUS[which(oset$DATA.RECID1!="")]<-1

oset<-cbind(oset, STATUS)
c<-dim(oset)[2]

logmodel<-data.frame()

for(i in 5:(c-1)){
  
  logmodel<-rbind(logmodel, summary(glm(STATUS~oset[,i], family=binomial(link=logit)))$coefficients[2, 1:4])
}

colnames(logmodel)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
row.names(logmodel)<-names(oset[,5:(c-1)])

write.table(logmodel, file=paste0(getwd(), "/", "logistic_regression_recurrence_v2.tsv"), sep="\t", dec=".", row.names=T, col.names=NA, quote=F)


