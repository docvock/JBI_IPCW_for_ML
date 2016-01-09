###################################
# File: ./graphcis.R
#		Purpose: Create graphics which display the length of followup and proportion of subjects where 
#							the adverse event status is unknown. 
#   Author: David M. Vock
#   Last Modified: January 8, 2016
#   Relies On: none
#   Files Created: ./followup.pdf 
#										./followup.eps
#										./prop_unknonwn.pdf 
#										./prop_unknonwn.eps 
####################################

setwd("//hifs00.ahc.umn.edu/Data/Carlson/Portal")
dat.base <- "Data/DataLock/DataLockv2.3.1/"

# libraries used
library(data.table)
library(dplyr)

#  Read the data
HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
#  Restrict Dataset to those over 40 and with no comorbidites
HP <- subset(HP.imp, age >= 40 & Comorbidity_All == 0) 

#pdf("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\followup.pdf",
setEPS()
postscript("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\followup.eps",
	height = 4, width = 4)
hist(HP$DaysToEvent_Fram/365, breaks = seq(from = 0, to =10, by = 0.5), 
	xlab = "Years from Index Date", ylab = "Number of Subjects in Cohort", main = "", col = "grey")
hist(HP$DaysToEvent_Fram[HP$CVDEvent_Fram == 1]/365, breaks = seq(from = 0, to =10, by = 0.5), 
	xlab = "Years from Index Date", ylab = "Number of Subjects in Cohort", main = "", col = "white", add = TRUE)
dev.off()

setEPS()
postscript("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\followup_v2.eps",
	height = 4, width = 8)
#pdf("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\followup_v2.pdf",
#	height = 4, width = 8)
par(mfrow = c(1, 2))
hist(HP$DaysToEvent_Fram/365, breaks = seq(from = 0, to =10, by = 0.5), 
	xlab = "Years from Index Date", ylab = "Number of Subjects", main = "Censored", col = "grey")
hist(HP$DaysToEvent_Fram[HP$CVDEvent_Fram == 1]/365, breaks = seq(from = 0, to =10, by = 0.5), 
	xlab = "Years from Index Date", ylab = "Number of Subjects", main = "CV Event", col = "grey")
dev.off()


HP <- mutate(HP, T.use = DaysToEvent_Fram, C.use = CVDEvent_Fram)

t <- seq(from = 0, to = 8, by = 0.5)
unknown.prop <- NULL
for (i in 1:length(t)) {
unknown <- 1-prop.table(table(ifelse(HP$T.use < (t[i]*365) & HP$C.use==0, 1, 0)))[1]
unknown.prop <- c(unknown.prop, unknown)
	print(i)
}

#pdf("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\prop_unknown#.pdf",
setEPS()
postscript("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\prop_unknown.eps",
	height = 4.5, width = 4.5)
plot(t, unknown.prop, type = "b", pch = 19, ylab = "Proportion with Event Status Unknown", 
	xlab = "Years from Index Date")
dev.off()

