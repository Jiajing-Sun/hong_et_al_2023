rm(list=ls())
set.seed(123456789 ) 
folder <- c("/Users/sunjiajing/Codes-empirical-analysis" ) 
setwd(folder) 
source("source-empirical.R") 
library(readxl)
library(forecast)
library(timeDate)
library(tseries)
library(sandwich)
library(xtable)

###################################################################  
# load data
temp.dat <- read_xlsx ("data.xlsx") 
dat <- temp.dat [ setdiff ((1: dim(temp.dat)[1]),which (temp.dat $weekdays=="Saturday" | temp.dat $weekdays=="Sunday") ) , ]

Dow_Jones <- dat$Dow_Jones_Index_Industrial_Average
SP<- dat$`S&P/TSX_Composite_Index` 
FTSE_100 <- dat$Index_FTSE_100
CAC_40 <- dat$France_Index_CAC_40
DAX<- dat$Germany_Index_DAX 

# rate of returns
Dow_Jones.return <- 	price2return(	dat$Dow_Jones_Index_Industrial_Average)
SP.return<- 	price2return(	dat$`S&P/TSX_Composite_Index` )
FTSE_100.return <- 	price2return(	dat$Index_FTSE_100)
CAC_40.return <- 	price2return(	dat$France_Index_CAC_40)
DAX.return<- 	price2return(	dat$Germany_Index_DAX )

#summary statistics
temp <- rbind(summary.statistics(Dow_Jones.return),
              summary.statistics(SP.return), 
              summary.statistics(FTSE_100.return),
              summary.statistics(CAC_40.return),
              summary.statistics(DAX.return)) 
################################################################### 
folder.new <- paste0(folder, "/output")
setwd(folder.new )

row.names <- c("Dow Jones" , "SP Composite", "FTSE 100" ,"CAC 40" , "DAX"  )
col.names <- c("min", "max", "average", "sd", "skewness", "kurtosis", "ADF", "ADF-P", "normality", "normality-P")
rownames(temp)<-row.names; colnames(temp) <- col.names 
#write.csv(temp, "summary_statistics_date_20221217.csv") 

temp_filename<- c("summary_statistics_empirical.tex")
table1 <- xtable(temp ,caption= paste0("Summary statistics of rates of returns for five major stock indices"), digits=4) 
print(table1 , file = temp_filename) 

table1


