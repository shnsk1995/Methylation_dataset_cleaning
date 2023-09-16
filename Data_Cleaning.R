#Load the required packages
library("GEOquery")
library(readxl)
library(httr)
library(stringr)

#Read the datasets information
methylation_data<- read_excel("Human_Methylation_Datasets.xlsx")


#for loop to fetch the GEO accesion number and download the raw data file and extract the metadata
for (i in 1:nrow(methylation_data)){
  
  dir.create(methylation_data$`Accession Number`[i])
  url <- methylation_data$`ftp link`[i]
  dest <- paste0(getwd(),"/",methylation_data$`Accession Number`[i],"/matrix.txt.gz")
  download.file(url,dest)
  file <- paste0(methylation_data$`Accession Number`[i],"/matrix.txt.gz")
  gse <- getGEO(filename = file,GSEMatrix=TRUE)
  write.csv(gse@phenoData@data,file = paste0(methylation_data$`Accession Number`[i],"/","metadata.csv"))
}


#Removing the datasets that have problem in downloading
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE42861",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE51032",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE210843",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE131989",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE104210",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE115278",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE114135",]

#Datasets with some problematic or wrong data
#GSE185445 - mother age and father age - ignored
#GSE213478 - Sex in numericals and age in groups - lower age saved

#Code to download a particular dataset manually
GSE <- "GSE42861"
dir.create(GSE)
file <- paste0(GSE,"/matrix.txt.gz")
gse <- getGEO(filename = file,GSEMatrix=TRUE)
write.csv(gse@phenoData@data,file = paste0(GSE,"/","metadata.csv"))


#Cleaning each dataset manually
#Saving the row number of a particular dataset to read the information for cleaning
k <- 41

#Read the data
data <- read.csv(paste0(methylation_data$`Accession Number`[k],'/metadata.csv'))

#Extract the sample_id and age
cleaned_data <- data.frame(data$X,data$characteristics_ch1.4)

#Assign column names
colnames(cleaned_data) <- c("sample_id","age")

#Extract only numerical values for age
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$age[i]
  value <- str_extract(value, "\\d+")
  cleaned_data$age[i] <- floor(as.numeric(value))
}

#Extract gender data
cleaned_data$gender <- data$characteristics_ch1.3

#Clean the gender data for uniformity
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$gender[i]
  value <- sub("Sex: ", "", value)
  #cleaned_data$gender[i] <- tolower(value)
  if(tolower(value) == "m"){
    cleaned_data$gender[i] <- "male"
  }else if (tolower(value) == "f"){
    cleaned_data$gender[i] <- "female"
  }
}

#Change unknown values to NA in any column
for(i in 1:nrow(cleaned_data)){
   if(cleaned_data$gender[i] == "unknown"){
     cleaned_data$gender[i] <- NA
   }
}

#Extract the type of data
cleaned_data$type <- data$type

#Extract the disease_state data
cleaned_data$disease_state <- data$characteristics_ch1

#Clean the disease_state data
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$disease_state[i]
  value <- sub("diagnosis: ", "", value)
  cleaned_data$disease_state[i] <- value
}

#Write the cleaned data to a csv file
write.csv(cleaned_data,file = paste0(methylation_data$`Accession Number`[k],"/","cleaned_metadata.csv"),row.names = FALSE)
