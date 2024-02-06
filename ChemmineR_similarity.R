### Similarity calculation using the ChemmineR maximum common substructure function
# Löffler et al. 2024, DOI: 
setwd("C:/Users/pllo0001/OneDrive - Sveriges lantbruksuniversitet/Skrivbordet/PhD/Review/Viewpoint")

library(readr) #data import
library(tidyverse) #pipe operator
library(data.table) #data handling
library(ChemmineR) #MCS computation

# data import and treatment --------------------
orig.data <- read_csv("C:/Users/pllo0001/OneDrive - Sveriges lantbruksuniversitet/Skrivbordet/PhD/Review/Viewpoint/PubChem_all_transformations_wExtraInfo.csv") %>% setDT()
data <- orig.data
data <- data[-c(4178,5652,5654)] #Exclude wrong SDFs from TP File
sdfset_parents <- read.SDFset("similarity_parents.sdf")
sdfset_parents <- sdfset_parents[-c(4178,5652,5654)]
sdfset_tps <- read.SDFset("similarity_tps.sdf")
valid.parents <- validSDF(sdfset_parents)
valid.tps <- validSDF(sdfset_tps)
valid <- as.data.frame(cbind(valid.parents, valid.tps))

false_rows <- which(!valid$valid.parents | !valid$valid.tps, arr.ind = TRUE)
print(false_rows)

ttest <- data.frame(MCS = numeric(7241)) #create empty dataframe for loop. This might take some time.
for(i in 1:7241){
  mcstest <- fmcs(sdfset_parents[i], sdfset_tps[i], au=2, bu=1)
  ttest[i,] <- mcstest@stats["Tanimoto_Coefficient"]
}

results <- cbind(data,ttest) %>% setDT()
#write.csv(results, "similarity_mcsTanimoto.csv")

### the file similarity_mcsTanimoto.csv is on Github so you can also start from here ---------------
results <- read_csv("similarity_mcsTanimoto.csv") %>% setDT()

ggplot(data = results, aes(x = MCS))+
  geom_histogram(aes(y =..density..), color = "black", fill = "white")+
  geom_vline(xintercept = 0.95, linetype = 1, color = "red", size = 1)+
  labs(x = "MCS tanimoto coefficient", y = "frequency count")+
  geom_density(alpha=.1, fill="#FF6666") +
  theme_bw()

results[MCS>=0.95, Similarity := as.character("similar")]
results[MCS<0.95, Similarity := as.character("dissimilar")]

table(results$Similarity)

print(sprintf("%d compounds were classified as similar, which is equivalent to %.2f%%.", nrow(results[results$MCS >= 0.95, ]), 100 * nrow(results[results$MCS >= 0.95, ]) / nrow(results)))
