## download the necessary packages for data acquisition and analysis

library(BiocManager)

library(Biostrings)

library(DECIPHER)

library(rentrez)

library(seqinr)

library(tidyverse)

library(pROC)

library(ggplot2)

library(randomForest)

##First, the entrez search function from the rentrez package is used to explore the COI sequences for both the class and the order to verify the number of sequences. web history is being used as there is a large amount of sequences in  each dataset.

Actinopterygii.sequences.raw <- entrez_search(db = "nuccore", term = "		
Actinopterygii[ORGN] AND COI[GENE] AND 500:1000[SLEN]", retmax = 10000, use_history = TRUE)

Gymnotiformes.sequences.raw <- entrez_search(db = "nuccore", term = "		
Gymnotiformes[ORGN] AND COI[GENE] AND 500:1000[SLEN]", retmax = 980, use_history = TRUE)

##Next, we use a for loop to obtain the sequences from the NCBI database and write those sequences to a fasta folder. although there is more than 10,000 sequences, guidelines for the assignment suggested a maximum of 10,000 sequences when assessing k-mer frequencies down the line. For Gymnotiformes, there are only 990 samples, but the sequence amount was set to a 1000 for conveniance. This code is adapted from Rentrez tutorial page, under the web history section (https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#use-a-web_history-object)

for( seq_start in seq(1,10000,50)){
  Fasta_Actinopterygii <- entrez_fetch(db="nuccore", web_history=Actinopterygii.sequences.raw$web_history,
                                       rettype="fasta", retmax=50, retstart=seq_start)
  cat(Fasta_Actinopterygii, file="Fasta_Actinopterygii_COI.fasta", append=TRUE)
  cat(seq_start+49, "sequences downloaded\r")
}


for( seq_start in seq(1,1000,50)){
  Fasta_Gymnotiformes <- entrez_fetch(db="nuccore", web_history=Gymnotiformes.sequences.raw$web_history,
                                      rettype="fasta", retmax=50, retstart=seq_start)
  cat(Fasta_Gymnotiformes, file="Fasta_Gymnotiformes_COI.fasta", append=TRUE)
  cat(seq_start+49, "sequences downloaded\r")
}

##Next, the fasta files produced are read in as DNA string sets from the biostrings package which allows for easier manipulation and visualization of the data

Actinopterygii_Sequences <- readDNAStringSet("Fasta_Actinopterygii_COI.fasta")

Gymnotiformes_Sequences <- readDNAStringSet("Fasta_Gymnotiformes_COI.fasta")

##Here, the nucletoides are oriented all in the same direction utilizing the orient nucleotides function from the DECIPHER package to maintain consistency

Actinopterygii_Sequences <- OrientNucleotides(Actinopterygii_Sequences)

Gymnotiformes_Sequences <- OrientNucleotides(Gymnotiformes_Sequences)

##Next, we specify some paramters to convert our current data into a data frame for each set of nucleotides. 2 columns are specified for gene source and sequence, then both are concatenated using the base R data.frame function.

Gene_source = names(Actinopterygii_Sequences)
sequence = paste(Actinopterygii_Sequences)
dfActinopterygii_Sequences <- data.frame(Gene_source, sequence)

Gene_source = names(Gymnotiformes_Sequences)
sequence = paste(Gymnotiformes_Sequences)
dfGymnotiformes_Sequences <- data.frame(Gene_source, sequence)

##remove objects that are no longer needed to clean up

rm(Actinopterygii.sequences.raw, Actinopterygii_Sequences, Gymnotiformes.sequences.raw, Gymnotiformes_Sequences, test)

##Next, since all Gymnotiformes have some electrogenesis capability, a column termed "electrogenesis" is added to the data set with all values set to "capable". the same is done for the Actinopterygii dataset, except all values are set to "uncapable".

dfGymnotiformes_Sequences <- dfGymnotiformes_Sequences %>% mutate(electrognesis = "Capable")

dfActinopterygii_Sequences <- dfActinopterygii_Sequences %>% mutate(electrognesis = "Not Capable")

##Finally, we combine the 2 dataframes together, with the "capable" values first, then remove any duplicates that may re appear in the datasets associated with the initial dataset because some Gymnotiformes may still be in the 10,000 sequences of Actinopterygii.

dfFused <- rbind(dfGymnotiformes_Sequences, dfActinopterygii_Sequences)

dfElectrogensis_Sequences <- dfFused[!duplicated(dfFused$Gene_source),]

##visualize the data in a table to see how many sequences are assoacited with the capability of electrogenesis.

table(dfElectrogensis_Sequences$electrognesis)

##For the analysis, most of the source code is adapted from script 9 from class as demonstrated by Professor Adamowicz. The generation of the k-mer frequencies was specifically adopted from script 9. Other features from script 9 are adopted with some modifications.(Adamowicz, 2021)


## First we separate the elements of the data frame and transform only the sequence column to a DNA string so that the base pairs in the string can be read in individually. 
dfElectrogensis_Sequences <- as.data.frame(dfElectrogensis_Sequences)
dfElectrogensis_Sequences$sequence <- DNAStringSet(dfElectrogensis_Sequences$sequence)

##A check to see the class of the sequences after changing them into DNA strings.
class(dfElectrogensis_Sequences$sequence)

##Here we count the number of base pairs in each sequence, and append new columns that highlight the count of each base pair in each column
dfElectrogensis_Sequences <- cbind(dfElectrogensis_Sequences, as.data.frame(letterFrequency(dfElectrogensis_Sequences$sequence, letters = c("A", "C","G", "T"))))

##Next, to find the proportions of each base pair, a new column is appended that highlights the base pair count divided by the sum of all base pair counts
dfElectrogensis_Sequences$Aprop <- (dfElectrogensis_Sequences$A) / (dfElectrogensis_Sequences$A + dfElectrogensis_Sequences$T + dfElectrogensis_Sequences$C + dfElectrogensis_Sequences$G)

dfElectrogensis_Sequences$Tprop <- (dfElectrogensis_Sequences$T) / (dfElectrogensis_Sequences$A + dfElectrogensis_Sequences$T + dfElectrogensis_Sequences$C + dfElectrogensis_Sequences$G)

dfElectrogensis_Sequences$Gprop <- (dfElectrogensis_Sequences$G) / (dfElectrogensis_Sequences$A + dfElectrogensis_Sequences$T + dfElectrogensis_Sequences$C + dfElectrogensis_Sequences$G)


##Now that the proportions are set, the dinucleotide frequency function from the bio strings package is used generate the first set of frequencies.
dfElectrogensis_Sequences <- cbind(dfElectrogensis_Sequences, as.data.frame(dinucleotideFrequency(dfElectrogensis_Sequences$sequence, as.prob = TRUE)))

##The same can be done for trinucloetide frequencies.
dfElectrogensis_Sequences <- cbind(dfElectrogensis_Sequences, as.data.frame(trinucleotideFrequency(dfElectrogensis_Sequences$sequence, as.prob = TRUE)))

##And finally, the same is done for K-mers of length 4.
dfElectrogensis_Sequences <- cbind(dfElectrogensis_Sequences, as.data.frame(oligonucleotideFrequency(x = dfElectrogensis_Sequences$sequence, width = 4, as.prob = TRUE)))

##Return the column back to a character class.
dfElectrogensis_Sequences$sequence <- as.character(dfElectrogensis_Sequences$sequence)

##visualize the electrogenesis again to see the amounts for sample number specifications
table(dfElectrogensis_Sequences$electrognesis)


##set seed to make the results reproducible. a verification dataset is created, containing 100 samples that are capable and 100 that aren't capable of electrogenesis.
set.seed(333)

dfVereification_Dataset <- dfElectrogensis_Sequences %>% group_by(electrognesis) %>% sample_n(100)


##verify that the desired sample number is obtained
table(dfValidation_Dataset$electrognesis)

##set seed to make the results reproducible. a practice dataset is created, containing 320 samples that are capable and 320 that aren't capable of electrogenesis. Samples that are in the verification dataset are excluded from practice dataset
set.seed(444)

dfPractice_Dataset <- dfElectrogensis_Sequences %>%
  filter(!Gene_source %in% dfVereification_Dataset$Gene_source) %>%
  group_by(electrognesis) %>%
  sample_n(320)


##verify that the desired sample number is obtained
table(dfPractice_Dataset$electrognesis)

##Here we create the classifier for electrogenesis using all the variables generated which include proportions, and k-mers of length 2,3 and 4. The responsible variable in this context is electrogenesis. n tree was specified as 50 because the practice data set is relatively small. 
Electrogenisis_classifier <- randomForest::randomForest(x = dfPractice_Dataset[, 8:346], y = as.factor(dfPractice_Dataset$electrognesis), ntree = 50, importance = TRUE)

##view the results of the classifier.
Electrogenisis_classifier

##error rate of 2.19%. acceptable performance.

##Now we test the model on the unseen data and see if still performs well with data not present in the practice dataset. View the predictions made
Model_Verification <- predict(Electrogenisis_classifier, dfVereification_Dataset[, c(2, 8:346)])
Model_Verification

##compare and visualize the results using a table.
table(observed = dfVereification_Dataset$electrognesis, predicted = Model_Verification)

##only 2 mistakes... pretty good!!


## For visualizations, first we plot the model to see how much trees are needed before the models error rate drops.
plot(Electrogenisis_classifier)

##Next, we can use the variable importance plot to verify which k-mers or proportions were the most important for sample classification.
varImpPlot(Electrogenisis_classifier)


##For the final visualization, an ROC or receiver operating characteristic curve for each set of classifier variables will be plotted to judge the performace of the model under different contexts

#First, we need a numerical representation of electrogenesis. we use the classifier to predict results of the verification dataset. setting the type to "prob" will give me the probability of electrogenesis for a given sample. norm votes gives fractions for use. proximity and nodes are set to false because they are not needed.code for generating this numeric data from binary data was adapted from WV View (http://www.wvview.org/spatial_analytics/random_forest/_site/index.html)
pred_test <- predict(Electrogenisis_classifier, dfVereification_Dataset[, c(2, 8:346)], index=2, type="prob", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)

##convert the results to a data frame
pred_test <- data.frame(pred_test)

##create an ROC curve using the verification dataset as reference. the value being tested for in this context is "capable" of electrogenesis 
Model_Verification_roc <- roc(dfVereification_Dataset$electrognesis, pred_test$Capable)

##plot model to verify effectiveness of the model
plot(Model_Verification_roc)

##Next we create a classifier with all variables as before except this time, only the nucleotide proportions will be used. All other varaibles are explained previously.
Electrogenisis_classifier_proportions <- randomForest::randomForest(x = dfPractice_Dataset[, 8:10], y = as.factor(dfPractice_Dataset$electrognesis), ntree = 50, importance = TRUE)

##a numerical representation of the new classifier.
pred_test <- predict(Electrogenisis_classifier_proportions, dfVereification_Dataset[, c(2, 8:10)], index=2, type="prob", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)

##conversion into a dataframe.
pred_test <- data.frame(pred_test)

##producing an ROC
Proportions_classifier_roc <- roc(dfVereification_Dataset$electrognesis, pred_test$Capable)

##creation of the new classifier using proportions and dinucleotides
Electrogenisis_classifier_dinucleotide <- randomForest::randomForest(x = dfPractice_Dataset[, 8:26], y = as.factor(dfPractice_Dataset$electrognesis), ntree = 50, importance = TRUE)

##a numerical representation of the new classifier
pred_test <- predict(Electrogenisis_classifier_dinucleotide, dfVereification_Dataset[, c(2, 8:26)], index=2, type="prob", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)

##conversion into a dataframe.
pred_test <- data.frame(pred_test)

##producing an ROC
Dinucleotide_classifier_roc <- roc(dfVereification_Dataset$electrognesis, pred_test$Capable)

##creation of the new classifier using proportions, dinucleotides, and trinucleotides.
Electrogenisis_classifier_trinucleotide <- randomForest::randomForest(x = dfPractice_Dataset[, 8:90], y = as.factor(dfPractice_Dataset$electrognesis), ntree = 50, importance = TRUE)

##a numerical representation of the new classifier
pred_test <- predict(Electrogenisis_classifier_trinucleotide, dfVereification_Dataset[, c(2, 8:90)], index=2, type="prob", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)

##conversion into a dataframe.
pred_test <- data.frame(pred_test)

##producing an ROC
Trinucleotide_classifier_roc <- roc(dfVereification_Dataset$electrognesis, pred_test$Capable)

##Finally, all ROC curves are plotted into one plot using the ggroc function from ggplot2. the ROC curves are specified and listed. The geom_abline with intercept and slope 1 serves as reference. Each of the different curves is represented by a different color. Both x and y axis are labelled and the plots text and legend are modified. General idea of the code was adapted from WV view and some modifications were made.(http://www.wvview.org/spatial_analytics/random_forest/_site/index.html)

ggroc(list(All = Model_Verification_roc,
           Proportions_only=Proportions_classifier_roc,
           with_dinucleotide=Dinucleotide_classifier_roc,
           with_trinucleotides=Trinucleotide_classifier_roc), lwd=1)+
  geom_abline(intercept = 1, slope = 1, color = "red", linetype = "dashed", lwd=1.2)+
  scale_color_manual(labels = c("K-mer length 4", "Proportions only", "K-mer length 2","K-mer length 3"), values= c("#88D04B", "#D65076", "#EFC050","#7FFFD4"))+
  labs(x="Specificity", y="Sensitivity")+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12))+
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(size=14))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.title=element_blank())

##This completes the analysis

