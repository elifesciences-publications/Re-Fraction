

applyReFraction <- function(protein.database, peptide.table, fdr.cutoff=0.01) {

############### load features from protein database 
rownames(protein.database) <- protein.database[,"pId"];

protein.weight <- as.numeric(protein.database[, "weight.kDa."]);
protein.length <- as.numeric(protein.database[, "length"]);
protein.tryptic <- as.numeric(protein.database[, "trypticCount"]);
protein.pI <- as.numeric(protein.database[, "pI"]);

names(protein.weight) <- protein.database[,"pId"];
names(protein.length) <- protein.database[,"pId"];
names(protein.tryptic) <- protein.database[,"pId"];
names(protein.pI) <- protein.database[,"pId"];
print("Step 1: database loaded");
   
   

############## step 3: load the peptide identification results from PM datasets
rownames(peptide.table) <- peptide.table[,"id"];
fractionIdxs <- grep("^Slice.\\d", colnames(peptide.table));
numFraction <- length(fractionIdxs);
print("Step 2: peptide identification table loaded");


############### step 4: creating modeling data

# substep 1: select training PSMs
#source("sampleSelector.R")
train.psm <- sampleSelector(peptide.table, fractionIdxs)
print("Step 3: training example selected");


# substep 2: construct modeling dataset
#source("trainDataConstractor.R")
train.fraction <- trainDataConstractor(train.psm$pId, train.psm$fraction, numFraction, 
protein.weight, 
protein.length, 
protein.tryptic, 
protein.pI)
print("Step 4: modeling dataset created");

# visualization 1: boxplot of features by fractions
#source("visualization.R")
featureBoxPlot(train.fraction, numFraction);

############### step 5: creating the model and classifying all the identified proteins
#source("trainSVM.R")
svm <- trainSVM(train.fraction, numFraction);

#### step 5.1 classify all the identified proteins to gel fractions
#source("batchClassify.R")
allProtein <- proteinClassification(peptide.table, svm$svm_model,
protein.weight, 
protein.length, 
protein.tryptic, 
protein.pI)
print("Step 5: classification applied");   

############### step 6: applying Re-Fraction   
#source("fractionLearning.R")

print("applying Re-Fraction on peptide dataset");
peptide.assignment <- fractionLearning(peptide.table, numFraction, allProtein$regressed)
print("Step 6: Re-Fraction applied");

# visualization 2: peptide plot
peptideBarPlot(peptide.assignment$range1, peptide.assignment$original);

############### step 7: unique protein identifications
#source("determineProteinID.R")
protein.assignment <- determineProteinID(peptide.assignment$peptide.ReFraction.table, peptide.assignment$range1, peptide.assignment$original, fdr.cutoff)

#visualization 3: protein plot
proteinBarPlot(protein.assignment$original, protein.assignment$fraction);

# Re-Fraction results
result <- list()
result$peptide.ReFraction.table <- peptide.assignment$peptide.ReFraction.table
result$deterministic.protein.table <- protein.assignment$deterministic.protein.table

return(result)

}