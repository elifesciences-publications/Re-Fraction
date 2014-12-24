
##############################################
#    creating the training data              #
##############################################

trainDataConstractor <- function(peptide.unique.fractions.pId, peptide.unique.fractions, numfraction, 
protein.weight, 
protein.length, 
protein.tryptic, 
protein.pI) {

   # grouping protein features to individual fractions
   fraction.weight <- list();
   fraction.length <- list();
   fraction.tryptic <- list();
   fraction.pI <- list();

   length(fraction.weight) <- numfraction;
   length(fraction.length) <- numfraction;
   length(fraction.tryptic) <- numfraction;
   length(fraction.pI) <- numfraction;


   stp <- 100;
   recorder <- 0;
   for (i in 1:length(peptide.unique.fractions.pId)){
      recorder <- recorder + 1;
      pId <- peptide.unique.fractions.pId[i];

      sId <- NULL;
      opt.sId <- as.numeric(which.max(peptide.unique.fractions[i,]));
      r <- sum(!is.na(peptide.unique.fractions[i, ]));
      
       
      if (r < 3) { # if the range is 1 or 2, select the optimal "maximum"
         sId <- opt.sId;
      } else if (max(peptide.unique.fractions[i,], na.rm=T) != median(as.numeric(peptide.unique.fractions[i,]), na.rm=T)) {
         sId <- opt.sId;
      } else {
         sId <- opt.sId + 1;
      }
   
      # feature 1 weight
      tmp <- fraction.weight[[sId]];
      tmp <- c(tmp, protein.weight[pId]);
      fraction.weight[[sId]] <- tmp;
   
      # feature 2 length
      tmp <- fraction.length[[sId]];
      tmp <- c(tmp, protein.length[pId]);
      fraction.length[[sId]] <- tmp;
   
      # feature 3 tryptic count
      tmp <- fraction.tryptic[[sId]];
      tmp <- c(tmp, protein.tryptic[pId]);
      fraction.tryptic[[sId]] <- tmp;
      
      # feature 4 pI
      tmp <- fraction.pI[[sId]];
      tmp <- c(tmp, protein.pI[pId]);
      fraction.pI[[sId]] <- tmp;
   
      if(recorder == stp) {
         print(i);            
         recorder <- 0;
      }       
   }


   # training data cleansing
   fraction <- list();
   for (i in 1:numfraction) {
      fraction[[i]] <- summary(as.factor(names(fraction.weight[[i]])), maxsum = Inf);
   }

   pId <- unique(peptide.unique.fractions.pId)

   for (i in 1:length(pId)) {
      assignment <- c();
      for(j in 1:numfraction){
         assignment <- c(assignment, fraction[[j]][pId[i]]);
      }

      count <- sum(!is.na(assignment))

      if (count > 1) {
         best.sId <- which.max(assignment);
         del.idx <- which(!is.na(assignment));
         for (j in 1:length(del.idx)){
            idx <- del.idx[j];
      
            if (idx != best.sId) {
               del <- which(names(fraction.weight[[idx]]) == pId[i]);
               fraction.weight[[idx]] <- fraction.weight[[idx]][-c(del)];
               fraction.length[[idx]] <- fraction.length[[idx]][-c(del)];
               fraction.tryptic[[idx]] <- fraction.tryptic[[idx]][-c(del)];
               fraction.pI[[idx]] <- fraction.pI[[idx]][-c(del)];        
            }
         }
      }        
   }

   train.raw <- list();
   train.raw$fraction.weight <- fraction.weight;
   train.raw$fraction.length <- fraction.length;
   train.raw$fraction.tryptic <- fraction.tryptic;
   train.raw$fraction.pI <- fraction.pI;
   return(train.raw);
}