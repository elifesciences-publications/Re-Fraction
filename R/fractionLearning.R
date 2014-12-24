# refraction peptide assignment

fractionLearning <- function (peptide, numSlice, proteinRegressed) {

   #####################
   # adjusting coding  #
   # peptide.assignment <- fractionLearning(peptide.table, numFraction, allProtein$regressed)
   #peptide <- peptide.table
   #numSlice <- numFraction
   #proteinRegressed <- allProtein$regressed
   
   
   peptide.slices <- peptide[, grep("^Slice.\\d", colnames(peptide))];

   # peptide original assignment
   peptide.original.assignment <- as.character(peptide[, "Proteins"]);
   names(peptide.original.assignment) <- peptide[,"id"];

   # improve peptide assignment based on protein weight information
   peptide.assignment.range1 <- list();
   peptide.unique.range1 <- c();
   peptide.predict <- list();

   count <- 0;
   stp <- 1000;
   recorder <- 0;

   for (i in 1:nrow(peptide)) {
      # a step recorder for printing
      recorder <- recorder + 1;

      prots <- as.character(peptide[i, "Proteins"]);
      pIds <- unlist(strsplit(prots, ";"))
   
      # if the assignment is to multiple proteins
      selected.pIds.range1 <- c();
      
      for (j in 1:length(pIds)) {
         # if the peptide is assigned to a conteminated protein, keep the assignment   
         if (regexpr("CON__", pIds[j]) > 0) {
            selected.pIds.range1 <- c(selected.pIds.range1, pIds[j]);
         } else if (regexpr("REV__", pIds[j]) > 0) {
            # if the peptide is assigned to reverse protein, remove the assignment (do nothing)
         } else {
            # SVM prediction instance count
            count <- count + 1;
            
            # using the pre-classified results from batch prediction to update the currnet assignment
            preds <- proteinRegressed[pIds[j]]
            
            peptide.predict[[count]] <- paste(i, pIds[j], preds, sep=";"); 
            preds <- round(preds);
            
            # determine if the predicted range contain spectra
            range1 <- preds;
            if (range1 < 1) {
               range1 <- 1;
            } 
            if (range1 > numSlice) {
               range1 <- numSlice;
            }
            
            if(sum(peptide.slices[i, range1], na.rm=T) > 0) {
               selected.pIds.range1 <- c(selected.pIds.range1, pIds[j]);
            }
         }
      }
      
      # if we successfully assigned the peptide using range 1
      if (length(selected.pIds.range1) > 0) {
         peptide.assignment.range1[i] <- paste(selected.pIds.range1, collapse = ";");
         if (length(selected.pIds.range1) == 1){
            peptide.unique.range1 <- c(peptide.unique.range1, selected.pIds.range1);
         }
      } else { # else we use the original assignment;
         peptide.assignment.range1[i] <- paste(pIds, collapse = ";");
      }
   
      if(recorder == stp) {
         print(i);            
         recorder <- 0;
      }       
   }

   peptide.assignment.range1 <- unlist(peptide.assignment.range1)
   names(peptide.assignment.range1) <- gsub(" ", "", peptide[, "id"])
   
   # combining ReFraction updated peptide assignment to the peptide table
   peptide.ReFraction.table <- cbind(peptide, ReFraction=peptide.assignment.range1)

   result <- list();
   result$predict <- peptide.predict;
   result$original <- peptide.original.assignment;
   result$range1 <- peptide.assignment.range1;
   result$extra.unique.range1 <- peptide.unique.range1;
   result$peptide.ReFraction.table <- peptide.ReFraction.table
   return(result);
}
