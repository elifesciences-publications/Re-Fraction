###################################################
#   selecting PSMs for creating training data     #
###################################################

sampleSelector <- function(peptide, slice.idx) {

   # remove the reverse matches and the contemination matches
   peptide.filtered <- peptide[grep("REV", invert=T, peptide[, "Proteins"]),]
   peptide.filtered <- peptide.filtered[grep("CON", invert=T, peptide.filtered[, "Proteins"]),]

   # get the uniquely assigned protein IDs
   peptide.assignment <- as.character(peptide.filtered[, "Proteins"]);
   peptide.unique.idx <- lapply(strsplit(peptide.assignment, ";"), length) == 1;
   peptide.unique.pId <- peptide.assignment[peptide.unique.idx];
   names(peptide.unique.pId) <- rownames(peptide.filtered)[peptide.unique.idx]

   peptide.slices <- peptide.filtered[peptide.unique.idx, slice.idx];

   # get the proteins that are identified from: (1) no more than 3 slices, and (2) a continue range
   peptide.range.slices <- peptide.slices[(rowSums(!is.na(peptide.slices)) <= 3),]
   peptide.count.slices <- peptide.range.slices[apply(peptide.range.slices, 1, max, na.rm=T) > 1, ]

   selected <- c();
   for (i in 1:nrow(peptide.count.slices)) {
      idx <- which(!is.na(peptide.count.slices[i,]))
   
      if (length(idx) == 1) {
         selected <- c(selected, rownames(peptide.count.slices)[i]);
      } else {
         flag <- 0;
         for (j in 1:(length(idx) - 1)){
            if((idx[j] + 1) != idx[j + 1]) {
               flag <- 1;
            }
         }
      
         if (flag == 0) {
            selected <- c(selected, rownames(peptide.count.slices)[i]);
         }
      }
   }

   peptide.unique.slices.pId <- peptide.unique.pId[selected]
   peptide.unique.slices <- peptide.slices[selected,]

   train.psm <- list();
   train.psm$selected <- selected;
   train.psm$pId <- peptide.unique.slices.pId;
   train.psm$fraction <- peptide.unique.slices;

   return(train.psm);
}