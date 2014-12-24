# deterministic protein identification and filtering

determineProteinID <- function (peptide, peptide.slice.assignment, peptide.original.assignment, fdr.cutoff=0.01, fdr.stp=20) {
   
   # indexing the unique peptide assignment
   peptide.slice.unique.idx <- lapply(strsplit(peptide.slice.assignment, ";"), length) == 1;
   peptide.original.unique.idx <- lapply(strsplit(peptide.original.assignment, ";"), length) == 1;

   # unique protein identifications
   peptide.slice.unique <- peptide.slice.assignment[peptide.slice.unique.idx]
   peptide.original.unique <- peptide.original.assignment[peptide.original.unique.idx]

   # "protein.slice.unique.info" is used to capture peptide id from the peptide table
   protein.slice.unique.info <- split(names(peptide.slice.unique), peptide.slice.unique)
   protein.slice.unique <- unique(peptide.slice.unique);
   protein.original.unique <- unique(peptide.original.unique);

   protein.slice.peptideCount <- summary(as.factor(peptide.slice.unique), maxsum = Inf)
   protein.original.peptideCount <- summary(as.factor(peptide.original.unique), maxsum = Inf)

   # using a two peptide identification for protein filtering
   protein.slice.multiPeptide <- names(protein.slice.peptideCount)[(protein.slice.peptideCount > 1)]
   protein.original.multiPeptide <- names(protein.original.peptideCount)[(protein.original.peptideCount > 1)]

   # get the peptide PEP 
   peptide.slice.pep <- cbind(peptide[peptide.slice.unique.idx, "PEP"], peptide.slice.unique)
   rownames(peptide.slice.pep) <- names(peptide.slice.unique)
   colnames(peptide.slice.pep) <- c("PEP", "pId");
   peptide.original.pep <- cbind(peptide[peptide.original.unique.idx, "PEP"], peptide.original.unique)
   rownames(peptide.original.pep) <- names(peptide.original.unique)
   colnames(peptide.original.pep) <- c("PEP", "pId");

   # sort the peptide by PEP
   peptide.slice.pep.sort <- peptide.slice.pep[order(peptide.slice.pep[,"PEP"]),]
   peptide.original.pep.sort <- peptide.original.pep[order(peptide.original.pep[,"PEP"]),]

   # calculate the PEP for protein
   protein.slice.PEP <- c();
   for (i in 1:length(protein.slice.multiPeptide)){
      prot <- peptide.slice.pep.sort[which(peptide.slice.pep.sort[,"pId"] == protein.slice.multiPeptide[i]), "PEP"];
      protein.slice.PEP <- c(protein.slice.PEP, prod(as.numeric(prot)));
   }
   names(protein.slice.PEP) <- protein.slice.multiPeptide;
   protein.slice.sort <- protein.slice.PEP[order(protein.slice.PEP)];

   protein.original.PEP <- c();
   for (i in 1:length(protein.original.multiPeptide)){
      prot <- peptide.original.pep.sort[which(peptide.original.pep.sort[,"pId"] == protein.original.multiPeptide[i]), "PEP"];
      protein.original.PEP <- c(protein.original.PEP, prod(as.numeric(prot)));
   }
   names(protein.original.PEP) <- protein.original.multiPeptide;
   protein.original.sort <- protein.original.PEP[order(protein.original.PEP)];


   # filtering the protein identification by FDR cutoff
   index <- 0;
   count <- 0;
   stp <- fdr.stp;
   decoy <- 0;
   for (i in 1:length(protein.slice.sort)) {
      count <- count + 1;
      if(regexpr("REV__", names(protein.slice.sort)[i]) > 0) {
         decoy <- decoy + 1;
      }
   
      if (stp == count){
         count <- 0;
         fdr <- 2 * decoy / i;
         index <- i;
         if (fdr > fdr.cutoff){
            break();
         }
      }
   }

   protein.slice.filtered <- protein.slice.sort[1:index]
   
   # filtering the protein identification by FDR 0.01
   index <- 0;
   count <- 0;
   stp <- fdr.stp;
   decoy <- 0;
   for (i in 1:length(protein.original.sort)) {
      count <- count + 1;
      if(regexpr("REV__", names(protein.original.sort)[i]) > 0) {
         decoy <- decoy + 1;
      }
   
      if (stp == count){
         count <- 0;
         fdr <- 2 * decoy / i;
         index <- i;
         if (fdr > fdr.cutoff){
            break();
         }
      }
   }

   protein.original.filtered <- protein.original.sort[1:index]

   # filtering the reverse and contenminations
   protein.slice.filtered <- protein.slice.filtered[grep("REV", invert=T, names(protein.slice.filtered))]
   protein.slice.filtered <- protein.slice.filtered[grep("CON", invert=T, names(protein.slice.filtered))]

   protein.original.filtered <- protein.original.filtered[grep("REV", invert=T, names(protein.original.filtered))]
   protein.original.filtered <- protein.original.filtered[grep("CON", invert=T, names(protein.original.filtered))]

   
   # create the deterministic protein table as final result
   deterministic.protein.table <- cbind(protein=names(protein.slice.filtered), PEP=protein.slice.filtered, peptides=sapply(protein.slice.unique.info[names(protein.slice.filtered)], function(x){paste(x, collapse=";")}))
   
   result <- list();
   result$deterministic.protein.table <- deterministic.protein.table
   result$fraction <- protein.slice.filtered;
   result$original <- protein.original.filtered;
   result$fraction.pepSort <- peptide.slice.pep.sort;
   result$original.pepSort <- peptide.original.pep.sort;
   
   return(result);
}