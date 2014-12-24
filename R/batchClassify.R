# refraction peptide assignment

proteinClassification <- function (peptide, svm_model,
protein.weight,
protein.length,
protein.tryptic,
protein.pI) {

   all.proteins <- unique(unlist(strsplit(as.character(peptide[,"Proteins"]), ";")))
   all.proteins <- all.proteins[-grep("CON__", all.proteins)]
   all.proteins <- all.proteins[-grep("REV__", all.proteins)]

   all.proteins.df <- data.frame(cbind(protein.weight[all.proteins], protein.length[all.proteins], protein.tryptic[all.proteins], protein.pI[all.proteins]))
   rownames(all.proteins.df) <- all.proteins 
 
   protein.regressed <- as.numeric(predict(svm_model, all.proteins.df));
   names(protein.regressed) <- all.proteins 



   result <- list();
   result$regressed <- protein.regressed;
   return(result);
}

