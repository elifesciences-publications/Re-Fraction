
# feature visualization
featureBoxPlot <- function(fraction, numFraction) {
   dev.new(width=9, height=7)
   par(mfrow=c(2,2), mai=c(0.45,0.65,0.45,0.25), cex=0.75)

   boxplot(fraction$fraction.weight, col=rainbow(numFraction+1), names=paste("F", 1:numFraction, sep=""), xlab="Fraction", ylab="Mass (kDa)", main="(a)")
   boxplot(fraction$fraction.length, col=rainbow(numFraction+1), names=paste("F", 1:numFraction, sep=""), ylab="Length (aa)", main="(b)")
   boxplot(lapply(fraction$fraction.tryptic, log2), col=rainbow(numFraction+1), names=paste("F", 1:numFraction, sep=""), ylab="Log2(Number of Peptides)", main="(c)")
   boxplot(fraction$fraction.pI, col=rainbow(numFraction+1), names=paste("F", 1:numFraction, sep=""), ylab="Isoelectric Point", main="(d)")
}

# peptide visualization
peptideBarPlot <- function (peptide.assignment.fraction, peptide.assignment.original, prefix) {

   peptide.slice.unique.idx <- lapply(strsplit(peptide.assignment.fraction, ";"), length) == 1;
   peptide.original.unique.idx <- lapply(strsplit(peptide.assignment.original, ";"), length) == 1;

   unique.count <- matrix(nrow=2, ncol=2)
   unique.count[1, 1] <- sum(peptide.original.unique.idx);
   unique.count[1, 2] <- sum(peptide.slice.unique.idx);
   unique.count[2, 1] <- length(peptide.original.unique.idx) - sum(peptide.original.unique.idx);
   unique.count[2, 2] <- length(peptide.slice.unique.idx) - sum(peptide.slice.unique.idx);

   # plot absolute value
   dev.new(width=3.5, height=5)
   barplot(unique.count, ylim=c(0,20000), ylab="Number of Peptides", names=c("Original", "Re-Fraction"), col=c("darkblue","red"), legend=c("Unique Peptide", "Shared Peptide"), horiz = F, border=NA)
   text(0.7, unique.count[1, 1]/2, unique.count[1, 1],cex = 1, col="white")
   text(0.7, unique.count[1, 1] + unique.count[2, 1]/2, unique.count[2, 1], cex = 1, col="white")
   text(1.9, unique.count[1, 2]/2, unique.count[1, 2],cex = 1, col="white")
   text(1.9, unique.count[1, 2] + unique.count[2, 2]/2, unique.count[2, 2], cex = 1, col="white")
 
   # the ratio of unique peptide assignment 
   dev.new(width=3.5, height=5)
   peptide.assignment.ratio <- c(sum(peptide.original.unique.idx) / length(peptide.original.unique.idx), sum(peptide.slice.unique.idx) / length(peptide.slice.unique.idx)) 
   peptide.assignment.ratio <- round(peptide.assignment.ratio * 100); 
   barplot(peptide.assignment.ratio, names=c("Original", "Re-Fraction"), ylab="Percentage of Unique Peptide (%)", col=c("darkblue", "red"), border=NA)
   text(0.7, peptide.assignment.ratio[1]/2, paste(round(peptide.assignment.ratio[1]), "%", sep = ""), cex=1, col="white")
   text(1.9, peptide.assignment.ratio[2]/2, paste(round(peptide.assignment.ratio[2]), "%", sep = ""), cex=1, col="white")
}

# protein visualization
proteinBarPlot <- function (protein.original, protein.fraction) {
   dev.new(width=3.5, height=5)
   barplot(c(length(protein.original), length(protein.fraction)), names=c("Original", "Re-Fraction"), ylab="Number of Deterministic Protein Identifications", col=c("darkblue", "red"), border=NA)
   text(0.7, length(protein.original)/2, length(protein.original), cex=1, col="white")
   text(1.9, length(protein.fraction)/2, length(protein.fraction), cex=1, col="white")
}