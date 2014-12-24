\name{applyReFraction}
\alias{applyReFraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function of Re-Fraction
}
\description{
This is the main function that apply Re-Fraction algorithm for deterministic protein identifications.
}
\usage{
ReFraction(protein.database, peptide.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein.database}{
The protein database that was used for peptide identification. This database should have been processed to 
extract protein properties including: protein id, protein mass, protein length, number theoretic tryptic peptides, and protein pI. See the supplied example database for detail.
There is a function for extracting protein properties from fasta file. See the example at the bottom.
}
  \item{peptide.table}{
ReFraction accepts peptide table that including fraction information and protein assignments. For fraction information, peptide table
should contain columns with names start with "Slice" end with a number and separated by ".".  An example of such a column would be
"Slice.1". For protein assignments, there should be a column with name "Proteins" and contains one or more protein identifiers that the peptide is assigned to. 
The identifier should be separated by ";". Column names are case sensitive. The detailed structure of such a peptide table could be found in the example dataset. The
detail of this dataset could be found in refernce [2].
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{determine.proteins}{This is a list of deterministic protein identifications after
applying Re-Fraction. Proteins in the list are ranked by their posterior error probability (PEP).
The smaller the PEP the more confident the protein identified.
}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
[1] Pengyi Yang, Sean J. Humphrey, Daniel J. Fazakerley, Matthew J. Prior, Guang Yang, David E. James, and Jean Yee-Hwa Yang, 
Re-Fraction: a machine learning approach for deterministic identification of protein homologs and splice variants in large-scale MS-based proteomics, 
Journal of Proteome Research, http://dx.doi.org/10.1021/pr300072j

[2] Matthew J. Prior, Mark Larance, Robert T. Lawrence, Jamie Soul, Sean Humphrey, James Burchfield, Carol Kistler, Jonathon R. Davey, Penelope J. La-Borde, Michael Buckley, Hiroshi Kanazawa, Robert G. Parton, Michael Guilhaus, and David E. James,
Quantitative proteomic analysis of the adipocyte plasma membrane,
Journal of Proteome Research, 2011, 10(11), 4970-4982
}
\author{
Pengyi Yang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load the example protein database that used for peptide identification
data("ipiMOUSEv385")
# load the example peptide table
data("peptideTable")

# to have a view of how extracted protein database should look like
protein.database[1:10,]

# to have a view of how peptide table should look like. There are extra columns than what are required by Re-Fraction. The required columns are: "id", "Proteins", "PEP", and all the columns start wtih "Slice".
peptide.table[1:10,]

# perform Re-Fraction algorithm for deterministic protein identifications
ReFraction.results <- applyReFraction(protein.database, peptide.table, fdr.cutoff=0.01)

# display deterministic protein identifications; assigned peptide ids match those in input peptide.table
ReFraction.results$deterministic.protein.table[1:15,]

# display peptide table; the last column contains the refined peptide assignment after applying ReFraction
ReFraction.results$peptide.ReFraction.table[1:15,]

# The above is an example run on the dataset used in reference [1].
# For processing your dataset, you need to firstly extracting protien properties from the database you used in peptide identification.
# Use the following function for this:
?extractDatabase

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
