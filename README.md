# Bioinformatics-Tidbits
A collection of scripts and files used whilst bioinformagicianing

**BioNetwork_2.R** - A script to produce a network graph of the MCRI bioinformatics collaborations with other groups based on a google docs spreadsheet.

**DEXY.R** - A script to extract Differential Exon Usage using DexSeq from a counts matrix (rows exons, columns samples) then summarised to per gene q values.

**DIY.html** - An html file (based on a Rmarkdown Rstudio notebook) which is a tutorial on how make your own ALL classifier, in this example i base the classifier of a Lund dataset.

**Diff_Ex_Skeleton.Rmd** - A cheeky R markdown skeleton i made for a _standard_ Differential Expression RNA seq analysis, this was made in the early days so it probably isn't that advanced.

**Expression Outliers.R** - A quick script i made to explore Expression Outliers for various interesting genes in the ALL samples we though were more unusual from the MCRI cohort.

**Express Yourself** - My first ever shiny app for looking into expression outliers as a web app instead of statically as in the above script.

**Making_AllSorts.Rmd** - A script to make the AllSorts Classifier (really messy, not cleaned)

**Making_AllSorts.html** - A script to make the AllSorts Classifier (to view it, preview the html on: https://rawgit.com/Quarkins/Bioinformatics-Tidbits/master/Making_AllSorts.html

**Combined Dec2016** - A script which combines all the counts and target data of samples measured at MCRI and then classifies them and constructs RData objects used by the Express Yourself R Shiny App
