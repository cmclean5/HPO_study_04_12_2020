#clean R space
rm(list=ls())

#Current version of study
version="04_12_2020"

#change 'myDir' to where package is downloaded too.
myDir   <- "/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/"

#Set the absolute path to this working directory
mainDir <- sprintf("%s/HPO_study_%s",myDir,version)

#Get Path to all top-level directories in folder
DIRS <- list.dirs(mainDir,recursive=F)

#---dataDIR, where to read PPI files
dataDIR <- DIRS[grepl("networks",DIRS)]

#---subDIR for loading and storing PPI Graphs
subDIR <- list.files(path=dataDIR)
subDIR <- subDIR[grep(".csv",subDIR)]
subDIR <- unlist(strsplit(subDIR,".csv"))

#---Location for randomisation files
myRndmDir <- "/disk/scratch/WORK/DATA/"
rndDIR    <- vector(length=1)
rndDIR[1] <- sprintf("%s/HPO_study_%s/Random_Model_Pairs/",myRndmDir,version)

#HPOperGene=FALSE
HPOperGene=TRUE

if( HPOperGene ){

Ngenes=77

#---model short names
pid  <- seq(0,(Ngenes-1),1)

#---model short names
ptype = c(sprintf("LIT:%gx",pid),sprintf("DDD:%gx",pid))

#---model long names
phenl = c(sprintf("literature_model:%gx",pid),sprintf("ddd_model:%gx",pid))

pid   = seq(0,(length(ptype)-1),1)

} else {
    
#---model short names
pid  <- vector(length=2);
pid[1]   = 1
pid[2]   = 2

#---model short names
ptype  <- vector(length=2);
ptype[1]   = "LIT";
ptype[2]   = "DDD";

#---model long names
phenl  <- vector(length=2);
phenl[1]   = "literature_model"
phenl[2]   = "ddd_model"

}
    
#---For Bonferroni correction
alpha <- vector(length=3)
alpha[1] <- 0.05
alpha[2] <- 0.01
alpha[3] <- 0.001

stars    <- vector(length=3)
stars[1] <- "*"
stars[2] <- "**"
stars[3] <- "***"

COLLAPSE <- vector(length=2)
COLLAPSE[1] <- ";"
COLLAPSE[2] <- "&"
c=1

#--- Set 
#--- Pre-load Graph of interest, stored in file 'graphs.csv'
pramFILES <- DIRS[grepl("parameterFiles",DIRS)]
Graph <- read.table(sprintf("%s/graphs.csv",pramFILES),header=F,sep="\t",quote="")
S     <- as.vector(Graph[which(as.vector(Graph[,1]) == 1)[1],2])
S     <- match(S,subDIR)
#S     <- grep(S,subDIR)
StudyName <- as.character(Graph[which(as.vector(Graph[,1]) == 1)[1],3])

#---Set clip of interest
#--- clip study
clipStudy <- matrix(NA,ncol=2,nrow=2)
clipStudy[1,1] = "forClip_4920_HPterms"; clipStudy[1,2] = "fullClip";
clipStudy[2,1] = "forClip_2548_HPterms"; clipStudy[2,2] = "litANDdddClip";

## clip choice
V=1

#---Set role-up study
roleUpStudy <- "roleUp"


#---All annotation types
Anno     <- vector(length=5)
Anno[1] = "gene"
Anno[2] = "disease.name"
Anno[3] = "literature_model"
Anno[4] = "ddd_model"
Anno[5] = "OV"


#--- TAX ID
#TAXid    <- vector(length=3)
#TAXid[1] <- "10090" #Mouse
#TAXid[2] <- "9606"  #Human
#TAXid[3] <- "7227"  #Fly

#GOonto    <- vector(length=3)
#GOonto[1] <- "MF"
#GOonto[2] <- "BP"
#GOonto[3] <- "CC"

#---WIDTH and HEIGHT for plots (in px)
WIDTH=480
HEIGHT=480

#---WIDTH and HEIGHT for plots (in inches)
WIDTHin=6.4
HEIGHTin=6.4

#---WIDTH and HEIGHT for plots (in cm)
WIDTHcm=12.7
HEIGHTcm=12.7

#---WIDTH and HEIGHT for plots (in mm)
WIDTHmm=127
HEIGHTmm=127

## for numerical stability   
SMALL = 1e-100

# define constants
EPSILON = .Machine$double.eps
#----

#set required R libraries
#source("http://www.bioconductor.org/biocLite.R")

library(igraph);
library(lattice);
library(ggplot2);
library(cowplot);
library(purrr);
library(stringr);
library(Matrix);   ## sparse matrix manipulation
##library(rdetools); ## radial basis function kernel, convert distance to similarity
library(kernlab);  ## radial basis function kernel, convert distance to similarity
library(wordspace); ## sparse matrix distance calculation
library(proxyC);    ## sparse matrix distance calculation
#library(MASS);
#library(gtable);
#library(grid);
#library(ggrepel);
#library(scales);
#library(clusterCons); ##commented out only because working on ubuntu home machine!
#library(plyr);
#library(VennDiagram);
#library(Vennerable);
#library(biomaRt);
#library(latex2exp);
#library(knitr);
#library(poweRlaw);
#library(WriteXLS);
#library(gdata);
library(methods);
#library(ggpubr);
#library(DBI);
#library(aricode);
library(reshape2);
#library(svglite);

#library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(org.Dm.eg.db)
#library(topGO)

#set default R options 
options(stringsAsFactors=F)

#printSessionInfo <- TRUE
printSessionInfo <- FALSE

if( printSessionInfo ){
    print(sessionInfo())
}

#printPackageInfo <- TRUE
printPackageInfo <- FALSE

if( printPackageInfo ){

cat("//----------------------------------------------\n")
cat("// Package Name       : Human synaptomsome  \n")
cat("// Package Version    : ", version, "\n")
cat("// Package Description: Structural, Functional and Disease analysis of\n") 
cat("//                    : of model PSP and Presynaptic Human PPI networks\n")
cat("//                    : for the Human Brain Project SGA1. \n")
cat("// Date               : 2019\n")
cat("// Author             : Colin D Mclean <Colin.D.Mclean@ed.ac.uk>\n")
cat("// Copyright (C) 2016 Colin Mclean \n")
cat("//----------------------------------------------\n")
cat("//      Package Description\n")
cat("//----------------------------------------------\n")
cat("// The package is split into serveral subpackages to build, cluster and analysis \n")
cat("// Human PSP and Presynaptic PPI network models:\n")
cat(" script(s)      : buildPPInetworks.R, filterPPInetworks.R, Rclustering.R\n")
cat(" script(s)      : Cppclustering.R, recluster.R, validateClustering.R\n")    
cat("//----------------------------------------------\n")
cat("//      subPackages \n")
cat("//----------------------------------------------\n")
cat("[1] Annotations : a store of the raw annotations (disease and functional) files, and\n")
cat("                  add this annotation data to the graphs.\n")
cat(" script(s)      : addAnnotation2Graphs.R\n")
cat("[2] Clustering  : a store of the clustering files for each algorithm and each PPI network model.\n")
cat("[3] Consensus   : build and plot the consensus matrices for each algorithm, and calculate the Bridgeness values.\n")
cat(" script(s)      : buildConsensusCDF.R, plotConsensusCDF.R, Bridgeness.R\n ")
cat("[4] DiseasePairs: calculate the co-occurance of disease pairs on a PPI network model.\n")
cat(" script(s)      : calObservedDiseasePairs.R, calRANDOMDiseasePairs.R, plots.R \n")
cat(" script(s)      : mergeRESULTS.R, groupDiseaseResults.R \n")
cat("[5] EnrichmentPackage: Stand alone package to reading in clustering and annotation files to calculate the network and clusteral enrichment.\n")
cat(" script(s)      : submitClustEnrch.sh, submitOverlapEnrch.sh, enrichmentAnalysis.R, sigmoidFit.R \n")
cat("[6] EntropyRate : calculate the graph entropy.\n")
cat(" script(s)      : PertubationEntropy.R\n")    
cat("[7] GO          : GO enrichment for given PPI network model using topGO.\n")
cat(" script(s)      : enrichmentGO.R\n")    
cat("[8] Graphs      : store of graphical PPI network models in .gml format.\n")
cat("[9] parameterFiles: set of parameter files to pass into the EnrichmentPackage.\n")
cat("[10] POWERlawFIT: Calculate vertex centrality measures, and Powerlaw fit to log of network degree distribution.\n")    
cat(" script(s)      : Centrality.R, ranCentrality.R\n")    
cat("//----------------------------------------------\n")
cat("//      GNU General Public Licenses v3 \n")
cat("//----------------------------------------------\n")
cat("// This program is free software: you can redistribute it and/or modify it \n")
cat("// under the terms of the GNU General Public License as published by the \n")
cat("// Free Software Foundation, either version 3 of the License, or (at your \n")
cat("// option) any later version.\n")
cat("//\n")
cat("// This program is distributed in the hope that it will be useful, but \n")
cat("// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY \n")
cat("// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n")
cat("//\n")
cat("// You should have received a copy of the GNU General Public License \n")
cat("// (GNU_GPL_v3)  along with this program.  If not, see <http://www.gnu.org/licenses/>.\n")
cat("//----------------------------------------------\n")
cat("//      Funding Acknowledgement\n")
cat("//----------------------------------------------\n")
cat("// This open source software code was developed in part or in whole in the Human Brain Project, funded from the European Union’s Horizon 2020\n")
cat("// Framework Programme for Research and Innovation under the Specific Grant Agreement No. 720270 (Human Brain Project SGA1).\n")
cat("//-----------------------------------------------\n")
cat("/////////////////////////////////////////////////\n")
cat("//-----------------------------------------------\n")

}
    
#---Print Graph of interest
cat("\n")
cat("\n")
cat("***********\n")
cat(" Graph/file is: ", subDIR[S], " set in 'graphs.csv' \n")
cat("***********\n")
cat("\n")
cat("\n")
