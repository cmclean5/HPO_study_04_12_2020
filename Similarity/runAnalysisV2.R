## directory paths and load required libraries
source('../setUp.R')

## see ../setUp.R for libraries required
#library(ontologyIndex)
#library(ggplot2)
#library(igraph)

## put all our functions here and load them
source('loadFunctions.R')

## set random number seed
set.seed(1)

## directories
annodir <- DIRS[grepl("Annotations",DIRS)]
#annodir <- DIRS[grepl("Annotations",DIRS) & !grepl("safe",DIRS)]
grdir   <- DIRS[grepl("Graphs",DIRS)]
obodir  <- DIRS[grepl("OBO",DIRS)]

## load HP ontology from latest OBO file
## may first have to download it into the OBO directory, i.e.,
## 1) cd HPO_study_30_04_2021/OBO
## 2) wget http://purl.obolibrary.org/obo/hp.obo
onto <- get_ontology(sprintf("%s/hp.obo",obodir),
                     propagate_relationships = "is_a",
                     extract_tags = "minimal")
#data(hpo)

## remove all obsolete hp terms
onto <- rm.obs.terms(ONTO=onto);

## find root node in onto
ROOTS <- root.terms(onto);
##-----------------

## add clip HP terms
clip <- read.delim(sprintf("%s/clip_HPO_terms.csv",annodir),sep="\t",header=F)[[1]]
onto <- add.hp.clip(ONTO=onto,CLIP=clip)
##-----------------

## add hp term level
onto <- get.hp.levels(ONTO=onto);
##-----------------

## add hp term number of leaf nodes
onto <- get.hp.leafs(ONTO=onto);
##-----------------

## define our gene model labels
models    <- vector(length=2)
models[1] <- "ddd"
models[2] <- "metamap"

Nmodels   <- length(models);

## build/use subontology (optional)
SUB=FALSE
#SUB=TRUE


## read-in model HP term per gene annotation data
LITset = readRDS(sprintf("%s/LIT.Gene.HPterms.RDS",annodir))
DDDset = readRDS(sprintf("%s/DDD.Gene.HPterms.RDS",annodir))

genes  = names(LITset)
Ngenes = length(genes)

## remove any redundant HP terms for sets
for( i in 1:Ngenes ){
    LITset[[i]] = get.terms.in.onto(ONTO=onto,TARGET=LITset[[i]]);
    DDDset[[i]] = get.terms.in.onto(ONTO=onto,TARGET=DDDset[[i]]);
}
#----------------------------------------

## add all unique ancestor HP terms to each gene model in each set
ANCset      <- list()
ANCset[[1]] <- add.ancestor.terms(ONTO=onto,SUB=SUB,SET=LITset)
ANCset[[2]] <- add.ancestor.terms(ONTO=onto,SUB=SUB,SET=DDDset)

## get IC for each HP term, and ancestor HP terms, in gene models...
ROOT.IC = list()
IC      = list()

for( m in 1:Nmodels ){

    anno.lab = models[m];

    ic = get.model.IC.2(ONTO=onto, SET=ANCset[[m]],
                        MODEL.LABEL=anno.lab, PRINT=FALSE, SUB=SUB);

    cat("------\n")
    
    ROOT.IC[[m]]      = ic$ROOT.IC
    IC     [[m]]      = ic$IC.anc

    names(ROOT.IC)[m] = anno.lab
    names(IC)[m]      = anno.lab

}


## find all unique HP terms in each set
lit.HP <- c()
ddd.HP <- c()

for( i in 1:Ngenes ){
    lit.HP     = c(lit.HP,LITset[[i]])
    ddd.HP     = c(ddd.HP,DDDset[[i]])
}

lit.HP       = unique(lit.HP)
ddd.HP       = unique(ddd.HP)
geneModel.HP = unique(c(lit.HP,ddd.HP))
#----------------------------------------


## pre-calculate MICA scores between all HP terms in ddd and metamap... 
## this will take so time ~15-20mins per ANNO set.
PRECAL = list()
PRECAL[[1]] = precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC, ANNO.SET=1 )[[1]]
PRECAL[[2]] = precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC, ANNO.SET=2 )[[1]]

## NOTE: A metric which seems to work is:
##    1) first calculating the average IC with each gene model,
##    2) then calculating the js_distance (jsd) between these two models.

## heat map scaling factor, i.e show how similar two gene models are, for jsd.
ScalingFactor=vector(length=4)
ScalingFactor[1] = 5
ScalingFactor[2] = 50
ScalingFactor[3] = 100
ScalingFactor[4] = 200


## calculate the similarity between all DDDvDDD gene models, using precalculate MICA scores...
## this will be fast
ddd.sim = cal.sim.matrix( MODEL.A=DDDset, MODEL.B=DDDset, PRECAL=PRECAL, IC=IC,
                          ANNO.SET.A=1, ANNO.SET.B=1,
                          TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def" )

## set which jsd to plot 
#ddd.jsd = ddd.sim$jsd.ic
ddd.jsd = ddd.sim$jsd.mica

## heat-map of distance between average IC value between two gene models: no scaled
ddd.gp      = list()
ddd.gp[[1]] = heat.map(DF=ddd.jsd,leg.tit="distance",SCALE=NULL)

## heat-map of distance between average IC value between two gene models: scaled
for( i in 1:length(ScalingFactor) ){
    ddd.gp[[i+1]] = heat.map(DF=ddd.jsd,SCALE=ScalingFactor[i])
}


## calculate the similarity between all LITvLIT gene models, using precalculate MICA scores...
## this will be fast
lit.sim = cal.sim.matrix( MODEL.A=LITset, MODEL.B=LITset, PRECAL=PRECAL, IC=IC,
                          ANNO.SET.A=2, ANNO.SET.B=2,
                          TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def")

## set which jsd to plot 
#lit.jsd = lit.sim$jsd.ic
lit.jsd = lit.sim$jsd.mica

## heat-map of distance between average IC value between two gene models: no scaled
lit.gp      = list()
lit.gp[[1]] = heat.map(DF=lit.jsd,leg.tit="distance",SCALE=NULL)

## heat-map of distance between average IC value between two gene models: scaled
for( i in 1:length(ScalingFactor) ){
    lit.gp[[i+1]] = heat.map(DF=lit.jsd,SCALE=ScalingFactor[i])
}


## calculate the similarity between all DDDvLIT gene models, using precalculate MICA scores...
## this will be fast
Com.sim = cal.sim.matrix( MODEL.A=DDDset, MODEL.B=LITset, PRECAL=PRECAL, IC=IC,
                          ANNO.SET.A=1, ANNO.SET.B=2,
                         TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def")


## set which jsd to plot 
#com.jsd = Com.sim$jsd.ic
com.jsd = Com.sim$jsd.mica
#com.jsd = Com.sim$jsd.avg

## heat-map of distance between average IC value between two gene models: no scaled
com.gp      = list()
com.gp[[1]] = heat.map(DF=com.jsd,leg.tit="distance",SCALE=NULL,fill.diag=FALSE)

## heat-map of distance between average IC value between two gene models: scaled
for( i in 1:length(ScalingFactor) ){
    com.gp[[i+1]] = heat.map(DF=com.jsd,SCALE=ScalingFactor[i],fill.diag=FALSE)
}
