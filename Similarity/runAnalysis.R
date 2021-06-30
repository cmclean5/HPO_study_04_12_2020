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

## now read-in frequency annotation file
df   = read.delim(sprintf("%s/hpo_freq_all_combined.csv",annodir),sep=",",header=T)

X  = as.character(df[,6])
GN = sapply(X,function(x) strsplit(x,",")[[1]][1])

Anno = cbind(as.character(df[,2]),as.numeric(df[,4]),as.numeric(df[,5]),as.character(GN))
colnames(Anno) <- c(colnames(df)[c(2,4,5)],"gene")
colnames(Anno)[which(colnames(Anno)=="ddd_freq")] = models[1]
colnames(Anno)[which(colnames(Anno)=="metamap_freq")] = models[2]
Anno = as.data.frame(Anno)

## add gene model, i.e. ddd and metamap, frequency data to ontology object
for( m in 1:length(models) ){
    onto = rm.onto.annotation(ONTO=onto, LABEL=models[m])
    onto = add.onto.annotation(ONTO=onto, ANNO=Anno, HPLABEL="hpo_id", LABEL=models[m])
    onto = rough.fixNA(ONTO=onto, ANNO.LAB=models[m], NEW.ANNO.LAB=sprintf("%s.%s",models[m],"median"),type="median")
    onto = level.fixNA(ONTO=onto, LEVEL.LAB="levels", ANNO.LAB=models[m], NEW.ANNO.LAB=sprintf("%s.%s",models[m],"levels"),type="mean")
}
##------------------------

## build hpo graph object (optional)
gg  <- build.hpo.graph(ONTO=onto);
##------------------------

## build/use subontology (optional)
SUB=FALSE
#SUB=TRUE
if( SUB ){
    onto = sub.ontology(ONTO=onto,  TARGETS=clip);
    sub <- get.clip.hpo.graph(GG=gg,TARGETS=clip);    
}
##-----------------


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


## find all unique ancestor HP terms for each set
LITanc <- list()
DDDanc <- list()
for( i in 1:Ngenes ){

    if( SUB ){ hp.terms = get.sub.ancestors(onto, LITset[[i]])
    } else { hp.terms   = get.ancestors(onto, LITset[[i]]) }
    
    LITanc[[i]] = unique(hp.terms)
                   
    if( SUB ){ hp.terms = get.sub.ancestors(onto, DDDset[[i]])
    } else { hp.terms   = get.ancestors(onto, DDDset[[i]]) }
        
    DDDanc[[i]] = unique(hp.terms)
                         
    names(LITanc)[i] = names(LITset)[i]
    names(DDDanc)[i] = names(DDDset)[i]
}
#----------------------------------------

## find all unique HP terms in each set
lit.HP <- c()
ddd.HP <- c()

lit.anc.HP <- c()
ddd.anc.HP <- c()

for( i in 1:Ngenes ){
    lit.HP     = c(lit.HP,LITset[[i]])
    lit.anc.HP = c(lit.anc.HP,LITanc[[i]])
    ddd.HP     = c(ddd.HP,DDDset[[i]])
    ddd.anc.HP = c(ddd.anc.HP,DDDanc[[i]])
}

lit.HP     = unique(lit.HP)
lit.anc.HP = unique(lit.anc.HP)
ddd.HP     = unique(ddd.HP)
ddd.anc.HP = unique(ddd.anc.HP)

geneModel.HP     = unique(c(lit.HP,ddd.HP))
geneModel.anc.HP = unique(c(lit.anc.HP,ddd.anc.HP,geneModel.HP))
#----------------------------------------

## get IC for each HP term, and ancestor HP terms, in gene models...
ROOT.IC = list()
IC.anc  = list()

for( m in 1:Nmodels ){

    anno.lab = sprintf("%s.%s",models[m],"levels")

    ic = get.model.IC(ONTO=onto, HP.TERMS=geneModel.anc.HP, ROOT=ROOTS[1],
                      MODEL.LABEL=anno.lab, PRINT=FALSE, SUB=SUB);

    cat("------\n")
    
    ROOT.IC[[m]]      = ic$ROOT.IC
    IC.anc[[m]]       = ic$IC.anc

    names(ROOT.IC)[m] = anno.lab
    names(IC.anc)[m]  = anno.lab

}

## also get the IC for each gene without any annotation
anno.lab="no.anno";

ic = get.model.IC(ONTO=onto, HP.TERMS=geneModel.anc.HP, ROOT=ROOTS[1],
                      MODEL.LABEL=anno.lab, PRINT=FALSE, SUB=SUB);

ROOT.IC[[(Nmodels+1)]]      = ic$ROOT.IC
IC.anc[[(Nmodels+1)]]       = ic$IC.anc

names(ROOT.IC)[(Nmodels+1)] = anno.lab
names(IC.anc)[(Nmodels+1)]  = anno.lab

##----------------------------------------------------


## pre-calculate MICA scores between all HP terms in ddd and metamap... 
## this will take so time ~15-20mins per ANNO set.
PRECAL = list()
PRECAL[[1]] = NULL;#precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC.anc, ANNO.SET=1 )[[1]]
PRECAL[[2]] = NULL;#precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC.anc, ANNO.SET=2 )[[1]]
PRECAL[[3]] = precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC.anc, ANNO.SET=3 )[[1]]

##--- TESTS ---##
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
ddd.sim = cal.sim.matrix( MODEL.A=DDDset, MODEL.B=DDDset, PRECAL=PRECAL, IC=IC.anc,
                          ANNO.SET.A=3, ANNO.SET.B=3,
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
lit.sim = cal.sim.matrix( MODEL.A=LITset, MODEL.B=LITset, PRECAL=PRECAL, IC=IC.anc,
                          ANNO.SET.A=3, ANNO.SET.B=3,
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
Com.sim = cal.sim.matrix( MODEL.A=DDDset, MODEL.B=LITset, PRECAL=PRECAL, IC=IC.anc,
                          ANNO.SET.A=3, ANNO.SET.B=3,
                          TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def")

## set which jsd to plot 
#com.jsd = Com.sim$jsd.ic
com.jsd = Com.sim$jsd.mica
#com.jsd = Com.sim$jsd.avg

## heat-map of distance between average IC value between two gene models: no scaled
com.gp      = list()
com.gp[[1]] = heat.map(DF=com.jsd,leg.tit="distance",SCALE=NULL,fill.diag=TRUE)

## heat-map of distance between average IC value between two gene models: scaled
for( i in 1:length(ScalingFactor) ){
    com.gp[[i+1]] = heat.map(DF=com.jsd,SCALE=ScalingFactor[i],fill.diag=TRUE)
}


##---- DONE ------------------------------------

