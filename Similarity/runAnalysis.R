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
annodir <- DIRS[grepl("Annotation",DIRS)]
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
    onto = rough.fixNA(ONTO=onto, ANNO.LAB=models[m], NEW.ANNO.LAB=sprintf("%s.%s",models[m],"median"),MED=TRUE)
    onto = level.fixNA(ONTO=onto, LEVEL.LAB="levels", ANNO.LAB=models[m], NEW.ANNO.LAB=sprintf("%s.%s",models[m],"levels"),MED=FALSE)
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
PRECAL[[1]] = precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC.anc, ANNO.SET=1 )[[1]]
PRECAL[[2]] = precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC.anc, ANNO.SET=2 )[[1]]
PRECAL[[3]] = precal.sim( ONTO=onto, HP.SET=geneModel.HP, IC=IC.anc, ANNO.SET=3 )[[1]]

## TESTS (uncomment line 'test1 = ...' etc to run) ##

## 1) calculate the similarity between LIT and DDD gene models, using precalculate MICA scores... this will be fast
#test1 = pheno.sim.precal(MODEL.A=DDDset[[1]],MODEL.B=LITset[[1]], PRECAL=PRECAL, IC=IC.anc, ANNO.SET.A=1, ANNO.SET.B=2)

## 2) otherwise, we can calculate the similarity between a DDD and LIT gene models given the
## HP terms... this will be slower
#test2 = pheno.sim( ONTO=onto, MODEL.A=DDDset[[1]], MODEL.B=DDDset[[1]], IC=IC.anc, ANNO.SET.A=1, ANNO.SET.B=1 )

## 3) calculate similarity between two LIT gene models not using annotation data
#test3 = pheno.sim.precal(MODEL.A=LITset[[1]],MODEL.B=LITset[[2]], PRECAL=PRECAL, IC=IC.anc, ANNO.SET.A=3, ANNO.SET.B=3)

##---- DONE ------------------------------------

## calculate the similarity between all LITvDDD gene models, using precalculate MICA scores...
## this will be fast
sim.max = matrix(NA,ncol=length(DDDset),nrow=length(LITset))
colnames(sim.max) = names(DDDset)
rownames(sim.max) = names(LITset)

sim.avg = matrix(NA,ncol=length(DDDset),nrow=length(LITset))
colnames(sim.avg) = names(DDDset)
rownames(sim.avg) = names(LITset)

jsd = matrix(NA,ncol=length(DDDset),nrow=length(LITset))
colnames(jsd) = names(DDDset)
rownames(jsd) = names(LITset)


for( i in 1:length(LITset) ){
    for( j in 1:length(DDDset) ){
        res = pheno.sim.precal(MODEL.A=LITset[[i]],MODEL.B=DDDset[[j]], PRECAL=PRECAL, IC=IC.anc, ANNO.SET.A=2, ANNO.SET.B=1)
        sim.max[i,j] = res$sim.max
        sim.avg[i,j] = res$sim.avg
        jsd[i,j]     = res$jsd
    }
}
##---- DONE ------------------------------------

## calculate the similarity between all LITvDDD gene models, not using annotation data, using precalculate MICA scores...
sim.max2 = matrix(NA,ncol=length(DDDset),nrow=length(LITset))
colnames(sim.max2) = names(DDDset)
rownames(sim.max2) = names(LITset)

sim.avg2 = matrix(NA,ncol=length(DDDset),nrow=length(LITset))
colnames(sim.avg2) = names(DDDset)
rownames(sim.avg2) = names(LITset)

for( i in 1:length(LITset) ){
    for( j in 1:length(DDDset) ){
        res = pheno.sim.precal(MODEL.A=LITset[[i]],MODEL.B=DDDset[[j]], PRECAL=PRECAL, IC=IC.anc, ANNO.SET.A=3, ANNO.SET.B=3)
        sim.max2[i,j] = res$sim.max
        sim.avg2[i,j] = res$sim.avg
    }
}
##---- DONE ------------------------------------

