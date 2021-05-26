source('../setUp.R')

library("SteinerNet")
##library("pcSteiner")

##load functions
source("SteinerFunctions.R")

#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("DDD_dataset",DIRS)]

#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")    
}

plotdir <- "PLOTS"

##---Read-in orginal annotation file
ANNO   <- read.delim("../Annotations/hpo_lists_per_disease.csv",sep=",",header=T)
GN     <- ANNO[match(ANNO[,1],pid[1:length(ANNO[,2])]),2]
GN_PLT <- sprintf("%s_%d",GN,pid[1:Ngenes])
#DIS  <- ANNO[match(ANNO[,1],pid[1:length(ANNO[,2])]),3]

#---load corresponding graph which was used to build the consensus matrices from 
gg    <- igraph::read.graph(sprintf("%s/%s.gml",OUT[1],subDIR[S]),format="gml")

N = length(V(gg))
E = length(E(gg))

Lmax = as.numeric(max(V(gg)$level))
Targ = 7

oo1 = roleUpTerms(GG=gg,level.max=Lmax,level.target=Lmax)
oo2 = roleUpTerms(GG=gg,level.max=Lmax,level.target=Targ)

xx1 = graph.frac.annoAtLevels(GG=gg)
xx2 = frac.annoAtLevels(MAP=oo2)

write.table(xx1, "HPOlevels.csv",sep="\t",col.names=T,row.names=F,quote=F)
write.table(xx2, sprintf("HPOlevels_atL%d.csv",Targ),sep="\t",col.names=T,row.names=F,quote=F)

xx3 = map.annoAtLevels.perGene(MAP=oo1, level.max=Lmax)

gp1 = plot.annoAtLevels(DF=xx3$LIT    ,XLAB=GN_PLT,YLAB="LIT per gene")
plot.save(GPLOT=gp1,PLOTDIR=plotdir,FILENAME="HPOlevels_LIT_perGene",
          pHEIGHT=2*HEIGHT,pWIDTH=WIDTH,UNITS="px")


gp2 = plot.annoAtLevels(DF=xx3$DDD    ,XLAB=GN_PLT,YLAB="DDD per gene")
plot.save(GPLOT=gp2,PLOTDIR=plotdir,FILENAME="HPOlevels_DDD_perGene",
          pHEIGHT=2*HEIGHT,pWIDTH=WIDTH,UNITS="px")


gp3 = plot.annoAtLevels(DF=xx3$LIT_DDD,XLAB=GN_PLT,YLAB="LITvDDD per gene")
plot.save(GPLOT=gp3,PLOTDIR=plotdir,FILENAME="HPOlevels_LITvDDD_perGene",
          pHEIGHT=2*HEIGHT,pWIDTH=WIDTH,UNITS="px")

xx4 = map.annoAtLevels.perGene(MAP=oo2, level.max=Lmax)

gp4 = plot.annoAtLevels(DF=xx4$LIT    ,XLAB=GN_PLT,YLAB="LIT per gene")
plot.save(GPLOT=gp4,PLOTDIR=plotdir,FILENAME=sprintf("HPOlevels_LIT_atL%d_perGene",Targ),
          pHEIGHT=2*HEIGHT,pWIDTH=WIDTH,UNITS="px")

gp5 = plot.annoAtLevels(DF=xx4$DDD    ,XLAB=GN_PLT,YLAB="DDD per gene")
plot.save(GPLOT=gp5,PLOTDIR=plotdir,FILENAME=sprintf("HPOlevels_DDD_atL%d_perGene",Targ),
          pHEIGHT=2*HEIGHT,pWIDTH=WIDTH,UNITS="px")

gp6 = plot.annoAtLevels(DF=xx4$LIT_DDD,XLAB=GN_PLT,YLAB="LITvDDD per gene")
plot.save(GPLOT=gp6,PLOTDIR=plotdir,FILENAME=sprintf("HPOlevels_LITvDDD_atL%d_perGene",Targ),
          pHEIGHT=2*HEIGHT,pWIDTH=WIDTH,UNITS="px")


xx5 = map.annoAtLevels(MAP=oo1,level.max=Lmax)
CN  = levels(factor(oo1[,5]))
CN  = CN[CN!=""]
pd1 = plot.annoAtLevels(DF=xx5$MODEL,XLAB=CN,YLAB="MODEL")
plot.save(GPLOT=pd1,PLOTDIR=plotdir,FILENAME="HPOlevels_Model",
          pHEIGHT=HEIGHT,pWIDTH=2*WIDTH,UNITS="px")


xx6 = map.annoAtLevels(MAP=oo2,level.max=Lmax)
CN  = levels(factor(oo2[,5]))
CN  = CN[CN!=""]
pd2 = plot.annoAtLevels(DF=xx6$MODEL,XLAB=CN,YLAB="MODEL")
plot.save(GPLOT=pd2,PLOTDIR=plotdir,FILENAME=sprintf("HPOlevels_atL%d_Model",Targ),
          pHEIGHT=HEIGHT,pWIDTH=2*WIDTH,UNITS="px")


##---Write role-up graph to level 'Targ' to to file
oo2 = as.data.frame(oo2)
gg2 = gg
gg2 = removeVertexTerm(gg2,"OV")
gg2 = set.vertex.attribute(gg2,"OV",V(gg2),oo2$OV)
gg2 = removeVertexTerm(gg2,"OVG")
gg2 = set.vertex.attribute(gg2,"OVG",V(gg2),oo2$OVG)

#igraph::write.graph(gg, sprintf("%s/roleUp/%s_l%d.gml",OUT[1],subDIR[S],Targ), format="gml")
