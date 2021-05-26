source('../setUp.R')

library("SteinerNet")
##library("pcSteiner")

##load functions
source("SteinerFunctions.R")

#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("DDD_dataset",DIRS)]

## clip choice, see value for parameter 'V' setUp.R

#---load corresponding graph which was used to build the consensus matrices from 
gg    <- igraph::read.graph(sprintf("%s/%s.gml",OUT[1],subDIR[S]),format="gml")
clipF <- igraph::read.graph(sprintf("%s/%s/HPOclip_fragmented.gml",OUT[1],clipStudy[V,2]),format="gml")

N = length(V(gg))
E = length(E(gg))

V(gg)$prizes <- rep(1,N)
E(gg)$costs  <- rep(1,E)

root_node = V(gg)$name[which(V(gg)$level==0)]
max_depth = max(as.numeric(V(gg)$level)) ## max_depth == 14

#---Read in HP clipping node - i.e. the terminal nodes of Steiner tree
# forClip_4920_HPterms.csv, 4920 unique terms between DDD_HPterms_16_12_2020, and LIT.
HPc <- read.delim(sprintf("%s/%s.csv",OUT[2],clipStudy[V,1]),sep="\t",header=F)

sTerm <- getTerminals(gg,as.vector(HPc[[1]]))

## SteinerNet
nTrees <- 100
trees  <- buildSteinerTrees(gg,sTerm,nTrees)


clips  <- list()
lccs  <- matrix(NA,nrow=nTrees,ncol=8)
for( i in 1:nTrees ){
    ##clips[[i]] = connectFragments(GG=gg,CLIPF=clipF,TREES=trees,INDX=i)
    Tres       = connectFragments(GG=gg,CLIPF=clipF,TREES=trees,INDX=i)
    clips[[i]] = connectSmallFagments(GG=gg,CLIP=Tres,TARGETS=as.vector(HPc[[1]]))
    temp       = findLCCinfo(clips[[i]])
    ov         = table(V(clips[[i]])$OV)
    if( is.vector(temp) ){        
        lccs[i,1]  = i
        lccs[i,2]  = 1
        lccs[i,3]  = length(V(clips[[i]]))
        lccs[i,4]  = length(E(clips[[i]]))
        lccs[i,5]  = as.vector(ov[1])
        lccs[i,6]  = as.vector(ov[2])
        lccs[i,7]  = as.vector(ov[3])
        lccs[i,8]  = as.vector(ov[4])
    } else {
        lccs[i,1]  = i
        lccs[i,2]  = length(temp[,1])
        lccs[i,3]  = length(V(clips[[i]]))
        lccs[i,4]  = length(E(clips[[i]]))
        lccs[i,5]  = as.vector(ov[1])
        lccs[i,6]  = as.vector(ov[2])
        lccs[i,7]  = as.vector(ov[3])
        lccs[i,8]  = as.vector(ov[4])
    }
}

## choose 'best' clip, with max 'lit&ddd' terms, and min 'other' terms.
indx = lccs[which(lccs[,8]==max(lccs[,8])),]
if( length(indx[,1]) > 1 ){
    indx2 = which(lccs[indx[,1],5]==min(lccs[indx[,1],5]))
    indx  = indx[indx2,1]
}
#----


#tree_len <- unlist(lapply(trees, function(x) length(E(x[[1]]))))
#indx     <- which(tree_len == min(tree_len))

clip <- clips[[indx[1]]]


## pcSteiner
#res2 = pcs.tree(graph=gg, terminals=sTerm, lambda=1, root=root_node, depth=max_depth,eps=1e-3,max_iter=10)

##---Write clip.gml graph to file
igraph::write.graph(clip, sprintf("%s/HPOclip.gml",OUT[1]), format="gml")

##--- Clip 2, subgraph of HPO using clips nodes.
clip2 = clip.from.subgraph(GG=gg,CLIP=clip)

##--Write clip2.gml graph to file
igraph::write.graph(clip2, sprintf("%s/HPOclip2.gml",OUT[1]), format="gml")

inClip = match(V(gg)$term,V(clip)$term)
inClip = ifelse(is.na(inClip),0,1)
gg     = removeVertexTerm(gg,"clip")
gg     = set.vertex.attribute(gg,"clip",V(gg),inClip)

##---Write gg.gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",OUT[1],subDIR[S]), format="gml")


