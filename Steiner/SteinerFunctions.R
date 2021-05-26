removeVertexTerm <- function(GG,NAME){

    if( !is.null(igraph::get.vertex.attribute(GG,NAME)) ){
        GG <- igraph::remove.vertex.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.vertex.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.vertex.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}

findUnique <- function(X,decreasing=TRUE){

    if( grepl(COLLAPSE[c], X) ) {
        Xsplit = strsplit(X,COLLAPSE[c])[[1]]
        X      = unique(c(Xsplit))
        if( decreasing ){
            X      = X[order(X,decreasing=TRUE)]
        } else {
            X      = X[order(X,decreasing=FALSE)]
        }
        return(X)
    } else {                
        return(X)
    }    
            
}

findandCollapse <- function(X,decreasing=TRUE){

    if( grepl(COLLAPSE[c], X) ) {
        Xsplit = strsplit(X,COLLAPSE[c])[[1]]
        X      = unique(c(Xsplit))
        if( decreasing ){
            X      = X[order(X,decreasing=TRUE)]
        } else {
            X      = X[order(X,decreasing=FALSE)]
        }

        return(collapseAnnotation(X))
    } else {
        return(X)
    }
                 
}

collapseAnnotation <- function(X){
    if( length(X) == 1 ){ return(X) }    
    else {
        X = X[X!=""]
        if( length(X) == 1 ){ return(X) }
        else { return(paste(X,collapse=COLLAPSE[c])) }    
    }
}

mergeAnnotation <- function(X,Y,decreasing=TRUE){

    X = as.character(findUnique(X,decreasing=decreasing))
    Y = as.character(findUnique(Y,decreasing=decreasing))

    temp = unique(c(X,Y))
    if( decreasing ){
        temp = temp[order(temp,decreasing=TRUE)]
    } else {
        temp = temp[order(temp,decreasing=FALSE)]
    }

    return(collapseAnnotation(temp))
    
    #if( length(temp) == 1 ){ return(temp) }    
    #else {
    #    temp = temp[temp!=""]
    #    if( length(temp) == 1 ){ return(temp) }
    #    else { return(paste(temp,collapse=COLLAPSE[c])) }    
    #}

}


#---Find Largest CC
findLCC <- function(GG){

    dec <- igraph::decompose.graph(GG)
    d=1
    CC=length(V(dec[[1]]))
    for( i in 1:length(dec) ){
        if(length(V(dec[[i]])) > CC){
            d=i
            CC=length(V(dec[[i]]))
        }
    }   

    GG  <- igraph::decompose.graph(GG)[[d]]
    return(GG)

}

findLCCinfo <- function(GG){

    dec <- igraph::decompose.graph(GG)

    N   <- length(dec)  
    oo  <- matrix(NA,nrow=N,ncol=3)
    
    for( i in 1:N ){
        oo[i,1] = i
        oo[i,2] = length(V(dec[[i]]))
        oo[i,3] = length(E(dec[[i]]))
    }

    oo = oo[order(oo[,2],decreasing=T),]
    
    return(oo)
}

getTerminals <- function(GG,TARGETS){

    TARGETS <- as.vector(TARGETS)

    STERM = V(GG)$name[match(TARGETS,V(gg)$term)]
    STERM = STERM[!is.na(STERM)]

    return(STERM)
}

buildSteinerTrees <- function(GG,STERM,NTREES){

    TREES  <- list()
    for( i in 1:NTREES ){
        TREES[[i]] = steinertree(type="SP", terminals=STERM,graph=GG,color=F,merge=F)[[1]]
    }

    return(TREES)
    
}

connectSmallFagments <- function(GG,CLIP,TARGETS){

    ed = get.edgelist(GG,names=F)
    ed = cbind(V(GG)$term[ed[,1]],V(GG)$term[ed[,2]])

    CLIPterms <- V(CLIP)$term
    indx      <- match(TARGETS,CLIPterms)
    indx      <- ifelse(is.na(indx),TRUE,FALSE)
    Terms     <- TARGETS[indx]

    edExt <- c()
    for( i in 1:length(Terms) ){
        tmp = ed[which(ed[,1]==Terms[i] | ed[,2]==Terms[i]),]
        if( i == 1 ){
            edExt = tmp
        } else {
            edExt = rbind(edExt,tmp)
        }
    }
    
    edExt = cbind(as.vector(edExt[,1]),as.vector(edExt[,2]))
    edExt = edExt[!duplicated(edExt),]

    if( length(edExt[,1]) > 0 ){
    
        ids  = unique(c(as.vector(edExt[,1]),as.vector(edExt[,2])))
        Xids = cbind(seq(1,length(ids),1),ids)
        ed2  = cbind(Xids[match(edExt[,1],Xids[,2]),1],Xids[match(edExt[,2],Xids[,2]),1])

        gg2  = graph_from_edgelist(ed2,directed=F)
        gg2  = set.vertex.attribute(gg2,"term",V(gg2), Xids[match(V(gg2)$name,Xids[,1]),2] )    

        ggL      = list()
        ggL[[1]] = gg2
        
        xx = connectFragments(GG=GG,CLIPF=CLIP,TREES=ggL,INDX=1)

    } else {
        xx = CLIP
    }
        
    return(xx)
    
}

graph.frac.annoAtLevels <- function(GG){

    tot  = table(V(GG)$OV)
    N    = length(tot)
    
    maxL = max(V(GG)$level)

    oo   = matrix(0,ncol=(N+1), nrow=(maxL+1))
    colnames(oo) = c("level",names(tot))
    
    for( i in 1:(maxL+1) ){
        oo[i,1]       = (i-1)

        Levi = table(V(GG)$OV[V(GG)$level==(i-1)])

        indx = match(names(tot),names(Levi))

        for( j in 1:length(indx) ){
            if( !is.na(indx[j]) ){
                if( tot[j] != 0 ){
                    oo[i,(j+1)] = Levi[indx[j]]/tot[j]
                }
            }
        }
        
    }

    names(tot)[names(tot) == ""]="non"
    colnames(oo) = c("level",names(tot))
    
    return(oo)
    
}

frac.annoAtLevels <- function(MAP){

    MAP  = as.data.frame(MAP)    
     tot = table(MAP$OV)
    N    = length(tot)    
    maxL = max(as.numeric(MAP$level))

    oo   = matrix(0,ncol=(N+1), nrow=(maxL+1))
    colnames(oo) = c("level",names(tot))
    
    for( i in 1:(maxL+1) ){
        oo[i,1]       = (i-1)

        Levi = table(MAP$OV[as.numeric(MAP$level)==(i-1)])
        indx = match(names(tot),names(Levi))

        for( j in 1:length(indx) ){
            if( !is.na(indx[j]) ){
                if( tot[j] != 0 ){
                    oo[i,(j+1)] = Levi[indx[j]]/tot[j]
                }
            }
        }
        
    }

    names(tot)[names(tot) == ""]="non"
    colnames(oo) = c("level",names(tot))
    
    return(oo)    
    
}

map.annoAtLevels <- function(MAP,level.max){

    MAP          = as.data.frame(MAP)    
    LEVELS       = (level.max+1)

    tt1CN        = levels(factor(MAP$OV))
    tt1CN        = tt1CN[tt1CN!=""]
    N            = (length(tt1CN)+1)
    
    tt1          = matrix(0,ncol=N,nrow=LEVELS)
    colnames(tt1)= c("level",tt1CN)
    tt1[,1]      = seq(0,level.max,1)

    for(l in 1:LEVELS ){
        atLevel = MAP$OV[as.numeric(MAP$level)==(l-1)]
        for( i in 2:N ){
            tt1BIN   = grepl(colnames(tt1)[i],atLevel)
            tt1[l,i] = sum(tt1BIN)
        }
    }
        
    return(list("MODEL"=tt1))
    
}


map.annoAtLevels.perGene <- function(MAP,level.max){

    MAP          = as.data.frame(MAP)    
    LEVELS       = (level.max+1)
    N            = (Ngenes+1)

    tt1CN        = ptype[1:Ngenes]
    tt1          = matrix(0,ncol=N,nrow=LEVELS)
    colnames(tt1)= c("level",tt1CN)
    tt1[,1]      = seq(0,level.max,1)

    tt2CN        = ptype[(Ngenes+1):length(ptype)]
    tt2          = matrix(0,ncol=N,nrow=LEVELS)
    colnames(tt2)= c("level",tt2CN)
    tt2[,1]      = seq(0,level.max,1)

    tt3CN        = rep("",Ngenes)
    for( i in 1:Ngenes ){ tt3CN[i] = mergeAnnotation(tt1CN[i],tt2CN[i]) }
    tt3          = matrix(0,ncol=N,nrow=LEVELS)
    colnames(tt3)= c("level",tt3CN)
    tt3[,1]      = seq(0,level.max,1)
    
    for(l in 1:LEVELS ){
        atLevel = MAP$OVG[as.numeric(MAP$level)==(l-1)]
        for( i in 2:N ){
            tt1BIN   = grepl(colnames(tt1)[i],atLevel)
            tt1[l,i] = sum(tt1BIN)

            tt2BIN   = grepl(colnames(tt2)[i],atLevel)
            tt2[l,i] = sum(tt2BIN)

            tt3BIN   = tt1BIN & tt2BIN
            tt3[l,i] = sum(tt3BIN)
            
        }
    }
        
    return(list("LIT"=tt1, "DDD"=tt2, "LIT_DDD"=tt3))
    
}

norm.annoAtLevels <- function(MAP){

    DIM  = dim(MAP)
    NC   = DIM[2]
    sums = colSums(MAP[,2:NC])
    for( i in 2:NC ){
        if( sums[(i-1)] != 0 ){
            MAP[,i] = MAP[,i]/sums[(i-1)]
        }
    }

    return(MAP)
    
}

plot.annoAtLevels <- function(DF,XLAB=NULL,YLAB="",LEGENDtit="count"){

    DF           <- DF[,-1]
    if( is.null(XLAB) ){ XLAB <- seq(1,length(colnames(DF)),1) }
    colnames(DF) <- XLAB    
    DF           <- melt(DF)
    DF           <- as.data.frame(DF)

    Xbreaks      <- unique(DF$Var1)
    Xlabs        <- Xbreaks-1
    
    #LegendTit    <- "count"    
    #YLAB="LITvDDD per gene"
    
    gplot <- ggplot(DF, aes(x=factor(Var1), y=Var2)) +        
        geom_tile(aes(fill = value), colour = "white") + 
        scale_fill_distiller(palette = "YlGnBu", direction = 1, name=LEGENDtit) +
        labs(y=YLAB,x="HPO level")+
        theme(legend.key=element_blank())+
        scale_x_discrete(breaks = Xbreaks, labels = Xlabs) + 
        theme_minimal() +
        coord_equal()+
         theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=rel(0.7)),
              axis.text.y = element_text(face="bold",size=rel(1)),
              axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.title=element_text(face="bold",size=rel(1.0)),
              legend.text=element_text(face="bold",size=rel(1.0)),
              legend.position="right")

    return(gplot)
    
}

plot.save <- function(GPLOT=NULL,PLOTDIR="PLOTS",FILENAME="temp.png",
                      pHEIGHT=HEIGHT,pWIDTH=WIDTH,UNITS="px"){

    if( !is.null(GPLOT) ){
        
        str = sprintf("%s/%s.png",PLOTDIR,FILENAME)        
    
        #png(str,width=WIDTH,height=2*HEIGHT,units="px");
        png(str,width=pWIDTH,height=pHEIGHT,units=UNITS);
        print(GPLOT)
        dev.off()
    }
    
}

roleUpTerms <- function(GG,level.max=NULL,level.target=NULL){

    if( is.null(level.max)    ){ level.max = max(V(GG)$level) }

    if( is.null(level.target) ){ level.target = 10 }
    
    map <- cbind(seq(1,length(V(GG)),1),V(GG)$name,
                 V(GG)$term,V(GG)$level,
                 V(GG)$OV, V(GG)$OVG,
                 rep(0,length(V(GG))))

    colnames(map) = c("id","name","term","level","OV","OVG","visited")

    cat("level.max = ", level.max, ", level.target = ", level.target,"\n")

    if( level.target < level.max ){
    
    for( l in level.max:level.target ){

        atLevel = map[as.numeric(map[,4])==l,]            
        N       = dim(atLevel)[1]

        for( i in 1:N ){
        
            deg  = igraph::neighbors(GG,atLevel[i,2],mode="all")
            neig = map[as.vector(deg),]

            if( is.vector(neig) ){

                if( as.numeric(neig[4]) < l ){

                    X = map[as.numeric(atLevel[i,1]),5]
                    Y = map[as.numeric(neig[1]),5]
                    map[as.numeric(neig[1]),5] = mergeAnnotation(X,Y)                   
                        
                    X = map[as.numeric(atLevel[i,1]),6]
                    Y = map[as.numeric(neig[1]),6]
                    map[as.numeric(neig[1]),6] = mergeAnnotation(X,Y)

                }

            } else {

                if( length(neig[,1]) > 0 ){
                
                    for( n in 1:length(neig[,1])){
                        if( as.numeric(neig[n,4]) < l ){

                            X = map[as.numeric(atLevel[i,1]),5]
                            Y = map[as.numeric(neig[n,1]),5]
                            map[as.numeric(neig[n,1]),5] = mergeAnnotation(X,Y)

                            X = map[as.numeric(atLevel[i,1]),6]
                            Y = map[as.numeric(neig[n,1]),6]
                            map[as.numeric(neig[n,1]),6] = mergeAnnotation(X,Y)

                        }#if
                    }#anno
                }#if
            }#else        
            map[as.numeric(atLevel[i,1]),5] = ""
            map[as.numeric(atLevel[i,1]),6] = ""
            map[as.numeric(atLevel[i,1]),7] = 1
        }#i
    }#l

    }#if
        
    map[,5] = sapply(map[,5],findandCollapse)
    map[,6] = sapply(map[,6],findandCollapse)
    
    return(map)
}

#old code
#connectSmallFagments2 <- function(GG,CLIP,TARGETS,RUNS=1000,N=1){
#
#    RUN       <- 1
#    CLIPNEW   <- CLIP 
#    
#    CLIPterms <- V(CLIPNEW)$term
#    indx      <- match(TARGETS,CLIPterms)
#    indx      <- ifelse(is.na(indx),TRUE,FALSE)
#    Terms     <- TARGETS[indx]
#
#    ROOTmin   <- min(V(CLIP)$level)
#    ROOT      <- V(CLIP)$term[which(V(CLIP)$level==ROOTmin)]
#
#    Terms     <- c(ROOT,Terms)
#    
#    cat("run: ", RUN, ", No Terms: ", length(Terms), "\n")
#    
#    while( length(Terms) > 0 && RUN < RUNS ){         
#
#        if( length(Terms) > 5 ){
#            Terms <- sample(Terms,5)
#            Terms <- c(ROOT,Terms)
#        }
#        
#        sTerm     <- getTerminals(GG,Terms)
#        if( length(sTerm) > 0 ){
#            trees     <- buildSteinerTrees(GG,sTerm,N)
#            tree_len  <- unlist(lapply(trees, function(x) length(E(x[[1]]))))
#            indxMin   <- which(tree_len == min(tree_len))#
#
#if( length(indxMin) == 1 ){
#        
#                CLIPNEW   <- connectFragments(GG=GG,CLIPF=CLIPNEW,TREES=trees,INDX=indxMin)
#
#                CLIPterms <- V(CLIPNEW)$term
#                indx      <- match(TARGETS,CLIPterms)
#                indx      <- ifelse(is.na(indx),TRUE,FALSE)
#                Terms     <- TARGETS[indx]
#                
#                RUN = RUN + 1
#
#                cat("run: ", RUN, ", No Terms: ", length(Terms), "\n")
#            }
#        }
#    }
#        
#        
#    return(CLIPNEW)
#                         
#}


connectFragments <- function(GG,CLIPF,TREES,INDX){

    ##map edges from res to HP terms
    res  = TREES[[INDX]]#[[1]]
    Sed  = igraph::get.edgelist(res,names=F)
    SedM = cbind(V(res)$term[Sed[,1]],V(res)$term[Sed[,2]])

    ##map edges from clipF to HP terms
    Ced  = igraph::get.edgelist(CLIPF,names=F)
    CedM = cbind(V(CLIPF)$term[Ced[,1]],V(CLIPF)$term[Ced[,2]])

    ##merge HP edges from res and clipF 
    Xed  = rbind(SedM,CedM)
    ids  = unique(c(Xed[,1],Xed[,2]))
    Xids = cbind(seq(1,length(ids),1),ids)
    XedM = cbind(Xids[match(Xed[,1],Xids[,2]),1],Xids[match(Xed[,2],Xids[,2]),1])

    xx   = igraph::graph_from_edgelist(XedM,directed=F)
    xx   = igraph::simplify(xx,remove.multiple=T,remove.loops=T)
    xx   = set.vertex.attribute(xx,"term",V(xx), Xids[match(V(xx)$name,Xids[,1]),2] )
    Xind = match(V(xx)$term,V(GG)$term)
    xx   = set.vertex.attribute(xx,"label",          V(xx),V(GG)$label[Xind])
    xx   = set.vertex.attribute(xx,"gene",           V(xx),V(GG)$gene[Xind])
    xx   = set.vertex.attribute(xx,"diseasename",    V(xx),V(GG)$diseasename[Xind])
    xx   = set.vertex.attribute(xx,"literaturemodel",V(xx),V(GG)$literaturemodel[Xind])
    xx   = set.vertex.attribute(xx,"dddmodel",       V(xx),V(GG)$dddmodel[Xind])
    xx   = set.vertex.attribute(xx,"OV",             V(xx),V(GG)$OV[Xind])
    xx   = set.vertex.attribute(xx,"OVG",            V(xx),V(GG)$OVG[Xind])
    xx   = set.vertex.attribute(xx,"level",          V(xx),V(GG)$level[Xind])
  
    clip = findLCC(xx)

    return(clip)
    
}

clip.from.subgraph <- function(GG,CLIP){

    terms = V(CLIP)$term 
    ids   = cbind(V(GG)$name,V(GG)$term)

    indx  = match(ids[,2],terms)
    indx  = ifelse(!is.na(indx),TRUE,FALSE)
    
    xx    = igraph::induced.subgraph(gg,ids[indx,1])

    clip = findLCC(xx)

    return(xx)
    
}
