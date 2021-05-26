##setup
source("../setUp.R")

##REFS: [1] https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000641#s3
##      [2] https://dennyzhou.github.io/papers/LLGC.pdf
##      [3] https://www.pnas.org/content/115/27/E6375
##      [4] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6719725/#MOESM1
##      [5] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3806810/
##      [6] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3806810/
##      [7] http://dsd.cs.tufts.edu/
##      [8] https://github.com/reemagit/DSD/blob/master/DSD/calculator.py
##      [9] https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/

## Find Largest CC
findLCC <- function(gg){

    dec <- igraph::decompose.graph(gg)
    d=1
    CC=length(V(dec[[1]]))
    for( i in 1:length(dec) ){
        if(length(V(dec[[i]])) > CC){
            d=i
            CC=length(V(dec[[i]]))
        }
    }   

    gg  <- igraph::decompose.graph(gg)[[d]]
    return(gg)

}


## normalise the Adjacency matrix (Aij):
## Sij = Aij/sqrt(Dii * Djj)
## Example we="litcost"
normAij <- function(gg,we=NULL,label="term",makeSC=TRUE, back.edge.we=1e-2,
                    graph.simplify=FALSE, graph.LCC=FALSE,PRINT=FALSE,
                    type=c("r.norm","c.norm","s.norm") ){

    ##r.norm = row norm
    ##c.norm = column norm
    ##s.norm = symmetric norm
    
    type =  match.arg(type);
    
    if( graph.simplify == TRUE ){
        gg = igraph::simplify(gg,remove.multiple=T,remove.loops=T);        
        if( graph.LCC == TRUE ){
            gg = findLCC(gg);
        }
    }

    N  = length(V(gg))
    E  = length(E(gg))

    map = cbind(seq(1,N,1),V(gg)$name, get.vertex.attribute(gg,label,V(gg)) )
    ed  = get.edgelist(gg,names=T)

    ii  = map[match(ed[,1],map[,2]),1]
    jj  = map[match(ed[,2],map[,2]),1]

    if (is.null(we) ){
        edges = cbind(ii,jj,rep(1,E))
    } else {
        edges = cbind(ii,jj,get.edge.attribute(gg,we,E(gg)))
    }

    ## if network directed then make it strongly connected, i.e.
    ## edge directions were kept and low-weight back edges were added
    ## so that the network is strongly connected.
    if( is.directed(gg) == TRUE && makeSC == TRUE ){
        min.edge.we    = min(as.numeric(edges[,3]),na.rm=T);
        min.edge.we    = min.edge.we * back.edge.we;
        back.edges     = cbind(edges[,2],edges[,1],rep(min.edge.we))
        edges          = rbind(edges,back.edges)
    }

    ## build edge set
    edges = as.data.frame(edges)
    colnames(edges) = c("i","j","we")
    
    
    ## build sparse Adjacency matrix (set symmetric=TRUE for undirected networks) 
    if( is.directed(gg) == FALSE ){
        Aij = sparseMatrix(i=as.numeric(edges$i),j=as.numeric(edges$j),
                           x=as.numeric(edges$we),dims=c(N,N),
                           dimnames = list(map[,3],map[,3]),
                           symmetric=TRUE)
    } else {
        Aij = sparseMatrix(i=as.numeric(edges$i),j=as.numeric(edges$j),
                           x=as.numeric(edges$we),dims=c(N,N),
                           dimnames = list(map[,3],map[,3]),
                           symmetric=FALSE)
    }
        
    ## out-degree
    Dii = apply(Aij,1,sum)

    ## in-degree
    Djj = apply(Aij,2,sum)

    ## out-degree matrix
    #D1  = sparseMatrix(i=seq(1,N,1),j=seq(1,N,1),x=Dii^(-1/2), dims=c(N,N),
    #                   dimnames = list(map[,3], map[,3]),
    #                   symmetric=TRUE)

    ## in-degree matrix
    #D2  = sparseMatrix(i=seq(1,N,1),j=seq(1,N,1),x=Djj^(-1/2), dims=c(N,N),
    #                   dimnames = list(map[,3], map[,3]),
    #                   symmetric=TRUE)

    W  = NULL
    D1 = NULL
    D2 = NULL
    
    switch(type, 
           "r.norm"={
               ## row normalized adjacency matrix
               D1  = sparseMatrix(i=seq(1,N,1),j=seq(1,N,1),x=Dii^(-1), dims=c(N,N),
                       dimnames = list(map[,3], map[,3]),
                       symmetric=TRUE)
               W = D1 %*% Aij;
               if(PRINT){ cat("> r.norm.\n"); }
           },
           "c.norm"={
               ##  column normalized adjacency matrix
               D2  = sparseMatrix(i=seq(1,N,1),j=seq(1,N,1),x=Djj^(-1), dims=c(N,N),
                                  dimnames = list(map[,3], map[,3]),
                                  symmetric=TRUE)
               W = Aij %*% D2;
               if(PRINT){ cat("> c.norm.\n"); }
           },
            "s.norm"={
                ##  compute symmetric normalized adjacency matrix
                D1  = sparseMatrix(i=seq(1,N,1),j=seq(1,N,1),x=Dii^(-1/2), dims=c(N,N),
                       dimnames = list(map[,3], map[,3]),
                       symmetric=TRUE)
                
                D2  = sparseMatrix(i=seq(1,N,1),j=seq(1,N,1),x=Djj^(-1/2), dims=c(N,N),
                                  dimnames = list(map[,3], map[,3]),
                                  symmetric=TRUE)
                W = D1 %*% Aij %*% D2;
                if(PRINT){ cat("> s.norm.\n"); }
            }
           )
          

    ## get row,column pairs from sparse matrix
    df <- as.data.frame(summary(W))
    df$hp.term <- rownames(Aij)[df$i]
    df$hp.term <- colnames(Aij)[df$j]

    return(list(gg=gg,Aij=Aij,W=W,edges=df,type=type,Dout=Dii,Din=Djj))

}

propagate <- function(gg=NULL, W=NULL, SET=NULL, alpha=c(0.9), thres=1e-5,
                      max.its=400, PRINT=FALSE){

    ## REF [1] https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000641#s3
    
    S = NULL;
    
    if( !is.null(gg) && !is.null(W) && !is.null(SET) ){

        N = length(V(gg))
        E = length(E(gg))       
        
        Y = rep(0,N)
        Y[match(SET,V(gg)$term)] = 1

        S = rep(0,N)

        t=1
        diff=1
        while( t < max.its && diff > thres ){

            if( (t-1) == 0 ){
                Sold = Y
                S    = alpha * W %*% Y + (1-alpha) * Y;
            } else {
                Sold = S
                S    = alpha * W %*% S + (1-alpha) * Y; 
            }

            diff = norm2(S-Sold);
            if( PRINT ){ cat("t:", t, "diff: ", diff, "\n"); }
            t=t+1

        }
       
    }#if

    return(list(S=S,step=t,diff=diff))
    
}

propagate.classes <- function(gg=NULL, W=NULL, D1=NULL, D2=NULL,
                              SETS=NULL, alpha=c(0.9), thres=1e-5,
                              max.its=100, PRINT=FALSE, calCOST=FALSE ){

    ## REF [1] https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000641#s3
    
    S = NULL;
    
    if( !is.null(gg) && !is.null(W) && !is.null(SETS) &&
       !is.null(D1) && !is.null(D2) ){

        N = length(V(gg))
        E = length(E(gg))       

        ## Y is now a nxc matrix, where n= number of nodes, c=number of sets
        C  = length(SETS)
        Y  = matrix(0,nrow=N,ncol=C)
        S  = as(Y,"dgCMatrix")
        for( i in 1:C ){
            Y[match(SETS[[i]],V(gg)$term),i] = 1
        }

        ## out-degree
        Yi = apply(Y,1,sum)
        
        ## convert Y to sparse matrix
        Y = as(Y,"dgCMatrix")

        t=1
        diff=1
        while( t < max.its && diff > thres ){

            if( (t-1) == 0 ){
                Sold = Y
                S    = alpha * W %*% Y + (1-alpha) * Y;
            } else {
                Sold = S
                S    = alpha * W %*% S + (1-alpha) * Y; 
            }

            diff = norm2(S-Sold);
            cost = 0;
            if( calCOST ){ cost = S.cost(W=W,F=S,Yi=Yi,D1=D1,D2=D2); }
            if( PRINT ){ cat("t:", t, "diff: ", diff, "cost: ", cost, "\n"); }
            t=t+1

        }
       
    }#if

    return(list(S=S,Yi=Yi,step=t,diff=diff))
    
}


norm2 <- function(x) sqrt(sum(x^2))

smoothness <- function(x,ni,nj){
    return(x * norm2(ni-nj));
}

S.cost <- function(W=NULL,F=NULL,Yi=NULL,D1=NULL,D2=NULL,mu=0.5){

    cost=NA

    if( !is.null(W) && !is.null(F) && !is.null(Yi) &&
         !is.null(D1) && !is.null(D2) ){

        N = dim(F)[1] #number of nodes
        C = dim(F)[2] #number of classes
        
        ## out-degree
        Fi    = as.numeric(apply(F,1,sum))
        normi = Fi/diag(D1)
        
        ## in-degree
        ##Fj    = as.numeric(apply(F,2,sum))
        normj = Fi/diag(D2)

        smooth = 0
        tmp    = matrix(0,N,N)
        jj     = seq_along(tmp[1,])
        for( i in 1:N ){
            val = as.vector(sapply(normi[i]-normj[jj],norm2));
            tmp[i,jj] = val
            tmp[jj,i] = val
            jj        = jj[-1]            
        }
        
        tmp    = as(tmp,"dgCMatrix")
        smooth = smooth + sum( (W * tmp), na.rm=T)
        
        cost = 0.5 * smooth + mu * sum(sapply(Fi-Yi,norm2),na.rm=T)
    }

    return(cost);
    
}

## nieve radidal basis function (rdf) kernel
nieve.kernel <- function(x,sigma=1){ return( exp(-(x/sigma)) ); }

## calculate diffusion state distance (DSD) and similarity
calDSD <- function(gg,we=NULL,label="term",nRW=5,sigma=1,
                   type=c("r.norm","c.norm","s.norm") ){

    ## REF [1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6719725/#MOESM1

    type =  match.arg(type);
    
    Dist = NULL;
    Sim  = NULL;
    
    Anorm = normAij(gg,we=we,label=label,makeSC=TRUE,type=type);

    #gg = Anorm$gg
    #N  = length(V(gg))
    #E  = length(E(gg))

    
    ## weighted normalised Adjacency matrix
    P  = Anorm$W;
    N  = dim(P)[1];

    ## node's out-degree 
    Deg = Anorm$Dout
    
    ## node's out-degree 
    #Deg = apply(W,1,sum);

    #p.norm <- function(x,norm){ return(x/norm); }
    #P  = matrix(0,N,N)
    #for( i in 1:N ){
    #    P[i,] = as.vector(W[i,])/Deg[i];
    #    ##P[i,] = as.vector(sapply(W[i,],p.norm,norm=Deg[i]));
    #}

    #P[lower.tri(P)] = t(P)[upper.tri(P)]

    ## convert P to sparse matrix
    P = as(P,"dgCMatrix")

    if( nRW > 0 ){

        C = diag(N);
        C = as(C,"dgCMatrix");
        for( i in 1:nRW ){
          C = (C %*% P) + as(diag(N),"dgCMatrix");
        }
        ## DSD 
        ##Dist = as.matrix(stats::dist(as.matrix(C),method="manhattan")); ##too slow
        Dist = as.matrix(wordspace::dist.matrix(C, method="manhattan", as.dist=TRUE)); # Manhattan distance, works
        ##Dist = as.matrix(proxyC::dist(C, method="manhattan")) ## also works
        ## DSD to similarity, using radial basis function kernal    
        ##Sim = rdetools::rbfkernel(Dist) #crash
        ##rbf = rbfdot(sigma=sigma)
        ##Sim = kernelMatrix(x=Dist,kernel=rbf) #tiny values??
        Sim = matrix(0,N,N)
        for( i in 1:N ){ Sim[,i] = nieve.kernel(x=as.vector(Dist[,i]),sigma=sigma); }
        rownames(Sim) = rownames(Dist)
        colnames(Sim) = colnames(Dist)  
    } else {
        Pi  = as(diag(Deg),"dgCMatrix")
        Pi  = Pi/sum(Deg,na.rm=T);
        Inv = Matrix::triu(as(diag(N),"dgCMatrix") - P - Pi);##chol2inv requires a triangular matrix
        Inv = Matrix::chol2inv(Inv);        
        ## DSD 
        ##Dist = as.matrix(stats::dist(as.matrix(Inv),method="manhattan"));
        Dist = as.matrix(wordspace::dist.matrix(Inv, method="manhattan", as.dist=TRUE)); # Manhattan distance
        
        ## DSD to similarity, using radial basis function kernal    
        ##Sim  = rdetools::rbfkernel(Dist)
        ##rbf = rbfdot(sigma=sigma)
        ##Sim = kernelMatrix(x=Dist,kernel=rbf)
        Sim = matrix(0,N,N)
        for( i in 1:N ){ Sim[,i] = nieve.kernel(x=as.vector(Dist[,i]),sigma=sigma); }
        rownames(Sim) = rownames(Dist)
        colnames(Sim) = colnames(Dist)  
    }
    
    return(list(Dist=Dist,Sim=Sim));

}

## directories
annodir <- DIRS[grepl("Annotation",DIRS)]
grdir   <- DIRS[grepl("Graphs",DIRS)]

## read in HPO network
gg     = read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
N      = length(V(gg))
E      = length(E(gg))
terms  = V(gg)$term
levels = V(gg)$level

##---Read-in orginal annotation file
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
colnames(Anno)[which(colnames(Anno)=="ddd_freq")]     = models[1]
colnames(Anno)[which(colnames(Anno)=="metamap_freq")] = models[2]
Anno = as.data.frame(Anno)

## read-in model HP term per gene annotation data
LITset = readRDS(sprintf("%s/LIT.Gene.HPterms.RDS",annodir))
DDDset = readRDS(sprintf("%s/DDD.Gene.HPterms.RDS",annodir))

genes  = names(LITset)
Ngenes = length(genes)


Anorm = normAij(gg=gg,type="s.norm")
#res   = propagate.classes(gg=gg,W=Anorm$W,D1=Anorm$D1,D2=Anorm$D2,SETS=LITset,
#                        PRINT=TRUE)
