##
# Calculate each diease-pair overlap/seperation on a selected  
# synaptic PPI network models, based on analysis described in: 
# Menche, J. et al. Uncovering disease-disease relationships through the incomplete interactome.
# Science, 347, (6224):1257601 (2015).
##

source('../setUp.R')

#Overlap of Disease A and B in the interactome
# GG   => igraph network
# GDA  => gda data for this graph
# disA => name of disease A
# disA => name of disease B
# OO   => minimum shorest paths for each gda, and each disease
diseaseOverlap <- function(GG, GDA, disA, disB, OO){

#disease A genes 
IDS1  <- V(GG)$name[grepl(disA,GDA,fixed=T)]
NIDS1 <- length(IDS1)

#disease B genes 
IDS2  <- V(GG)$name[grepl(disB,GDA,fixed=T)]
NIDS2 <- length(IDS2)

    #disease A given B
    paths  <- igraph::shortest.paths(GG,IDS1,IDS2,weights=NA)
    dsA    <- as.numeric(as.vector(apply(paths,1,min)))   
        
    #disease B given A
    paths <- igraph::shortest.paths(GG,IDS2,IDS1,weights=NA) 
    dsB   <- as.numeric(as.vector(apply(paths,1,min)))
    
    #network-based separation between disease A and B 
    dAB <- (sum(dsA)+sum(dsB))/(NIDS1+NIDS2)

    #network-based localisation of disease A
    indA <- which(colnames(OO)==disA)    
    dA   <- mean(as.numeric(as.vector(OO[OO[,indA[1]]!=".",indA[1]])))

    #network-based localisation of disease B
    indB <- which(colnames(OO)==disB)    
    dB   <- mean(as.numeric(as.vector(OO[OO[,indB[1]]!=".",indB[1]])))

    #overlap between disease A and B
    sAB = as.numeric(dAB) - (as.numeric(dA)+as.numeric(dB))/2

    return(sAB)
    
}

#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
resdir <- subDIR[S]
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}
#---


#---load corresponding graph which was used to build the consensus matrices from 
if( subDIR[S] == "HPO" ){
    st1 = sprintf("%s/%s.gml",OUT[1],subDIR[S])
} else {
    st1 = sprintf("%s/%s/%s.gml",OUT[1],clipStudy[V,2],subDIR[S])
}

cat("loading graph: ", st1, "...")
gg <- igraph::read.graph(st1,format="gml")
cat("done.\n")

#--- Find all HP Associations
HPA <- V(gg)$OV

#--- The number of HPA's in graph
NN  <- length(which(HPA!=""))


#---Check
#---Remove Diseases with zero HPA's
remove <- c()
cat("Following Models have zero HPA's.\n")
for( p in 1:length(ptype) ){
    IDS <- V(gg)$name[grepl(ptype[p],HPA,fixed=T)]
    if( length(IDS) == 0 ){
        cat(dtype[p], " => ", length(IDS),"\n")
	remove <- c(remove,p)
     }
}

if( length(remove) == 0 ){
    cat(" => None.\n")
}

if( length(remove) > 0 ){
    pid   <- pid[-remove]
    ptype <- ptype[-remove]
    phenl <- phenl[-remove]
}    
#---     

#--- outfile
res     <- matrix(0 ,ncol=4, nrow=length(ptype))
colnames(res) <- c("Model","N","mean_ds","SD_ds")
res[,1] <- ptype


#store minimum shorest paths for each HPA, and each disease
oo <- matrix(".",nrow=NN,ncol=(length(ptype)+2))
colnames(oo) <- c("ID","HPO.ID",ptype)
oo[,1] <- V(gg)$name[HPA !=""]
oo[,2] <- V(gg)$term[HPA !=""]

#loop over each model
for( p in 1:length(ptype) ){

    IDS <- V(gg)$name[grepl(ptype[p],HPA,fixed=T)]
    N   <- length(IDS)
    #ds  <- rep(0,N)
     
    #for each HPA, find the minimum shortest path to next HPA (of the same model)
    XX=igraph::shortest.paths(gg,IDS,IDS,weights=NA)
    diag(XX) = NA
    ds   = apply(XX,1,min,na.rm=T)
    indX = match(names(ds),oo[,1])
    oo[indX,(2+p)] = as.vector(ds)
    #for( i in 1:N ){
    #
    #    ds[i] <- min(as.vector(igraph::shortest.paths(gg,IDS[i],IDS[-i],weights=NA)))
    #
    #    indX <- which(oo[,1]==IDS[i])
    #    oo[indX,(2+p)] <- ds[i]
    #}

    res[p,2] <- as.character(N)
    res[p,3] <- as.character(mean(ds))
    res[p,4] <- as.character(sd(ds))

}

#outfile, the model localisation info, mean shortest ds for each model
outfile <- file(sprintf("%s/%s_model_localisation.csv",resdir,subDIR[S]),"w")
write.table(res, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

#outfile, the HPA shortest distance info
outfile <- file(sprintf("%s/%s_HP_model_separation.csv",resdir,subDIR[S]),"w")
write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

#model-model overlap
DAB <- matrix(".",ncol=length(ptype),nrow=length(ptype))
colnames(DAB) <- ptype
rownames(DAB) <- ptype

#--- NOTE ---#
# DAB is bound by -dmax <= DAB <= dmax
# where dmax denotes the diameter of the network
# dmax <- diameter(gg,directed=F)

# calculate disease-disease overlap
for( i in 1:length(ptype) ){
    for( j in i:length(ptype) ){

        DAB[i,j] <- 0

        #X = unlist(strsplit(rownames(DAB)[i],":"))
        #Y = unlist(strsplit(rownames(DAB)[j],":"))
        
        if( (i!=j) ){
            DAB[i,j] <- diseaseOverlap(gg,HPA,rownames(DAB)[i],colnames(DAB)[j],oo)
        }
        
    }
}

#outfile, the model overlap info
outfile <- file(sprintf("%s/%s_model_separation.csv",resdir,subDIR[S]),"w")
write.table(DAB, file=outfile, append=T, row.names=T, col.names=T, sep="\t",quote=F);
close(outfile);



