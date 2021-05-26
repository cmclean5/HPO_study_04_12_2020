##
# Calculate each diease-pair overlap/seperation on a selected  
# synaptic PPI network models, based on analysis described in: 
# Menche, J. et al. Uncovering disease-disease relationships through the incomplete interactome.
# Science, 347, (6224):1257601 (2015).
##

source('../setUp.R')
library(WGCNA) #for it's qvalue calculation

strSplit <- function(X,AT=":",POS=1){

    X = strsplit(X,AT)    
    X = lapply(X,function(x,p=POS) x[p])

    return(unlist(X))
    
}

filter <- function(DF){

    X = rownames(DF)
    Y = colnames(DF)

    N = length(X)
    M = length(Y)
    
    tempX=strsplit(X,":")
    tempY=tempX
    
    X1 = unlist(lapply(tempX,function(x) x[1]))
    Y1 = unlist(lapply(tempX,function(x) x[2]))
    
    X2 = unlist(lapply(tempY,function(x) x[1]))
    Y2 = unlist(lapply(tempY,function(x) x[2]))    

    x=c()
    y=c()
    
    for( i in 1:N ){
        for( j in 1:M ){

            if( (i!=j) && (X1[i]!=X2[j]) && (Y1[i]==Y2[j]) ){
                if( length(which(x==i)) == 0 ){ x=c(x,i) }
                if( length(which(y==j)) == 0 ){ y=c(y,j) }      
            }

        }
    }

    retrun(cbind(x,y))

}

qscore <- function(zz,FDR){

    LL <- FDR[FDR[,1] < as.numeric(zz),2]

    if( length(LL) != 0 ){    
        return(LL[end(LL)[1]])
    }

    return(1)
}

actualOverlap <- function(HPA, SELECTION=1){

    oo <- NULL
    
    if( HPOperGene ){
    
        #oo <- matrix(NA,ncol=2,nrow=Ngenes)
        oo <- data.frame(a=character(),b=as.character(),c=as.character(),d=as.character())

        k=1
        for( i in 1:length(ptype) ){
            for( j in 1:length(ptype) ){

                if( i <= j ){

                    ov = sum(grepl(ptype[i],HPA) & grepl(ptype[j],HPA))
                    
                    X = unlist(strsplit(ptype[i],":"))
                    Y = unlist(strsplit(ptype[j],":"))
                    
                    #Str = c(X,Y)
                    #Str = paste(as.character(Str),collapse=COLLAPSE[c])
                    
                    if( SELECTION == 1 ){

                        if( (i!=j) && (X[1]!=Y[1]) && (X[2]==Y[2]) ){
                            tmp = data.frame(a=as.character(k),
                                             b=ptype[i],c=ptype[j],
                                             d=as.character(ov))
                            oo  = rbind(oo,tmp)
                            k=k+1
                        }
                    }#1
                        
                    if( SELECTION == 2 ){

                        if( (i!=j) && (X[1]!=Y[1]) ){
                            tmp = data.frame(a=as.character(k),
                                             b=ptype[i],c=ptype[j],
                                             d=as.character(ov))
                            oo  = rbind(oo,tmp)
                            k=k+1
                        } 
                    }#2

                }
            }
        }
                        
        #for( i in 1:Ngenes ){
        #
        #    litX = ptype[i]
        #    dddX = ptype[(i+Ngenes)]
        #
        #    Str  = c(litX,dddX)
        #
        #    oo[i,1] = paste(as.character(Str),collapse=COLLAPSE[c])
        #    oo[i,2] = sum(grepl(litX,HPA) & grepl(dddX,HPA))
        #
        #}
        #
    }

    return(oo)
    
}

#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]

#---load corresponding graph which was used to build the consensus matrices from 
if( subDIR[S] == "HPO" ){
    st1 = sprintf("%s/%s.gml",OUT[1],subDIR[S])
} else {
    st1 = sprintf("%s/%s/%s.gml",OUT[1],clipStudy[V,2],subDIR[S])
}

cat("loading graph: ", st1, "...")
gg <- igraph::read.graph(st1,format="gml")
cat("done.\n")

#--- Find all HP Associations (HPA)
HPA <- V(gg)$OVG


##---Read-in orginal annotation file
ANNO <- read.delim("../Annotations/hpo_lists_per_disease.csv",sep=",",header=T)
GN   <- ANNO[match(ANNO[,1],pid[1:length(ANNO[,2])]),2]
GN2  <- c(GN,GN)
DIS  <- ANNO[match(ANNO[,1],pid[1:length(ANNO[,2])]),3]
DIS2 <- c(DIS,DIS)

##---Check or create output dir called RESULTS
if( !file_test("-d","RESULTS") ){
    dir.create("RESULTS")
}

resdir <- sprintf("RESULTS/%s",subDIR[S])

if( !file_test("-d",resdir) ){
    dir.create(resdir)
}
#---

FILES <- vector(length=5)
FILES[1] <- "HP_model_separation_per_gene"  
FILES[2] <- "model_separation_per_gene"
FILES[3] <- "random_separation_per_gene"
FILES[4] <- "model_localisation_per_gene"
FILES[5] <- "actualOverlap_per_gene"

#For Model Sizes on the Ontology
loc <- read.table(sprintf("%s/%s_%s.csv",subDIR[S],subDIR[S],FILES[4]),sep="\t",header=T)

ds <- read.table(sprintf("%s/%s_%s.csv",subDIR[S],subDIR[S],FILES[1]),sep="\t",header=T)

RNDsub    = vector(length=2)
RNDsub[1] = "GENES"
RNDsub[2] = "LITvALLDDD"
rsub=2

if( subDIR[S] == "HPO" ){
    RANds <- read.table(sprintf("%s/%s/%s/%s/random_%s.csv",rndDIR[1],subDIR[S],RNDsub[rsub],subDIR[S],FILES[1]),sep="\t",header=T)
} else {
    RANds <- read.table(sprintf("%s/%s/%s/%s/random_%s.csv",rndDIR[1],subDIR[S],RNDsub[rsub],clipStudy[V,2],FILES[1]),sep="\t",header=T)
}
    
#Ovlap <- read.table(sprintf("%s/%s_%s.csv",subDIR[S],subDIR[S],FILES[5]),sep="\t",header=T)
Ovlap <- actualOverlap(HPA=HPA,SELECTION=rsub)

##comment out for moment
CN <- colnames(ds)[3:length(ds[1,])]

if( grepl(".",CN[1]) ){
CN <- gsub("\\.",":",CN)
}

res     <- matrix(0 ,ncol=10, nrow=length(CN))
colnames(res) <- c("ID","GENE.NAME","Model.long","Model","N","mean_ds","SD_ds","Ran_mean_ds","Ran_SD_ds","Utest.pvalue")

res[,1] <- pid[match(CN,ptype)]
res[,2] <- GN2
res[,3] <- strSplit(phenl[match(CN,ptype)])
res[,4] <- strSplit(CN)
res[,5] <- loc[match(CN,loc[,1]),2]


##significance of ds for each model
for( i in 1:length(CN) ){
     
     #gda matching indices
     indx <- ds[,(2+i)]!="."
     
     #gene ids
     ids <- ds[indx,1]

     #observed ds values
     DS <- as.numeric(as.vector(ds[indx,(2+i)]))
     res[i,6] <- as.numeric(mean(DS))
     res[i,7] <- as.numeric(sd(DS))

     indy <- match(ids,RANds[,1])
     
     #random ds values
     RDS <- as.numeric(as.vector(RANds[indy,(2+i)]))
     res[i,8]  <- as.numeric(mean(RDS))
     res[i,9]  <- as.numeric(sd(RDS))
     res[i,10] <- 1.0
    
    #compute wilcox test between observable ds and random ds, and store p.values, 
    #see (Menche et al., 2015).
    if( !is.infinite(DS) && !is.nan(DS) && !is.na(DS) &&  !is.infinite(RDS) && !is.nan(RDS) && !is.na(RDS) ){
        if( length(DS) != 0 && length(RDS) != 0 ){ 
            wt        <- wilcox.test(DS,RDS)
            res[i,10] <- as.numeric(wt$p.value)
        }
    }
}

outfile <- file(sprintf("%s/model_location_per_gene_sig.csv",resdir),"w");
write.table( res, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);


#--- load observed mean sAB, i.e. overlap between model A and B
sAB <- read.table(sprintf("%s/%s_%s.csv",subDIR[S],subDIR[S],FILES[2]),sep="\t",header=T)

if( grepl(".",rownames(sAB)) ){rownames(sAB) <- gsub("\\.",":",rownames(sAB))}
if( grepl(".",colnames(sAB)) ){colnames(sAB) <- gsub("\\.",":",colnames(sAB))}


#--- load random sAB values for each iteration
if( subDIR[S] == "HPO" ){
    RAW_sAB  <- read.table(sprintf("%s/%s/%s/%s/RAW_random_sAB_%s.csv",rndDIR[1],subDIR[S],RNDsub[rsub],subDIR[S],FILES[3]),sep="\t",header=T)
} else {
   RAW_sAB  <- read.table(sprintf("%s/%s/%s/%s/RAW_random_sAB_%s.csv",rndDIR[1],subDIR[S],RNDsub[rsub],clipStudy[V,2],FILES[3]),sep="\t",header=T)
}
    
#---no. of random permutations
perms <- length(RAW_sAB[1,])

RAN_sAB_mean <- matrix(".",ncol=3,nrow=length(RAW_sAB[,1]))
RAN_sAB_sd   <- matrix(".",ncol=3,nrow=length(RAW_sAB[,1]))

RAN_sAB_mean[,1] <- as.character(RAW_sAB[,1])
RAN_sAB_mean[,2] <- as.character(RAW_sAB[,2])

RAN_sAB_sd[,1] <- as.character(RAW_sAB[,1])
RAN_sAB_sd[,2] <- as.character(RAW_sAB[,2])

#---no: of levels for Bonferroni correction
Nlevels = 0;

for( i in 1:length(RAW_sAB[,1]) ){

    X = unlist(strsplit(RAW_sAB[i,1],":"))
    Y = unlist(strsplit(RAW_sAB[i,2],":"))

    ## select LITvDDD pairs
    #if( (X[1]!=Y[1]) && (X[2]==Y[2]) ){
    ## select LIT v ALL DDD pairs
    if( (X[1]!=Y[1]) ){        
        RAN_sAB_mean[i,3] <- as.character(mean(as.vector(unlist(RAW_sAB[i,3:perms]))))
        RAN_sAB_sd[i,3]   <- as.character(sd(as.vector(unlist(RAW_sAB[i,3:perms]))))
        Nlevels=Nlevels+1;
    }  

}

#--- Outfile for model-model separation/overlap
zsCN <-  c("ID","GENE.NAME","DIS","Model.long","Model","N","ID","GENE.NAME","DIS","Model.long","Model","N","actual.Overlap","sAB","Separated","Overlap","zScore","pvalue","Separation/Overlap.than.chance","Bonferroni","p.adjusted","q-value")
zs <- matrix(".", nrow=length(RAN_sAB_mean[,1]), ncol=length(zsCN))
colnames(zs) <- zsCN

for( i in 1:length(RAW_sAB[,1]) ){    
    
    X = unlist(strsplit(RAW_sAB[i,1],":"))
    Y = unlist(strsplit(RAW_sAB[i,2],":"))

    ## select only LITvDDD pairs
    #if( (!is.nan(as.numeric(RAN_sAB_sd[i,3]))) && (X[1]!=Y[1]) && (X[2]==Y[2]) ){

    ## select LIT v ALL DDD pairs
    if( (!is.nan(as.numeric(RAN_sAB_sd[i,3]))) && (X[1]!=Y[1]) ){

        isAB = which(colnames(sAB)==RAW_sAB[i,1])
        jsAB = which(rownames(sAB)==RAW_sAB[i,2])

        if( sAB[isAB,jsAB] == "." ){
            itmp = isAB
            isAB = jsAB
            jsAB = itemp 
        }

        zScore = 0
        
         #compute z-score, i.e. separation of mean sAB, against a randomised model (of the mean of sAB),   
         #see (Menche et al., 2015).
        if( as.numeric(RAN_sAB_sd[i,3]) != 0){
            zScore = (as.numeric(as.vector(sAB[isAB,jsAB])) - as.numeric(as.vector(RAN_sAB_mean[i,3])))/(as.numeric(as.vector(RAN_sAB_sd[i,3])))
        }


        #compute p.value from the normal distribution
        #See also http://www.cyclismo.org/tutorial/R/pValues.html     
        pval <- pnorm(-abs(zScore))
        pval <- 2 * pval

        zs[i,1]  <- as.character(pid[which(ptype==RAN_sAB_mean[i,1])])
        zs[i,2]  <- GN2[as.numeric(zs[i,1])+1]
        zs[i,3]  <- DIS2[as.numeric(zs[i,1])+1]
        zs[i,4]  <- as.character(phenl[which(ptype==RAN_sAB_mean[i,1])])    
        zs[i,4]  <- strSplit(zs[i,4])
        zs[i,5]  <- as.character(RAN_sAB_mean[i,1])
        zs[i,5]  <- strSplit(zs[i,5])
        zs[i,6]  <- as.character(loc[which(loc[,1]==RAN_sAB_mean[i,1]),2])

        zs[i,7]  <- as.character(pid[which(ptype==RAN_sAB_mean[i,2])])
        zs[i,8]  <- GN2[as.numeric(zs[i,7])+1]
        zs[i,9]  <- DIS2[as.numeric(zs[i,7])+1]
        zs[i,10] <- as.character(phenl[which(ptype==RAN_sAB_mean[i,2])])
        zs[i,10] <- strSplit(zs[i,10])
        zs[i,11] <- as.character(RAN_sAB_mean[i,2])
        zs[i,11] <- strSplit(zs[i,11])
        zs[i,12] <- as.character(loc[which(loc[,1]==RAN_sAB_mean[i,2]),2])
    
        #sAB, the model-model separation/overlap measure, on the Ontology
        zs[i,14] <- as.character(sAB[isAB,jsAB])

        #sAB > 0, implies separation
        zs[i,15] <- ifelse((as.numeric(zs[i,14]) > 0), "YES", ".")

        #sAB < 0, implies overlap
        zs[i,16] <- ifelse((as.numeric(zs[i,14]) < 0), "YES", ".")

        zs[i,17] <- as.character(zScore)
        zs[i,18] <- as.character(pval)

        #z-scores < 0 (>0), implies separation/overlap smaller (larger) than by chance  
        zs[i,19] <- ifelse((as.numeric(zs[i,17]) < 0), "Smaller", "larger")

        #Bonferroni correction for p.value
        temp <- "."
        for( x in 1:length(alpha) ){
            if(as.numeric(zs[i,18]) < as.numeric(alpha[x]/Nlevels)){ temp <- stars[x] }
        }

        zs[i,20] <- temp
        
    }#if  

}

#--- remove empty records
zs = zs[zs[,1]!=".",]


#--- add in the actual overlap
tmp1 = sapply(Ovlap[,2],function(x) strsplit(x,":")[[1]][2])
tmp1 = sapply(tmp1,function(x) strsplit(x,"x")[[1]][1])

tmp2 = sapply(Ovlap[,3],function(x) strsplit(x,":")[[1]][2])
tmp2 = sapply(tmp2,function(x) strsplit(x,"x")[[1]][1])

#zs[,13] = as.character(Ovlap[match(tmp,zs[,1]),4])

for( i in 1:length(zs[,1]) ){
    indxA = (as.numeric(zs[,1]))==as.numeric(tmp1[i])
    indxB = (as.numeric(zs[,7])-Ngenes)==as.numeric(tmp2[i])    
    zs[i,13] = as.character(Ovlap[indxA & indxB,4])
}


##--- FDR calculation ---##
# Note: In my study i've only 66 disease-disease pairs to consider when correcting for multiple hypothesis
#       corrections, and i focus on Bonerroni corrections @ 0.001, implying i expect to get <<1 false positive.
#       If you've > 100 pairs to consider, you'll need to use the FDR calculation below. 
#
##---                 ---##
#calFDR=0
calFDR=1

#p.adjusted
zs[,21] <- p.adjust(as.numeric(zs[,18]),method="BY")

#Else, will use WGCNA's qvalue calculation
zs[,22] <- WGCNA::qvalue(as.numeric(zs[,18]),smooth.df=Nlevels)$qvalue
    

if(calFDR){

    #--- pre-calculate each random z-score, using standard deviation (sd) obtained from each eddie job...
    #    you could use here the global sd obtained from all random permutations.

    # The Global random sd 
    #GLOBALsd <- mean(as.vector(unlist(RAW_sd)))
    
    tests <- matrix(".", nrow=length(RAW_sAB[,1]),ncol=(length(RAW_sAB[1,])-2))

    for( i in 1:length(RAW_sAB[,1]) ){    
    
    X = unlist(strsplit(RAW_sAB[i,1],":"))
    Y = unlist(strsplit(RAW_sAB[i,2],":"))

    ## select only LITvDDD pairs
    #if( (!is.nan(as.numeric(RAN_sAB_sd[i,3]))) && (X[1]!=Y[1]) && (X[2]==Y[2]) ){

    ## select LIT v ALL DDD pairs
    if( (!is.nan(as.numeric(RAN_sAB_sd[i,3]))) && (X[1]!=Y[1]) && (X[2]==Y[2]) ){
        test <- 0
        
        Mn   <- as.numeric(RAN_sAB_mean[i,3])
        Sd   <- as.numeric(RAN_sAB_sd[i,3])
            
        RsAB <- as.numeric(RAW_sAB[i,-c(1,2)])
            
        #store the random z-score 
        tests[i,] <- as.numeric((as.numeric(RsAB - Mn))/Sd)

    }

    }
      

    tests = tests[tests[,1]!=".",]
    
#--- save random zscores for plotting    
outfile <- file(sprintf("%s/random_zscores_per_gene.csv",resdir),"w");
write.table( tests, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);


}#if

outfile <- file(sprintf("%s/model_pairs_per_gene.csv",resdir),"w");
write.table( zs, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);
