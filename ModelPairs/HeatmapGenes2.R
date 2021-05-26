#run 'mergeRESULTS.R' first

source('../setUp.R')
library(corrplot)
#require(cowplot)

actualOverlap <- function(HPA,dddX,litX){

    ov <- c(0,0,0)
    ov[1] = sum(grepl(litX,HPA))
    ov[2] = sum(grepl(dddX,HPA))
    ov[3] = sum(grepl(litX,HPA) & grepl(dddX,HPA))

    return(ov)
    
}


SabExtract <- function(SabXX, HPA, ROUND=TRUE, SELECTION=1,REMOVEx=TRUE){

colnames(SabXX) <- gsub("\\.",":",colnames(SabXX))

N = length(colnames(SabXX))

    yy <- data.frame(a=as.character(),b=as.character(),c=as.character(),
                     d=as.character(),e=as.character(),f=as.character())
    
# calculate model-model overlap/separation
for( i in 1:N ){
    for( j in 1:N ){

        if( i < j ){

            ov = actualOverlap(HPA,litX=rownames(SabXX)[i],dddX=colnames(SabXX)[j])
            
            X = unlist(strsplit(rownames(SabXX)[i],":"))
            Y = unlist(strsplit(colnames(SabXX)[j],":"))

            if( REMOVEx ){
                model.A = sprintf("%s:%s",X[1],gsub("x","",X[2]))
                model.B = sprintf("%s:%s",Y[1],gsub("x","",Y[2]))
            } else {
                model.A = rownames(SabXX)[i]
                model.B = colnames(SabXX)[j]
            }
                
            if( SELECTION == 1 ){ 
        
                #select only LITvDDD pairs
                if( (i!=j) && (X[1]!=Y[1]) && (X[2]==Y[2]) ){       
            
                    tmp = data.frame(a=model.A,b=model.B,
                                     c=as.character(SabXX[i,j]),
                                     d=as.character(ov[1]),
                                     e=as.character(ov[2]),
                                     f=as.character(ov[3]))
                    yy = rbind(yy,tmp)
                }
            }#1

            if( SELECTION == 2 ){ 
        
                #select LIT V. any DDD pairs
                if( (i!=j) && (X[1]!=Y[1]) ){       
            
                    tmp = data.frame(a=model.A,b=model.B,
                                     c=as.character(SabXX[i,j]),
                                     d=as.character(ov[1]),
                                     e=as.character(ov[2]),
                                     f=as.character(ov[3]))
                    yy = rbind(yy,tmp)
                }
            }#2        
            
        }#i<j
    }
}

    if( ROUND ){ 
        yy[,3] = as.character(round(as.numeric(as.vector(yy[,3])),5))
    }
        
    return(yy)
    
}

rankSab <- function(DF,MODEL=NULL){

    if( !is.null(MODEL) ){

        indx <- grepl(MODEL,DF[,1])
        X    <- DF[indx,]
        X    <- X[order(as.numeric(X[,3]),decreasing=F),]
        V1   <- as.vector(X[,2])  ##ddd.model.name
        V2   <- as.numeric(X[,8]) ##ddd.model.size
        V3   <- as.numeric(X[,3]) ##sAB       
        V4   <- as.numeric(X[,9]) ##actual overlap
        V5   <- as.numeric(X[,4]) ##pvalue
        V6   <- as.numeric(X[,6]) ##qscore
        if( grepl("x",V1) ){ V1 = gsub("x","",V1); }
        RES  <- cbind(V1,V2,V3,V4,V5,V6)
    }

    return(RES)
    
}

#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]


#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

plotdir <- sprintf("PLOTS/%s",subDIR[S])

if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}

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
DIS  <- ANNO[match(ANNO[,1],pid[1:length(ANNO[,2])]),3]
Nlab <- seq(0,(length(GN)-1),1)
lab  <- sprintf("%s_%d",GN,Nlab)


FILES    <- vector(length=1)
FILES[1] <- "model_pairs_per_gene"

#--- load results
ff <- read.table(sprintf("RESULTS/%s/%s.csv",subDIR[S],FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

    df = cbind(sprintf("%s_%d",ff[,2],ff[,1]),sprintf("%s_%d",ff[,8],(as.numeric(ff[,7])-Ngenes)),
               ff[,14],ff[,18],ff[,21],ff[,22],ff[,6],ff[,12],ff[,13])
    df = as.data.frame(df)    
    colnames(df) <- c("V1","V2","sAB","pval","adj","qscore","V1.SIZE","V2.SIZE","OVLAP")

    df$V1      <- as.vector(factor(df$V1))
    df$V2      <- as.vector(factor(df$V2))
    df$sAB     <- as.vector(factor(df$sAB))
    df$pval    <- as.vector(factor(df$pval))
    df$adj     <- as.vector(factor(df$adj))
    df$qscore  <- as.vector(factor(df$qscore))
    df$V1.SIZE <- as.vector(factor(df$V1.SIZE))
    df$V2.SIZE <- as.vector(factor(df$V2.SIZE))
    df$OVLAP   <- as.vector(factor(df$OVLAP))

    df$V2 = factor(df$V2, levels=unique(df$V2))
    #df$V2 = factor(df$V2, levels=rev(levels(factor(df$V2))))
    df$V1 = factor(df$V1, levels=unique(df$V1))
    df$V1 = factor(df$V1, levels=rev(levels(factor(df$V1))))


    df2 = df
    MODELS   = as.vector(unique(df2$V1))
    MODELS.N = as.vector(df2$V1.SIZE[match(lab,df2$V1)])
    N        = length(MODELS)

    CN <- rep("",6*N)
    oo <- matrix("",nrow=N,ncol=length(CN))

    k1=1
    k2=2
    k3=3
    k4=4
    k5=5
    k6=6
for( i in 1:N ){
    if( grepl("x",MODELS[i]) ){ V1 = gsub("x","",MODELS[i]) }
    else { V1 = MODELS[i]; }
    CN[k1] = V1
    CN[k2] = V1
    CN[k3] = V1
    CN[k4] = V1
    CN[k5] = V1
    CN[k6] = V1
    oo[,k1:k6] = rankSab(DF=df2,MODEL=MODELS[i]) 
    k1=k1+6
    k2=k2+6
    k3=k3+6
    k4=k4+6
    k5=k5+6
    k6=k6+6
}

colnames(oo) = CN

#write.table(oo,"LITvDDDmodels_Sab.csv",sep="\t",row.names=F,col.names=T,quote=F)


CNrr = c("LIT.MODEL","LIT.SIZE","DIS","matching.Sab.score","matching.pvalue","matching.qscore","matching.DDD.SIZE","matching.DDD.rank","actual.overlap","actual.overlap.rank","highest.DDD.SIZE","lowest.Sab.score","pvalue","qscore","highest.DDD.model","actual.overlap","actual.overlap.rank","Overlap.Separate")
rr     = matrix("",ncol=length(CNrr),nrow=N)
colnames(rr) = CNrr
rr[,1] = MODELS
rr[,2] = MODELS.N
rr[,3] = DIS

for( i in 1:N ){
    Cindx   = grep(MODELS[i],colnames(oo))
    X       = oo[,Cindx[1]]
    X.N     = oo[,Cindx[2]]
    Sab     = oo[,Cindx[3]]
    OV      = oo[,Cindx[4]]
    PV      = oo[,Cindx[5]]
    QS      = oo[,Cindx[6]]
    Rindx   = which(X==MODELS[i])
    rr[i,4] = Sab[Rindx]
    rr[i,5] = PV[Rindx]
    rr[i,6] = QS[Rindx]
    rr[i,7] = X.N[Rindx]
    rr[i,8] = Rindx      ##order
    rr[i,9] = OV[Rindx]
    rr[i,10] = which(OV[Rindx]==OV[order(as.numeric(OV),decreasing=T)])[1]    

    rr[i,11] = X.N[1]
    rr[i,12] = Sab[1]
    rr[i,13] = PV[1]
    rr[i,14] = QS[1]
    rr[i,15] = X[1]
    rr[i,16] = OV[1]
    rr[i,17] = which(OV[1]==OV[order(as.numeric(OV),decreasing=T)])[1]    

    if( as.numeric(Sab[1]) < 0 ){
        rr[i,18] = "Overlap"
    } else {
        rr[i,18] = "Separate"
    }
    
}

rr = rr[order(as.numeric(rr[,8])),]

write.table(rr,"modelRanking_sig.csv",sep="\t",row.names=F,col.names=T,quote=F)

#--------------

    
    top = (max(as.numeric(as.vector(df$sAB))))
    bot = (min(as.numeric(as.vector(df$sAB))))
    #mid = ((top-bot)/2) - abs(bot)

    top = abs(bot)
    bot = bot
    mid = 0     

    df$sAB = ifelse(as.numeric(as.vector(df$sAB))>top,top,as.numeric(as.vector(df$sAB)))
    #df$sAB = ifelse(as.numeric(as.vector(df$sAB))<bot,bot,as.numeric(as.vector(df$sAB)))
    
    #sAB score which are no significant
    #df$sAB = ifelse(as.numeric(as.vector(df$pval)) < 0.05,as.numeric(as.vector(df$sAB)),NA)
    df$sAB = ifelse(as.numeric(as.vector(df$qscore)) < 0.05,as.numeric(as.vector(df$sAB)),NA)
    #df$sAB = ifelse(as.numeric(as.vector(df$adj)) < 0.05,as.numeric(as.vector(df$sAB)),NA)
    #------------

    BORDERCOL= "grey50"#"black"#grey60
    NACOL    = "grey30"#black"#grey90"

    # Create a ggheatmap
    gplot <- ggplot(df, aes(df$V2, df$V1, fill = as.numeric(df$sAB)))+
        labs(x="",y="",title="")+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid="white", 
                             midpoint = mid, limit = c(bot,top), space = "Lab", 
                             name="Sab",
                             na.value=NACOL) +
        #theme_minimal()+ # minimal theme
        theme_grey(base_size=10)+
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
        coord_fixed()+        

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(angle = 40, hjust = 1, face="bold", size=rel(0.6)),
        axis.text.y  = element_text(face="bold",size=rel(1)),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = BORDERCOL, fill=NA, size=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = BORDERCOL),
        axis.ticks = element_blank(),
        legend.justification = c("right","center"),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.position = "right",
        legend.direction = "vertical")+
    
        guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                     title.position = "top", title.hjust = 0.5))

        str = sprintf("%s/HeatMap%s_sig.png",plotdir, FILES[1])    
    
        png(str,width=WIDTH,height=HEIGHT,units="px");
        print(gplot)
        dev.off()


