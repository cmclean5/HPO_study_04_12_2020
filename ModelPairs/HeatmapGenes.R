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
        V1   <- as.vector(X[,2])
        V2   <- as.numeric(X[,5])
        V3   <- as.numeric(X[,3])        
        V4   <- as.numeric(X[,6])
        if( grepl("x",V1) ){ V1 = gsub("x","",V1); }
        RES  <- cbind(V1,V2,V3,V4)
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
FILES[1] <- "model_separation_per_gene"
#HPO_model_separation_per_gene.csv

#--- load results
ff <- read.table(sprintf("%s/%s_%s.csv",subDIR[S],subDIR[S],FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

df = SabExtract(SabXX=ff,HPA=HPA,SELECTION=2,REMOVEx=FALSE)
df = as.data.frame(df)
df[,3] = as.numeric(df[,3])
df[,4] = as.numeric(df[,4])
colnames(df) <- c("V1","V2","value","V1.SIZE","V2.SIZE","overlap")

map  = cbind(unique(df$V1),lab)
map  = rbind(map,cbind(unique(df$V2),lab))

df$V1 = map[match(df$V1,map[,1]),2]
df$V2 = map[match(df$V2,map[,1]),2]

df$V1      <- as.vector(factor(df$V1))
df$V2      <- as.vector(factor(df$V2))
df$value   <- as.vector(factor(df$value))
df$V1.SIZE <- as.vector(factor(df$V1.SIZE))
df$V2.SIZE <- as.vector(factor(df$V2.SIZE))
df$overlap <- as.vector(factor(df$overlap))
df$V2 = factor(df$V2, levels=unique(df$V2))
#df$V2 = factor(df$V2, levels=rev(levels(factor(df$V2))))
df$V1 = factor(df$V1, levels=unique(df$V1))
df$V1 = factor(df$V1, levels=rev(levels(factor(df$V1))))

df2 = df
#df2[,1] = sprintf("%sx",df$V1)
#df2[,2] = sprintf("%sx",df$V2)

MODELS   = as.vector(unique(df2$V1))
MODELS.N = as.vector(df2$V1.SIZE[match(lab,df2$V1)])
N      = length(MODELS)


CN <- rep("",4*N)
oo <- matrix("",nrow=N,ncol=length(CN))

k1=1
k2=2
k3=3
k4=4
for( i in 1:N ){
    if( grepl("x",MODELS[i]) ){ V1 = gsub("x","",MODELS[i]) }
    else { V1 = MODELS[i]; }
    CN[k1] = V1
    CN[k2] = V1
    CN[k3] = V1
    CN[k4] = V1
    oo[,k1:k4] = rankSab(DF=df2,MODEL=MODELS[i]) 
    k1=k1+4
    k2=k2+4
    k3=k3+4
    k4=k4+4
}

colnames(oo) = CN

write.table(oo,"LITvDDDmodels_Sab.csv",sep="\t",row.names=F,col.names=T,quote=F)

CNrr = c("LIT.MODEL","LIT.SIZE","DIS","matching.Sab.score","matching.DDD.SIZE","matching.DDD.rank","actual.overlap","actual.overlap.rank","highest.DDD.SIZE","lowest.Sab.score","highest.DDD.model","actual.overlap","actual.overlap.rank","Overlap.Separate")
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
    Rindx   = which(X==MODELS[i])
    rr[i,4] = Sab[Rindx]
    rr[i,5] = X.N[Rindx]
    rr[i,6] = Rindx
    rr[i,7] = OV[Rindx]
    rr[i,8] = which(OV[Rindx]==OV[order(as.numeric(OV),decreasing=T)])[1]    

    rr[i,9]  = X.N[1]
    rr[i,10] = Sab[1]
    rr[i,11] = X[1]
    rr[i,12] = OV[1]
    rr[i,13] = which(OV[1]==OV[order(as.numeric(OV),decreasing=T)])[1]    

    if( as.numeric(Sab[1]) < 0 ){
        rr[i,14] = "Overlap"
    } else {
        rr[i,14] = "Separate"
    }
    
}

rr = rr[order(as.numeric(rr[,6])),]

write.table(rr,"modelRanking.csv",sep="\t",row.names=F,col.names=T,quote=F)

#--------------

top = (max(as.numeric(as.vector(df$value))))
bot = (min(as.numeric(as.vector(df$value))))
#mid = ((top-bot)/2) - abs(bot)

top = abs(bot)
bot = bot
mid = 0     

BORDERCOL= "grey50"#"black"#grey60
NACOL    = "red"#grey90"

 # Create a ggheatmap
    gplot <- ggplot(df, aes(df$V2, df$V1, fill = as.numeric(df$value)))+
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

        str = sprintf("%s/HeatMap%s.png",plotdir, FILES[1])    
    
        png(str,width=WIDTH,height=HEIGHT,units="px");
        print(gplot)
        dev.off()


