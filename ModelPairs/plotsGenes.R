#run 'mergeRESULTS.R' first

source('../setUp.R')
#require(cowplot)

makeBold <- function(src, bolder) {
    #require(purrr)
    if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
    src_levels <- levels(src)                                 # retrieve the levels in their order
    temp <- bolder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
    if (all(temp)) {                                         # if so
        b_pos <- purrr::map_int(bolder, ~which(.==src_levels)) # then find out where they are
        b_vec <- rep("plain", length(src_levels))               # make'm all plain first
        b_vec[b_pos] <- "bold"                                  # make our targets bold
        b_vec                                                   # return the new vector
    } else {
        stop("All elements of 'bolder' must be in src")
    }
}


makePairs <- function(DF){

    TB = table(DF)

    pairs = rep("",length=length(DF))
    
    for( i in 1:length(TB) ){

        indx = which(DF==names(TB)[i])

        k=indx[1]
        for( j in 1:length(indx)){            
            pairs[k] = paste(sprintf("%s_%s",names(TB)[i],letters[j]))            
            k=k+1
        }
    }

    return(pairs)
    
}
    
#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

plotdir <- sprintf("PLOTS/%s",subDIR[S])

if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}


FILES    <- vector(length=2)
FILES[1] <- "model_pairs_per_gene"
FILES[2] <- "random_zscores_per_gene"

#--- load results
ff <- read.table(sprintf("RESULTS/%s/%s.csv",subDIR[S],FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

#--- load randomised zscores
tests <- read.table(sprintf("RESULTS/%s/%s.csv",subDIR[S],FILES[2]),sep="\t",header=F,stringsAsFactors=F,quote="", check.names=F)

#--set disease pair of interest
phenA = vector(length=1)
phenA[1] = "LIT"

phenB = vector(length=1)
phenB[1] = "DDD"

run1=0

if( run1 ){
    
GPLOTS <- list()
g=1

for( i in 1:length(phenA) ){

    pA = phenA[i]
    pB = phenB[i]
     
    #--- locate disease pair in results file
    indx = (ff[,3] == pA & ff[,7] == pB) | (ff[,3] == pB & ff[,7] == pA)

    #--- Plot 1: Z-score of disA versus disB for 1000 randomised studies on the network model.
    if( sum(indx) == 1 ){

        obsZs <- as.numeric(ff[indx,12])    
        DF    <- as.numeric(as.vector(tests[indx,]))
        DF    <- as.data.frame(DF)
        names(DF) <- "V1"

        MIN = sign(min(DF)) * (abs(min(DF))+1)
        MAX = sign(max(DF)) * (abs(max(DF))+1)

        if( abs(MIN) < 5 ){
            MIN = -5
        }

        if( abs(MAX) < 5 ){
            MAX = 5
        }
        
        if( obsZs <= MIN ){
            MIN = sign(obsZs) * (abs(obsZs)+1)
        }

        if( obsZs >= MAX ){
            MAX = sign(obsZs) * (abs(obsZs)+1)
        }
       
        
        BINS=50
        ARROW_Y_HEIGHT=350
        
        tit = sprintf("%sv%s",pA,pB)

        annoT = grobTree(textGrob(sprintf("z-score = %.3f",round(as.numeric(obsZs),4)),
                                  x=0.05, y=0.35, hjust=0,
                                  gp=gpar(col="black", fontface="bold", fontsize=20)))
    
        gplot1 <- ggplot(DF,aes(x=DF$V1))+geom_histogram(bins=BINS,colour="brown1",fill="cadetblue") +
            annotation_custom(annoT)+
            labs(x=c(expression(paste(" ",S[AB]))),y="Count",title=sprintf("%s Vs %s",as.character(gsub("_"," ",ff[indx,2])),as.character(gsub("_"," ",ff[indx,6]))))+
            ##labs(x=c(expression(paste("Separation ",S[AB]))),y="Count",title=sprintf("%s Vs %s",as.character(gsub("_"," ",ff[indx,2])),as.character(gsub("_"," ",ff[indx,6]))))+
            scale_x_continuous(limit = c(MIN, MAX))+
            theme(
                title=element_text(face="bold",size=rel(1.5)),
                axis.title.x=element_text(face="bold",size=rel(1.5)),
                axis.title.y=element_text(face="bold",size=rel(1.5)),
                legend.title=element_text(face="bold",size=rel(2.0)),
                legend.text=element_text(face="bold",size=rel(2.0)),
                legend.key=element_blank())+
            theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                  panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill = "white"),
                  panel.border = element_rect(linetype="solid",fill=NA))+
            geom_segment(data=data.frame(obsZs), aes(x=obsZs,y=ARROW_Y_HEIGHT,xend=obsZs,yend=0.1),size=2,colour="orange2",arrow=arrow(type = "closed"))
        
        png(sprintf("%s/%s_%s.png",plotdir, tit, subDIR[S]),width=WIDTH,height=HEIGHT,units="px")
        print(gplot1)
        dev.off()

        GPLOTS[[g]]      <- gplot1
        names(GPLOTS)[g] <- LETTERS[g]
        g=g+1
        
    }
}

p <- plot_grid(plotlist=GPLOTS, labels = "AUTO", label_size=50, label_fontface="bold")
ggsave(sprintf("%s/Plots.png",plotdir), p, width=20, height=20, device="png")

}


run2=1

if( run2 ){


sabIDX <- which(grepl("sAB",colnames(ff),fixed=T)==TRUE)

zsIDX  <- which(grepl("zScore",colnames(ff),fixed=T)==TRUE)

BonIDX <- which(grepl("Bonferroni",colnames(ff),fixed=T)==TRUE)

pvIDX  <- which(grepl("pvalue",colnames(ff),fixed=T)==TRUE)

adIDX  <- which(grepl("p.adjusted",colnames(ff),fixed=T)==TRUE)
    
qvIDX  <- which(grepl("q-value",colnames(ff),fixed=T)==TRUE)

DDpairs <- sprintf("%s_%g",ff[,2],ff[,1])
##DDpairs <- makePairs(ff[,2])

nIDX <- which(colnames(ff)=="N")

acIDX  <- which(grepl("actual.Overlap",colnames(ff),fixed=T)==TRUE)
    
oo     <- matrix(0, ncol=9, nrow=length(DDpairs) )
colnames(oo) <- c("pairs","sAB","zscore","pvalue","p.adjusted","qvalue","Nmin","overlap","label")
oo[,1] <- as.character(DDpairs)
oo[,2] <- as.character(ff[,sabIDX[1]])
oo[,3] <- as.character(ff[,zsIDX[1]])
oo[,4] <- as.character(ff[,pvIDX[1]])
oo[,5] <- as.character(ff[,adIDX[1]])
oo[,6] <- as.character(ff[,qvIDX[1]])
oo[,7] <- as.character(apply(ff[,nIDX],1,min))
oo[,8] <- as.character(ff[,acIDX[1]])
oo[,9] <- rep("",length(DDpairs))

#oo <- oo[oo[,2] != 0 & oo[,3] != 0,]
    
df <- as.data.frame(oo)

##overlap
df1 <- df[as.numeric(df$zscore) < 0,]

ReOrder = order(as.numeric(df1$pvalue),decreasing=T)
df1 = df1[ReOrder,]   
df1$pairs      <- as.vector(factor(df1$pairs))
df1$sAB        <- as.vector(factor(df1$sAB))
df1$zscore     <- as.vector(factor(df1$zscore))
df1$pvalue     <- as.vector(factor(df1$pvalue))
df1$p.adjusted <- as.vector(factor(df1$p.adjusted))
df1$qvalue     <- as.vector(factor(df1$qvalue))
df1$Nmin       <- as.vector(factor(df1$Nmin))
df1$overlap    <- as.vector(factor(df1$overlap))
df1$pairs      <- factor(df1$pairs, levels=df1$pairs)
    
write.table(df1,sprintf("%s_overlap.csv", subDIR[S]), sep="\t", col.names=T, row.names=F, quote=F)

##separation
df2 <- df[as.numeric(df$zscore) >= 0,]

ReOrder = order(as.numeric(df2$pvalue),decreasing=T)
df2 = df2[ReOrder,]   
df2$pairs      <- as.vector(factor(df2$pairs))
df2$sAB        <- as.vector(factor(df2$sAB))
df2$zscore     <- as.vector(factor(df2$zscore))
df2$pvalue     <- as.vector(factor(df2$pvalue))
df2$p.adjusted <- as.vector(factor(df2$p.adjusted))
df2$qvalue     <- as.vector(factor(df2$qvalue))
df2$Nmin       <- as.vector(factor(df2$Nmin))
df2$overlap    <- as.vector(factor(df2$overlap))
df2$pairs      <- factor(df2$pairs, levels=df2$pairs)
    
write.table(df2,sprintf("%s_separation.csv", subDIR[S]), sep="\t", col.names=T, row.names=F, quote=F)

  
}


##overlap
run3 = 1

if( run3 ){

    df1 = read.table(sprintf("%s_overlap.csv",subDIR[S]),sep="\t", header=T)
    df1$label = "Overlap"   

    df2 = read.table(sprintf("%s_separation.csv",subDIR[S]),sep="\t", header=T)
    df2$label = "Separation"  
    
    groups <- c("Over","Sep")   

    if( length(groups) == 1 ){
        DF = rbind(df2[,c(1,4,5,6,7,8,9)])
        df = DF
        colours <- c('royalblue2')
        shapes  <- c(17)
    }   
    
    if( length(groups) == 2 ){        
        DF = rbind(df1[,c(1,4,5,6,7,8,9)],df2[,c(1,4,5,6,7,8,9)])
        df = DF
        colours <- c('firebrick2','royalblue2')
        shapes  <- c(16,17)
    }   
     

    #---scale Nmin
    df$Nmin = 10*(df$Nmin/max(df$Nmin))

    #---scale actual.overlap
    Ovmin      = min(as.numeric(df$overlap))
    Ovmax      = max(as.numeric(df$overlap))
    df$overlap = 10*(as.numeric(df$overlap)-Ovmin)/(Ovmax-Ovmin)
    
    #Highlight Epi-XX pairs
    #target <- c("SCH-Epi","ASD-Epi","Epi-MS","BD-Epi")
    #target <- c("AD-HTN","HTN-PD","AD-PD")   

    plotSize=0
    plotNmin=0
    plotOvlp=1
    
    #plot p.adjusted
    ReOrder = order(as.numeric(df$p.adjusted),decreasing=T)
    df = df[ReOrder,]
    df$pairs      <- as.vector(factor(df$pairs))
    df$pvalue     <- as.vector(factor(df$pvalue))
    df$p.adjusted <- as.vector(factor(df$p.adjusted))
    df$qvalue     <- as.vector(factor(df$qvalue))
    df$label      <- as.vector(factor(df$label))
    df$overlap    <- as.vector(factor(df$overlap))
    df$Nmin       <- as.vector(factor(df$Nmin))
    df$pairs = factor(df$pairs, levels=unique(df$pairs))

    XMIN= 0
    XMAX=-log10(min(as.numeric(df$p.adjusted)))

    gplot <- ggplot(df,aes(y=as.factor(df$pairs),x=-log10(as.numeric(as.vector(df$p.adjusted))),group=df$label) )+
    {if(plotOvlp)geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(as.numeric(df$overlap)))}+
    {if(plotNmin)geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(as.numeric(df$Nmin)))}+
        {if(plotSize)geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(3))}+
        scale_size_identity()+
        labs(y="LITvDDD per gene",x="-log10(p.adjusted)")+
        theme(legend.key=element_blank())+
        geom_vline(xintercept=-log10(0.05),colour="grey25",linetype="longdash", alpha=0.5, size=rel(1),show.legend=F)+
        theme_bw()+
        scale_x_continuous(limit = c(XMIN, XMAX))+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"))+
        theme(axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.5)),
              #axis.text.y = element_text(face=makeBold(df$pairs,target),size=rel(1.5)),
              axis.text.y = element_text(face="bold",size=rel(1.5)),
              axis.title.x=element_text(face="bold",size=rel(2.0)),
              axis.title.y=element_text(face="bold",size=rel(2.0)),
              legend.title=element_text(face="bold",size=rel(1.7)),
              legend.text=element_text(face="bold",size=rel(1.7)),
              legend.position="top"
              )+
        guides(colour = guide_legend(override.aes = list(shape = shapes, size=rel(7))),
               size   = FALSE,
               shape  = FALSE)+
        scale_color_manual("",breaks=levels(factor(df$label)),values=c(colours))

    #ggsave(gplot,file=sprintf("%s/%s_scale.svg",pltdir,alg[1]), width=WIDTHcm, height=2*HEIGHTcm, units="cm");

    str = sprintf("%s/%s_model_pairs_overlaps_adj.png",plotdir, subDIR[S])    
    if(plotNmin){
        str = sprintf("%s/%s_model_pairs_overlapsNmin_adj.png",plotdir, subDIR[S])        
    }
    if(plotOvlp){
        str = sprintf("%s/%s_model_pairs_overlapsOvlp_adj.png",plotdir, subDIR[S])        
    }
    
    png(str,width=WIDTH,height=2*HEIGHT,units="px");
    print(gplot)
    dev.off()

    #plot qvalue
    ReOrder = order(as.numeric(df$qvalue),decreasing=T)
    df = df[ReOrder,]
    df$pairs      <- as.vector(factor(df$pairs))
    df$pvalue     <- as.vector(factor(df$pvalue))
    df$p.adjusted <- as.vector(factor(df$p.adjusted))
    df$qvalue     <- as.vector(factor(df$qvalue))
    df$label      <- as.vector(factor(df$label))
    df$overlap    <- as.vector(factor(df$overlap))
    df$Nmin       <- as.vector(factor(df$Nmin))
    df$pairs = factor(df$pairs, levels=unique(df$pairs))

    XMIN= 0
    XMAX=-log10(min(as.numeric(df$qvalue)))     

    
    
    gplot <- ggplot(df,aes(y=as.factor(df$pairs),x=-log10(as.numeric(as.vector(df$qvalue))),group=df$label) )+
        {if(plotOvlp)geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(as.numeric(df$overlap)))}+
        {if(plotNmin)geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(as.numeric(df$Nmin)))}+
        {if(plotSize)geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(3))}+
        scale_size_identity()+
        labs(y="LITvDDD per gene",x="-log10(q-value)")+
        theme(legend.key=element_blank())+
        geom_vline(xintercept=-log10(0.05),colour="grey25",linetype="longdash", alpha=0.5, size=rel(1),show.legend=F)+
        theme_bw()+
        scale_x_continuous(limit = c(XMIN, XMAX))+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"))+
        theme(axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.5)),
              #axis.text.y = element_text(face=makeBold(df$pairs,target),size=rel(1.5)),
              axis.text.y = element_text(face="bold",size=rel(1.5)),
              axis.title.x=element_text(face="bold",size=rel(2.0)),
              axis.title.y=element_text(face="bold",size=rel(2.0)),
              legend.title=element_text(face="bold",size=rel(1.7)),
              legend.text=element_text(face="bold",size=rel(1.7)),
              legend.position="top"
              )+
        guides(colour = guide_legend(override.aes = list(shape = shapes, size=rel(7))),
               size   = FALSE,
               shape  = FALSE)+
        scale_color_manual("",breaks=levels(factor(df$label)),values=c(colours))

    #ggsave(gplot,file=sprintf("%s/%s_scale.svg",pltdir,alg[1]), width=WIDTHcm, height=2*HEIGHTcm, units="cm");

    str = sprintf("%s/%s_model_pairs_overlaps.png",plotdir, subDIR[S])    
    if(plotNmin){
        str = sprintf("%s/%s_model_pairs_overlapsNmin.png",plotdir, subDIR[S])        
    }
    if(plotOvlp){
        str = sprintf("%s/%s_model_pairs_overlapsOvlp.png",plotdir, subDIR[S])        
    }
    
    png(str,width=WIDTH,height=2*HEIGHT,units="px");
    print(gplot)
    dev.off()
    
}

