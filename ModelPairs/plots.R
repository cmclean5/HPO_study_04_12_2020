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


#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

plotdir <- sprintf("PLOTS/%s",subDIR[S])

if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}


run1=1

if( run1 ){

FILES    <- vector(length=2)
FILES[1] <- "model_overlap_sig"
FILES[2] <- "random_zscores"

#--- load results
ff <- read.table(sprintf("RESULTS/%s/%s.csv",subDIR[S],FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

#--- load randomised zscores
tests <- read.table(sprintf("RESULTS/%s/%s.csv",subDIR[S],FILES[2]),sep="\t",header=F,stringsAsFactors=F,quote="", check.names=F)

#--set disease pair of interest
phenA = vector(length=1)
phenA[1] = "LIT"

phenB = vector(length=1)
phenB[1] = "DDD"
    
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

