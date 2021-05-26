
logBase <- function(x, BASE=0){

    x = as.numeric(x)

    if( BASE == 2 ){
        return(log2(x))
    } else if ( BASE == 10 ){
        return(log10(x))
    } 

    return(log(x))

}

entropy <- function(x, BASE=0){
    x   <- as.numeric(x)    
    if( is.na(x) ){ return(NA) }
    else if( x == 0 ){ return(0) }
    else { return(-x*logBase(x,BASE)) }    
}

## KL divergence KL(P || Q)
#kl_divergence <- function(p,q){    
#    p = as.numeric(p)
#    q = as.numeric(q)
#    if( is.na(p) || is.na(q) ){return(NA)}
#    else if(p==0 || q==0){ return(NA) }
#    else { return(p*log(p/q)) }
#}

## KL divergence KL(P || Q)
kl_divergence <- function(p,q){
    p = as.numeric(p)
    q = as.numeric(q)
    if( p==0 ){ p=p+SMALL; }
    if( q==0 ){ q=q+SMALL; }    
    if( is.na(p) || is.na(q) ){return(NA)}
    else { return(p*log(p/q)) }
}


## Jensen-Shannon divergence JS(P || Q)
js_divergence <- function(p,q){
    p = as.numeric(p)
    q = as.numeric(q)
    m = 0.5*(p+q)
    return( (0.5*(kl_divergence(p,m)) + 0.5*(kl_divergence(q,m))) )    
}

## Jensen-Shannon divergence JS (P || Q)
## See eqn 4, https://arxiv.org/pdf/1503.03512.pdf
js_divergence2 <- function(p,q){
    p = as.numeric(p)
    q = as.numeric(q)
    m = 0.5*(p+q);
    r = min(p,q)/m;
    if( is.na(m) || is.na(r) ){ return(NA) }
    term1 = ifelse( r==0,     0, r*log(r));
    term2 = ifelse( (2-r)==0, 0, (2-r)*log(2-r));
    return( m*0.5*(term1+term2) );
}

rough.fixNA <- function(ONTO,ANNO.LAB="",NEW.ANNO.LAB="",type=c("median","mean")){
## if data is continuous, use column median (MED==TRUE),
## or mean (MED==FALSE).
## Ref: https://cran.r-project.org/web/packages/RRF/RRF.pdf

    type =  match.arg(type);
    
    label.indx = which(names(ONTO)==ANNO.LAB)

    if( length(label.indx) != 0 ){

        if( NEW.ANNO.LAB == "" ){
            NEW.ANNO.LAB   = "temp"
            Nlabels = length(names(ONTO));
            new.label.indx = Nlabels+1;
        }       
        else {
            new.label.indx = which(names(ONTO)==NEW.ANNO.LAB);
            if( length(new.label.indx) != 0 ){
                ONTO       = rm.onto.annotation(ONTO,NEW.ANNO.LAB);
                level.indx = which(names(ONTO)==LEVEL.LAB)
                label.indx = which(names(ONTO)==ANNO.LAB)
            }
            Nlabels = length(names(ONTO));
            new.label.indx = Nlabels+1;
        }#else
        
        val = as.numeric(ONTO[[label.indx]])        
    
        if( type == "median" ){
            med  = median(val,na.rm=T)
            freq = ifelse(is.na(val),med,val)  
        } else {
            mean = mean(val,na.rm=T)
            freq = ifelse(is.na(val),mean,val)
        }

        ONTO[[new.label.indx]] = freq
        if( NEW.ANNO.LAB != "" ){
            names(ONTO)[new.label.indx] = NEW.ANNO.LAB; 
        }
        
    }#if
        
    return(ONTO)
    
}

## use median, or mean, of annotations on a given level to impute missing HP annotation
## values.
level.fixNA <- function(ONTO,LEVEL.LAB="levels",ANNO.LAB="", NEW.ANNO.LAB="",type=c("median","mean")){
    
    level.indx = which(names(ONTO)==LEVEL.LAB)
    label.indx = which(names(ONTO)==ANNO.LAB)

    if( length(level.indx) != 0 && length(label.indx) != 0 ){

        if( NEW.ANNO.LAB == "" ){
            NEW.ANNO.LAB   = "temp"
            Nlabels = length(names(ONTO));
            new.label.indx = Nlabels+1;
            ##new.label.indx = label.indx;
        }       
        else {
            new.label.indx = which(names(ONTO)==NEW.ANNO.LAB);
            if( length(new.label.indx) != 0 ){
                ONTO       = rm.onto.annotation(ONTO,NEW.ANNO.LAB);
                level.indx = which(names(ONTO)==LEVEL.LAB)
                label.indx = which(names(ONTO)==ANNO.LAB)
            }
            Nlabels = length(names(ONTO));
            new.label.indx = Nlabels+1;
        }#else
        
        df = cbind(ONTO[[level.indx]],ONTO[[label.indx]])

        Nlevels = sum(!is.na(unique(df[,1])))
        
        levels = matrix(NA,nrow=Nlevels,ncol=2)
        levels[,1] = seq(0,(Nlevels-1),1)

        for( l in 0:(Nlevels-1) ){
            if( type == "median" ){
                levels[(l+1),2] = median(df[df[,1]==l,2],na.rm=T)
            } else {
                #levels[(l+1),2] = ceiling(mean(df[df[,1]==l,2],na.rm=T))
                levels[(l+1),2] = mean(df[df[,1]==l,2],na.rm=T)
            }
        }

        if( type == "median" ){
            NA.level = median(levels[,2],na.rm=T)
        } else {
            #NA.level = ceiling(mean(levels[,2],na.rm=T))
            NA.level = mean(levels[,2],na.rm=T)
        }
        freq     = levels[match(df[,1],levels[,1]),2]
        freq     = ifelse(is.na(freq),NA.level,freq)
        freq     = ifelse(is.na(df[,2]),freq,df[,2])
        
        ONTO[[new.label.indx]] = freq
        if( NEW.ANNO.LAB != "" ){
            names(ONTO)[new.label.indx] = NEW.ANNO.LAB; 
        }
        
    }#if

    return(ONTO)
    
}

add.hp.clip <- function(ONTO,CLIP){

    indx = match(ONTO$id,CLIP)
    indx = ifelse(is.na(indx),0,1)
    ONTO$clip = indx
    names(ONTO$clip) = ONTO$id

    return(ONTO)
}

get.hp.subsumers <- function(ONTO){

    ONTO$Nsubsumers = rep(NA,length(ONTO$id))
    names(ONTO$Nsubsumers) = ONTO$id
    
    for( i in 1:length(ONTO$id) ){

        if( SUB ){
            anc = get.sub.acendants( ONTO, ONTO$id[i] );
        } else {
            anc = get.descendants( ONTO, ONTO$id[i] );
        }

        ONTO$Nsubsumers[i] = length(as.vector(anc));
        
    }

    
    return(ONTO)
    
}


is.leaf.node <- function(ONTO,TARGET){
    children = ONTO$children[match(TARGET,ONTO$id)][[1]]
    return(ifelse(length(children)==0,TRUE,FALSE))
}

get.hp.leafs <- function(ONTO, SUB=FALSE){

    ONTO$Nleafs = rep(NA,length(ONTO$id))
    names(ONTO$Nleafs) = ONTO$id
    
    for( i in 1:length(ONTO$id) ){

        if( SUB ){
            desc = get.sub.decendants( ONTO, ONTO$id[i] );
        } else {
            desc = get.descendants( ONTO, ONTO$id[i] );
        }

        ONTO$Nleafs[i] = sum(as.vector(unlist(sapply(desc,is.leaf.node,ONTO=ONTO))))
        
    }

    return(ONTO)
    
}

get.onto.children <- function(ONTO,TARGET,LEVEL){
    children = ONTO$children[which(ONTO$id==TARGET)][[1]]
    indx     = match(children,ONTO$id)
    indx     = indx[!is.na(indx)]
    levels   = ONTO$levels[indx]
    levels   = ifelse( is.na(levels), LEVEL, levels)
    ONTO$levels[indx] = levels
    return(list(ONTO=ONTO,children=children))
}

get.next.level <- function(ONTO, ROOT, LEVEL ){

    res <- get.onto.children(ONTO,ROOT,LEVEL)
    ONTO     = res$ONTO
    CHILDREN = res$children
    
    while( length(CHILDREN) > 0 && sum(is.na(ONTO$levels)) > 0 ){     
        ONTO     = get.next.level(ONTO,CHILDREN[1],LEVEL+1)
        CHILDREN = CHILDREN[-1]
    }

    return(ONTO);
}

get.hp.levels <- function(ONTO){

    ONTO$levels = rep(NA,length(ONTO$id))
    names(ONTO$levels) = ONTO$id
    
    ROOT = root.terms(ONTO=ONTO)[1]

    LEVEL = 0
    
    ONTO$levels[which(ONTO$id==ROOT)] = LEVEL

    ONTO = get.next.level(ONTO,ROOT,LEVEL+1)

    return(ONTO)
    
}

## remove obsolete hp terms in onto
rm.obs.terms <- function(ONTO){

    Nlabels = length(names(ONTO))

    for( i in 1:Nlabels ){ ONTO[[i]] = ONTO[[i]][!ONTO$obsolete] }

    return(ONTO)

}

## return all root nodes
root.terms <- function(ONTO){

    Nparents = sapply(ONTO$parents, function(x) length(x) )
    roots    = names(Nparents)[Nparents == 0]
    obs      = names(onto$obsolete)[onto$obsolete]
    indx     = match(roots,obs)
    indx     = ifelse(is.na(indx),TRUE,FALSE)
    return(roots[indx])
    
}

## return all obsolete terms in ontology
obs.terms <- function(ONTO){ return(ONTO$obsolete) }

## remove annotation vector from ontology object
rm.onto.annotation <- function(ONTO, LABEL=NULL){

    if( !is.null(LABEL) ){
        indx = which(names(ONTO)==LABEL)
        if( length(indx) != 0 ){
            ONTO = ONTO[-indx]
        }
    }

    return(ONTO)
}

## add annotation vector from ontology object
add.onto.annotation <- function(ONTO, ANNO, HPLABEL="hpo_id", LABEL=NULL ){

    hpo.indx   = which(colnames(ANNO)==HPLABEL)
    label.indx = which(colnames(ANNO)==LABEL)

    if( length(hpo.indx) != 0 && length(label.indx) ){          
        indx        = match(ONTO$id,ANNO[,hpo.indx])
        ONTO$freq   = as.numeric(as.vector(ANNO[indx,label.indx]))      
        names(ONTO)[which(names(ONTO)=="freq")] = LABEL;
    }

    return(ONTO)
    
}

get.terms.in.onto <- function(ONTO,TARGET){
    indx = match(TARGET,ONTO$id,nomatch=NA)
    indx = indx[!is.na(indx)]
    return(as.vector(ONTO$id[indx]))
}

# get all ancestor terms to the given target term
get.ancestors <- function(ONTO,TARGET){
    ## example TRAGET = geneModel.anc.HP[1]
    return(get.sub.ancestors(SUBONTO=ONTO,TARGET=TARGET));
}


## get all decendant terms to the given target term
get.descendants <- function(ONTO,TARGET){
    ## example TRAGET = geneModel.anc.HP[1]
    return(get.sub.decendants(SUBONTO=ONTO,TARGET=TARGET));
}

get.term.frequency <- function( ONTO, TARGET, LABEL=NULL, PRINT=FALSE, SUB=FALSE ){
## Ref: https://www.cell.com/ajhg/fulltext/S0002-9297(09)00399-1
## The frequency of a term is defined as the proportion of objects that are annotated by the term or any of its descendent terms.

    freq   = 0;    
    OPTION = 0;
    Aindx  = 0;
    
    if( is.null(LABEL) ){
        cat('no annotation provided, calculate frequency based on term decendants.\n')
        OPTION = 0;
    } else {
        Aindx = which(names(ONTO)==LABEL)
        if( length(Aindx) == 0 ){
            if(PRINT){
                cat('no annotation provided, calculate frequency based on term decendants.\n')
            }
            OPTION = 0;
        } else {
            if(PRINT){
                cat('annotation ', LABEL, ' provided.\n')
            }
            OPTION = 1;
        }
    }

    Tindx = match(TARGET,ONTO$id,nomatch=NA)

    if( is.na(Tindx) ){
        if(PRINT){
            cat('term not found in ontology.')
        }
    } else {
    
        if( OPTION == 0 ){
            if( SUB ){
                desc = get.sub.decendants( ONTO, TARGET );
            } else {
                desc = get.descendants(ONTO, TARGET);
            }
            freq = 1 + length(desc);
        }


        if( OPTION == 1 ){
            freq = ifelse(is.na(ONTO[[Aindx]][Tindx]),0,as.numeric(ONTO[[Aindx]][Tindx]))
            if( SUB ){
                desc = get.sub.decendants( ONTO, TARGET );
            } else {
                desc  = get.descendants(ONTO, TARGET);
            }
            Dindx = match(desc,ONTO$id);
            Dfreq = ONTO[[Aindx]][Dindx];
            Dfreq = sum(Dfreq,na.rm=TRUE);            
            freq  = freq + ifelse(is.na(Dfreq),0,Dfreq);
        }

    }
        
    return(as.vector(freq));
            
}

get.model.IC <- function(ONTO, HP.TERMS, ROOT, MODEL.LABEL, SUB=FALSE, PRINT=FALSE,
                         PARM=list(k=0.5)){

    ## set any parameters needed
    k = PARM[[which(names(PARM)=="k")]][1]
    
    ## get frequency for each HP term in geneModel and geneModel.anc
    ROOT.IC = list()
    IC.anc  = list()

    cat("> get IC for each HP term, and ancestors terms, in gene model:", MODEL.LABEL, "...")
          
    ## (1) first get model frequency for root node, i.e. the normalisation factor
    ROOT.IC[[1]] = get.term.frequency( ONTO=ONTO, TARGET=ROOT,
                                      LABEL=MODEL.LABEL, PRINT=PRINT, SUB=SUB );
    names(ROOT.IC)[1] = MODEL.LABEL   

    ## (2) next get model frequency for all HP, and HP ancestor, terms
    tmp = sapply(HP.TERMS,get.term.frequency,
                 ONTO=ONTO,LABEL=MODEL.LABEL,PRINT=PRINT,SUB=SUB);

    ## (3) store HP term freqeuncies in IC.anc
    xx           = matrix(NA,ncol=3,nrow=length(HP.TERMS))
    colnames(xx) = c("HP.term","levels","Freq")
    xx[,1]       = HP.TERMS
    if( !is.na(match("levels",names(onto),nomatch=NA)) ){
        xx[,2]   = as.numeric(ONTO$levels[match(xx[,1],ONTO$id)])
    }
    xx[match(xx[,1],names(tmp)),3] = as.numeric(tmp)

    IC.anc[[1]]      = xx;
    names(IC.anc)[1] = MODEL.LABEL;

    cat("done.\n")

    cat("> 1) get probability for each HP term in each gene model...")

    ## (4) get probability (prob) for each gene model HP term and ancentor HP term
    indx  = which(colnames(IC.anc[[1]])=="prob.def")
    if( length(indx) != 0 ){ IC.anc[[1]] = IC.anc[[1]][,-indx]; }
    indx  = which(colnames(IC.anc[[1]])=="Freq")
    freq  = as.numeric(IC.anc[[1]][,indx])
    prob  = (freq+1)/(as.numeric(ROOT.IC[[1]])+1)
    IC.anc[[1]] = cbind(IC.anc[[1]],prob)
    colnames(IC.anc[[1]])[dim(IC.anc[[1]])[2]] = "prob.def"

    cat("done.\n")

    cat("> 2) get information content for each HP term in each gene model...")

    ## (5) get information content (IC) for each gene model HP term and ancentor HP term
    indx  = which(colnames(IC.anc[[1]])=="IC.def")
    if( length(indx) != 0 ){ IC.anc[[1]] = IC.anc[[1]][,-indx]; }
    indx  = which(colnames(IC.anc[[1]])=="prob.def")
    prob  = as.numeric(IC.anc[[1]][,indx])
    ic    = -log(prob)
    IC.anc[[1]] = cbind(IC.anc[[1]],ic)
    colnames(IC.anc[[1]])[dim(IC.anc[[1]])[2]] = "IC.def"

    cat("done.\n")

    cat("> 3) get information content with relative depth for each HP term in each gene model...")
    ## REF: [1] https://link.springer.com/article/10.1007/s10844-017-0479-y
    ##      [2] https://www.sciencedirect.com/science/article/pii/S1532046415002877#b0065
    ##      [3] https://www.sciencedirect.com/science/article/pii/S0950705120306948#b8
    indx  = which(colnames(IC.anc[[1]])=="IC.deg")
    if( length(indx) != 0 ){ IC.anc[[1]] = IC.anc[[1]][,-indx]; }
    indx  = which(colnames(IC.anc[[1]])=="Freq")
    freq  = as.numeric(IC.anc[[1]][,indx])
    indx  = which(colnames(IC.anc[[1]])=="levels")
    level = as.numeric(IC.anc[[1]][,indx])
    ic.deg  = k*(1 - (log(freq+1)/log(as.numeric(ROOT.IC[[1]])+1)));
    ic.deg  = ic.deg + (1-k) * (log(level+1)/log(max(level,na.rm=T)+1));
    IC.anc[[1]] = cbind(IC.anc[[1]],ic.deg)
    colnames(IC.anc[[1]])[dim(IC.anc[[1]])[2]] = "IC.deg"

    cat("done.\n")
    
    #cat("> 3) get entropy for each HP term in each gene model...")

    ## (6) get entropy for each gene model HP term and ancentor HP term
    #indx  = which(colnames(IC.anc[[1]])=="Entropy")
    #if( length(indx) != 0 ){ IC.anc[[1]] = IC.anc[[1]][,-indx]; }
    #indx  = which(colnames(IC.anc[[1]])=="prob.def")
    #prob  = as.numeric(IC.anc[[1]][,indx])
    #ent   = sapply(prob,entropy)
    #IC.anc[[1]] = cbind(IC.anc[[1]],ent)
    #colnames(IC.anc[[1]])[dim(IC.anc[[1]])[2]] = "Entropy"

    #cat("done.\n")

    return(list(IC.anc=IC.anc[[1]],ROOT.IC=ROOT.IC[[1]]))
    
}

## return the Maximum Information Content Ancestor, given HP terms
## HP.A and HP.B
MICA <- function(ONTO, HP.A, HP.B, IC, ANNOindx, HPindx, ICindx){

    mica = NA
    anc  = get_ancestors(ONTO,c(HP.A,HP.B))
    indx = match(anc,IC[[ANNOindx]][,HPindx])
    
    if( length(indx) != 0 ){
        set        = as.numeric(IC[[ANNOindx]][indx,ICindx])
        set        = set[!is.na(set)]
        mica       = as.numeric(set[which.max(set)])
    }

    return(mica)
    
}

## return index of the Maximum Information Content Ancestor,
## given HP terms HP.A and HP.B
MICA.INDX <- function(ONTO, HP.A, HP.B, IC, ANNOindx, HPindx, ICindx){

    mica = NA
    anc  = get_ancestors(ONTO,c(HP.A,HP.B))
    indx = match(anc,IC[[ANNOindx]][,HPindx])
    
    if( length(indx) != 0 ){
        set        = as.numeric(IC[[ANNOindx]][indx,ICindx])
        names(set) = indx
        set        = set[!is.na(set)]
        mica       = as.numeric(names(set)[which.max(set)])
    }

    return(mica)
    
}


precal.sim <- function(ONTO, HP.SET, IC, ANNO.SET=1, PRINT=TRUE,
                       TERM.LAB="HP.term",IC.LAB="IC.def"){

    Sim     = NULL       
    N       = length(HP.SET)
    Nfrac   = ceiling((20 * N)/100)
    
    HPindx  = which(colnames(IC[[ANNO.SET]])==TERM.LAB)#"HP.term") 
    ICindx  = which(colnames(IC[[ANNO.SET]])==IC.LAB)#"IC")    
    
    if( length(N) > 0 ){

        cat("> precal ", names(IC)[ANNO.SET], "...\n")
        
        Sim = matrix(NA,nrow=N,ncol=N)
        rownames(Sim) = HP.SET;
        colnames(Sim) = HP.SET;            

        ii = seq_along(Sim[1,])
        for( i in 1:N ){
             val = as.vector(sapply(rownames(Sim)[ii],
                                    MICA.INDX,
                                    ONTO=ONTO,
                                    HP.B=colnames(Sim)[i],
                                    ANNOindx=ANNO.SET,
                                    IC=IC,HPindx=HPindx,ICindx=ICindx))

             Sim[ii,i] = val
             Sim[i,ii] = val
             ii        = ii[-1]

             if( i %% Nfrac == 0 ){
                 if( PRINT ){ cat("> col: ", i, " (", ceiling((i/N)*100), "%)\n"); }
             }
        }
        
        cat("> done.\n")
        
    }#if
       
    return(list(Sim=Sim))    
    
}

## compute the symmetric similarity (Sim_max and Sim_avg)
## between two HP term sets, MODEL.A and MODEL.B, using MICA
pheno.sim <- function(ONTO, MODEL.A, MODEL.B, IC, ANNO.SET.A=1, ANNO.SET.B=1,
                      TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def"){
    
    ##REF: [1] https://www.cell.com/ajhg/fulltext/S0002-9297(19)30147-8#secsectitle0020 (PUBMED:31104773)
    ##     [2] https://www.cell.com/ajhg/fulltext/S0002-9297(09)00399-1 (PUBMED:19800049)
    ##     [3] https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-302#Sec9 (PUBMED:16776819)
    
    Sim.A   = NULL
    Sim.B   = NULL
    sim.max = NULL
    sim.avg = NULL
    jsd     = NULL
    
    Na = length(MODEL.A)
    Nb = length(MODEL.B)

    HPindx = which(colnames(IC[[ANNO.SET.A]])==TERM.LAB)#"HP.term")
    ICindx = which(colnames(IC[[ANNO.SET.A]])==IC.LAB)#"IC")
    PRindx = which(colnames(IC[[ANNO.SET.A]])==PR.LAB)#"prob")
    
    if( length(Na) > 0 && length(Nb) > 0 ){

        if( ANNO.SET.A == ANNO.SET.B ){
        
            Sim.A = matrix(NA,nrow=Na,ncol=Nb)
            rownames(Sim.A) = MODEL.A;
            colnames(Sim.A) = MODEL.B;

            for( i in 1:Nb ){
                Sim.A[,i] = as.vector(sapply(rownames(Sim.A),
                                             MICA,
                                             ONTO=ONTO,
                                             HP.B=colnames(Sim.A)[i],
                                             ANNOindx=ANNO.SET.A,
                                             IC=IC,HPindx=HPindx,ICindx=ICindx))
            }
             
            row.sum = sum(apply(Sim.A,1,max,na.rm=T))
            col.sum = sum(apply(Sim.A,2,max,na.rm=T))

            ## patient similarity based on MICA
            sim.max = 0.5 * ( row.sum + col.sum )
            sim.avg = 0.5 * ( (1/Na) * row.sum + (1/Nb) * col.sum )
    
            return(list(Sim=Sim.A,sim.max=sim.max,sim.avg=sim.avg,jsd=jsd))

        } else {

            tmp.A = matrix(NA,nrow=Na,ncol=Nb)            
            rownames(tmp.A) = MODEL.A;
            colnames(tmp.A) = MODEL.B;

            tmp.B = matrix(NA,nrow=Na,ncol=Nb)
            rownames(tmp.B) = MODEL.A;
            colnames(tmp.B) = MODEL.B;
            
            for( i in 1:Nb ){
                
                tmp.A[,i] = as.vector(sapply(rownames(tmp.A),
                                             MICA.INDX,
                                             ONTO=ONTO,
                                             HP.B=colnames(tmp.A)[i],
                                             ANNOindx=ANNO.SET.A,
                                             IC=IC,HPindx=HPindx,ICindx=ICindx))

                tmp.B[,i] = as.vector(sapply(rownames(tmp.B),
                                             MICA.INDX,
                                             ONTO=ONTO,
                                             HP.B=colnames(tmp.B)[i],
                                             ANNOindx=ANNO.SET.B,
                                             IC=IC,HPindx=HPindx,ICindx=ICindx))
                
            }

            Sim.A = apply(tmp.A,2,get.precal.data,IC=IC.anc[[ANNO.SET.A]],CINDX=ICindx)
            Sim.B = apply(tmp.B,2,get.precal.data,IC=IC.anc[[ANNO.SET.B]],CINDX=ICindx)
            
            Avg   = 0.5 * (Sim.A + Sim.B);
            
            row.sum = sum(apply(Avg,1,max,na.rm=T))
            col.sum = sum(apply(Avg,2,max,na.rm=T))

            ## patient similarity based on MICA
            sim.max = 0.5 * ( row.sum + col.sum )
            sim.avg = 0.5 * ( (1/Na) * row.sum + (1/Nb) * col.sum )

            ## jenson-shannon distance
            Sim.A = apply(tmp.A,2,get.precal.data,IC=IC.anc[[ANNO.SET.A]],CINDX=PRindx)
            Sim.B = apply(tmp.B,2,get.precal.data,IC=IC.anc[[ANNO.SET.B]],CINDX=PRindx)
            js    = sum(js_divergence(p=as.vector(Sim.A),q=as.vector(Sim.B)),na.rm=T)
            jsd   = sqrt(js)

            return(list(Sim=Avg,sim.max=sim.max,sim.avg=sim.avg,jsd=jsd))
            
        }#else


    }##length

    return(list(Sim=Sim.A,sim.max=sim.max,sim.avg=sim.avg,jsd=jsd))    
}

get.precal.data <- function( IC, RINDX, CINDX ){

    return(as.numeric(IC[RINDX,CINDX]))
    
}

## compute the symmetric similarity (Sim_max and Sim_avg)
## between two HP term sets, MODEL.A and MODEL.B, using precalulated MICA
pheno.sim.precal <- function(MODEL.A, MODEL.B, PRECAL, IC, ANNO.SET.A=1, ANNO.SET.B=1,
                             TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def"){

    Sim.A   = NULL
    Sim.B   = NULL
    sim.max = NULL
    sim.avg = NULL
    jsd     = NULL
    
    Na = length(MODEL.A)
    Nb = length(MODEL.B)

    HPindx = which(colnames(IC[[ANNO.SET.A]])==TERM.LAB)#"HP.term")
    ICindx = which(colnames(IC[[ANNO.SET.A]])==IC.LAB)#"IC")
    PRindx = which(colnames(IC[[ANNO.SET.A]])==PR.LAB)#"prob")
    
    if( length(Na) > 0 && length(Nb) > 0 ){

        if( ANNO.SET.A == ANNO.SET.B ){

            Rind  = match(MODEL.A,rownames(PRECAL[[ANNO.SET.A]]))
            Cind  = match(MODEL.B,colnames(PRECAL[[ANNO.SET.A]]))
            Sim.A = PRECAL[[ANNO.SET.A]][Rind,Cind]
            Sim.A = apply(Sim.A,2,get.precal.data,IC=IC.anc[[ANNO.SET.A]],CINDX=ICindx)
            
            row.sum = sum(apply(Sim.A,1,max,na.rm=T))
            col.sum = sum(apply(Sim.A,2,max,na.rm=T))

            ## patient similarity based on MICA
            sim.max = 0.5 * ( row.sum + col.sum )
            sim.avg = 0.5 * ( (1/Na) * row.sum + (1/Nb) * col.sum )       
            
            return(list(Sim=Sim.A,sim.max=sim.max,sim.avg=sim.avg,jsd=jsd))

        } else {

            Rind  = match(MODEL.A,rownames(PRECAL[[ANNO.SET.A]]))
            Cind  = match(MODEL.B,colnames(PRECAL[[ANNO.SET.A]]))
            tmp.A = PRECAL[[ANNO.SET.A]][Rind,Cind]
            Sim.A = apply(tmp.A,2,get.precal.data,IC=IC.anc[[ANNO.SET.A]],CINDX=ICindx)

            Rind  = match(MODEL.A,rownames(PRECAL[[ANNO.SET.B]]))
            Cind  = match(MODEL.B,colnames(PRECAL[[ANNO.SET.B]]))
            tmp.B = PRECAL[[ANNO.SET.B]][Rind,Cind]
            Sim.B = apply(tmp.B,2,get.precal.data,IC=IC.anc[[ANNO.SET.B]],CINDX=ICindx)

            Avg   = 0.5 * (Sim.A + Sim.B);
            
            row.sum = sum(apply(Avg,1,max,na.rm=T))
            col.sum = sum(apply(Avg,2,max,na.rm=T))

            ## patient similarity based on MICA
            sim.max = 0.5 * ( row.sum + col.sum )
            sim.avg = 0.5 * ( (1/Na) * row.sum + (1/Nb) * col.sum )

            ## jenson-shannon distance
            Sim.A = apply(tmp.A,2,get.precal.data,IC=IC.anc[[ANNO.SET.A]],CINDX=PRindx)
            Sim.B = apply(tmp.B,2,get.precal.data,IC=IC.anc[[ANNO.SET.B]],CINDX=PRindx)
            js    = sum(js_divergence(p=as.vector(Sim.A),q=as.vector(Sim.B)),na.rm=T)     
            jsd   = sqrt(js)

            return(list(Sim=Avg,sim.max=sim.max,sim.avg=sim.avg,jsd=jsd))

        }#else

    }#length

    return(list(Sim=Sim.A,sim.max=sim.max,sim.avg=sim.avg,jsd=jsd))
    
}


get.avg.prob <- function(MODEL.A, PRECAL, IC, ANNO.SET.A=1,
                         TERM.LAB="HP.term",IC.LAB="IC.def",PR.LAB="prob.def"){

    Sim.A    = NULL
    avg.prob = NULL
    
    Na = length(MODEL.A)

    HPindx = which(colnames(IC[[ANNO.SET.A]])==TERM.LAB)#"HP.term")
    ICindx = which(colnames(IC[[ANNO.SET.A]])==IC.LAB)#"IC")
    PRindx = which(colnames(IC[[ANNO.SET.A]])==PR.LAB)#"prob")
    
    if( length(Na) > 0 ){

        Rind  = match(MODEL.A,rownames(PRECAL[[ANNO.SET.A]]))
        Cind  = match(MODEL.A,colnames(PRECAL[[ANNO.SET.A]]))
        tmp.A = PRECAL[[ANNO.SET.A]][Rind,Cind]
        Sim.A = apply(tmp.A,2,get.precal.data,IC=IC.anc[[ANNO.SET.A]],CINDX=PRindx)

        avg.prob = sum(as.vector(Sim.A),na.rm=T)/length(Sim.A)
        
        return(list(avg.prob=avg.prob))

    }#length

    return(list(avg.prob=avg.prob))
}

    
##--------------------------------------------------------
## functions for dealing with a subontology
##--------------------------------------------------------

sub.ontology <- function(ONTO, TARGETS){

    SUBONTO = list()

    if( length(TARGETS) == 1 ){    
        desc = get.descendants(ONTO, TARGETS)
    } else {
        desc = TARGETS
    }
        
    Nlabels = length(names(ONTO))
    indx    = match(desc,names(ONTO$id))
    
    for( i in 1:Nlabels ){
        SUBONTO[[i]]      = ONTO[[i]][indx]
        names(SUBONTO)[i] = names(ONTO)[i]
    }

    return(SUBONTO)
    
}

get.sub.ancestors <- function( SUBONTO, TARGET ){
    ace  = get_ancestors(SUBONTO,TARGET)
    indx = match(ace,SUBONTO$id,nomatch=NA)
    indx = indx[!is.na(indx)]
    return(as.vector(SUBONTO$id[indx]))

}

get.sub.decendants <- function( SUBONTO, TARGET ){
    dec  = get_descendants(SUBONTO,TARGET)
    indx = match(dec,SUBONTO$id,nomatch=NA)
    indx = indx[!is.na(indx)]
    return(as.vector(SUBONTO$id[indx]))

}


##--------------------------------------------------------
## functions for dealing with hpo graph
##--------------------------------------------------------

##arrange df vars by position
## 'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
    ##REF: https://stackoverflow.com/questions/5620885/how-does-one-reorder-columns-in-a-data-frame
  
    ##stop if not a data.frame (but should work for matrices as well)
    stopifnot(is.data.frame(data))

    ##sort out inputs
    data.nms <- names(data)
    var.nr <- length(data.nms)
    var.nms <- names(vars)
    var.pos <- vars
    ##sanity checks
    stopifnot( !any(duplicated(var.nms)), 
               !any(duplicated(var.pos)) )
    stopifnot( is.character(var.nms), 
               is.numeric(var.pos) )
    stopifnot( all(var.nms %in% data.nms) )
    stopifnot( all(var.pos > 0), 
               all(var.pos <= var.nr) )

    ##prepare output
    out.vec <- character(var.nr)
    out.vec[var.pos] <- var.nms
    out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
    stopifnot( length(out.vec)==var.nr )

    ##re-arrange vars by position
    data <- data[ , out.vec]
    return(data)
}


removeVertexTerm <- function(GG,NAME){

    if( !is.null(igraph::get.vertex.attribute(GG,NAME)) ){
        GG <- igraph::remove.vertex.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.vertex.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.vertex.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}


build.hpo.graph <- function(ONTO=onto){

    ## mapping HP terms to integers
    N   = length(ONTO$id)    
    ids = cbind(as.vector(ONTO$id), seq(1,N,1))

    ## read-in onto edges and map HP terms to integers
    edges_df = data.frame(from=as.integer(),to=as.integer())
    for( i in 1:N ){
        edges    = ONTO$children[i][[1]]
        E        = length(edges)
        from     = rep(names(ONTO$children)[i],E)
        from     = as.integer(ids[match(from,ids[,1]),2])
        to       = edges
        to       = as.integer(ids[match(to,ids[,1]),2])
        tmp      = data.frame(from=from,to=to)
        edges_df = rbind(edges_df,tmp)
    }           

    ## build node attribute data frame for df
    rm.lab = c("parents","children","ancestors")
    nodes_df = ONTO[-match(rm.lab,names(ONTO))]
    names(nodes_df)[which(names(nodes_df)=="id")]   = "term"
    names(nodes_df)[which(names(nodes_df)=="name")] = "label"
    nodes_df$name                                   = as.integer(ids[,2])
    nodes_df$id                                     = as.integer(ids[,2])
    nodes_df = as.data.frame(nodes_df)   
    rownames(nodes_df) = NULL
    nodes_df = arrange.vars(nodes_df, c("name"=1))
    ##----
           
    cat("> build hpo network in igraph...")
        
    ## building graph (using igraph)
    gg = igraph::graph_from_data_frame(d=edges_df, vertices=nodes_df, directed=FALSE)
    gg = igraph::simplify(gg,remove.multiple=TRUE,remove.loops=TRUE)
    
    cat("done.\n")

    return(gg);

}

get.clip.hpo.graph <- function(GG, TARGETS=NULL){

    sub = NULL
    
    if( is.null(TARGETS) ){##try the clip
        TARGETS = V(GG)$term[V(GG)$clip=="1"]
    }

    ed   = get.edgelist(GG,names=T)
    vids = match(TARGETS,V(GG)$term,nomatch=NA)
    vids = vids[!is.na(vids)]

    ## can't find any targets
    if( length(vids) == 0 ){ return(sub); }
    
    sub  = igraph::induced.subgraph(GG,vids=vids)
    
    return(sub);
    
}
