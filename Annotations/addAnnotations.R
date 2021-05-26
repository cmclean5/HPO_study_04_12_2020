##clean
rm(list=ls())

source('../setUp.R')

removeVertexTerm <- function(GG,NAME){

    if( !is.null(igraph::get.vertex.attribute(GG,NAME)) ){
        GG <- igraph::remove.vertex.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.vertex.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.vertex.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}

getHPterms <- function(X){

    temp=strsplit(as.vector(X),split=",")
    temp=lapply(temp,function(x) gsub("\\[","",x))
    temp=lapply(temp,function(x) gsub("\\]","",x))
    temp=lapply(temp,function(x) gsub("\\'","",x))
    temp=lapply(temp,function(x) gsub(" ","",x))
    temp=unlist(temp)
    
    return(temp)
}

matches <- function(X,IDS){

    IDindx = match(X,IDS[,2])
    IDindx = IDindx[!is.na(IDindx)]

    return(length(IDindx))
}

mergeAnnotation <- function(X,Y){

    X = as.character(X)
    Y = as.character(Y)

    temp = unique(c(X,Y))

    if( length(temp) == 1 ){ return(temp) }    
    else {
        temp = temp[temp!=""]
        if( length(temp) == 1 ){ return(temp) }
        else { return(paste(temp,collapse=COLLAPSE[c])) }    
    }

}

addAnnotation <- function(SET,IDS,X,VATTS="",VALS=""){

    IDindx = match(X,IDS[,2])
    IDindx = IDindx[!is.na(IDindx)]
    
    for( v in 1:length(VATTS) ){

        newStr = VALS[v]

        if( newStr != "" ){
        
            Sindx = which(names(SET)==VATTS[v])

            if( length(Sindx) != 0 ){
            
                for( i in 1:length(IDindx) ){

                    oldStr = SET[[Sindx]][IDindx[i]]                   
                    
                    if( oldStr == "" ) {
                        SET[[Sindx]][IDindx[i]] <- as.character(newStr)
                    } else if( grepl(COLLAPSE[c], oldStr) ){
                        oldStr = strsplit(oldStr,COLLAPSE[c])[[1]]
                        oldStr = unique(c(newStr,oldStr))
                        SET[[Sindx]][IDindx[i]] <- paste(as.character(oldStr),collapse=COLLAPSE[c])
                    } else {
                        oldStr = unique(c(newStr,oldStr))
                        SET[[Sindx]][IDindx[i]]  <- paste(as.character(oldStr),collapse=COLLAPSE[c])
                    }
                }
            }          
        }
    }

    return(SET)
    
}

    
## read in annotation file
ann = read.delim("hpo_lists_per_disease.csv",sep=",",header=T)

grdir = "../Graphs"

## read in HPO network
gg  = read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
N   = length(V(gg))
ids = cbind(V(gg)$name,V(gg)$term)
    
annos    = vector(length=6)
annos[1] = "gene"
annos[2] = "disease.name"
annos[3] = "literature_model"
annos[4] = "ddd_model"
annos[5] = "OV"
annos[6] = "OVG"

annoSET  = list()
annoSET2 = list()

for( i in 1:6 ){
    val = rep("",length=N)  

    annoSET[[i]]       = val
    names(annoSET)[i]  = annos[i]

    annoSET2[[i]]      = val
    names(annoSET2)[i] = annos[i]
    
}


Nann=length(ann[,1])

LITset=rep(0,Nann)
DDDset=rep(0,Nann)

GNindx=which(colnames(ann)==annos[1])
DISindx=which(colnames(ann)==annos[2])
LITindx=which(colnames(ann)==annos[3])
DDDindx=which(colnames(ann)==annos[4])

all.lit = c()
all.ddd = c()

## for Similarity analysis
LIT.Gene.HPterms = list()
DDD.Gene.HPterms = list()

for( i in 1:Nann ){

    ID  = as.numeric(ann[i,1])
    X1  = getHPterms(ann[i,LITindx[1]])
    X2  = getHPterms(ann[i,DDDindx[1]])
    GN  = as.character(ann[i,GNindx[1]])
    DIS = as.character(ann[i,DISindx[1]])

    ## for Similarity analysis
    LIT.Gene.HPterms[[i]] = X1
    names(LIT.Gene.HPterms)[i] = sprintf("%s:%gx",GN,ID)
    DDD.Gene.HPterms[[i]] = X2
    names(DDD.Gene.HPterms)[i] = sprintf("%s:%gx",GN,ID)
    ##------
    
    all.lit = c(X1,all.lit)    
    all.ddd = c(X2,all.ddd)
    
    LITset[i] = as.numeric(matches(X=X1,IDS=ids))
    DDDset[i] = as.numeric(matches(X=X2,IDS=ids))

    annoSET = addAnnotation(SET=annoSET,IDS=ids,X=X1,VATTS=annos[c(1,2,3)],VALS=c(GN,DIS,sprintf("LIT:%gx",ID)))
    annoSET = addAnnotation(SET=annoSET,IDS=ids,X=X2,VATTS=annos[c(1,2,4)],VALS=c(GN,DIS,sprintf("DDD:%gx",ID)))
    
    annoSET2 = addAnnotation(SET=annoSET2,IDS=ids,X=X1,VATTS=annos[c(1,2,3)],VALS=c(GN,DIS,"LIT"))
    annoSET2 = addAnnotation(SET=annoSET2,IDS=ids,X=X2,VATTS=annos[c(1,2,4)],VALS=c(GN,DIS,"DDD"))
    
}

## for Similarity analysis
saveRDS(LIT.Gene.HPterms,"./LIT.Gene.HPterms.RDS",compress=T)
saveRDS(DDD.Gene.HPterms,"./DDD.Gene.HPterms.RDS",compress=T)

all.lit = unique(all.lit)
all.ddd = unique(all.ddd)
write.table(all.lit,"../DDD_dataset/all_lit_HPterms.csv",sep="\t",col.names=F,row.names=F,quote=F)
write.table(all.ddd,"../DDD_dataset/all_ddd_HPterms.csv",sep="\t",col.names=F,row.names=F,quote=F)

for( i in 1:N ){
    annoSET[[5]][i] =  mergeAnnotation(annoSET2[[3]][i],annoSET2[[4]][i])
    annoSET[[6]][i] =  mergeAnnotation(annoSET[[3]][i],annoSET[[4]][i])
}

for( i in 1:6 ){
    gg = removeVertexTerm(gg,annos[i])
    gg = set.vertex.attribute(gg,annos[i],V(gg),as.character(annoSET[[i]]))
}

##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), format="gml")


## output vertex attribute table
oo = matrix("",nrow=N,ncol=(4+length(annos)))
colnames(oo) = c("name","HPO.term","HPO.label","HPO.level",annos)

oo[,1] = V(gg)$name
oo[,2] = V(gg)$term
oo[,3] = V(gg)$label
oo[,4] = V(gg)$level
for(i in 1:length(annos)){
    oo[,(i+4)] = get.vertex.attribute(gg,annos[i],V(gg))
}

oo = oo[order(as.numeric(oo[,1])),]
write.table(oo,"vertex_attributues.csv",sep="\t",col.names=T,row.names=F,quote=F)
