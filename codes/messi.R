library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="configuration file name", metavar="character"),
  make_option(c("-t", "--train"), type="character", default=NULL, 
              help="configuration file name of the training", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
if (is.null(opt$train)){
  print_help(opt_parser)
  stop("One configuration or train file must be supplied (input file).n", call.=FALSE)
}}

source('messi_functions.R')

#Training the model

if (is.null(opt$train) == FALSE){
  source(opt$train)
  merclust(dat=dat,cum=TRUE,gal=TRUE,rank=FALSE,relaxed=FALSE,name.groups=name.groups,name.gal=name.gal,folder=folder,ntotal=0)

  dat_cum<-read.table(name.groups,header=T)
  id_mer<-1:length(dat_cum$ngroup)
  id_rel<-1:length(dat_cum$ngroup)
  for(i in 1:length(dat_cum$ngroup)){
    id_mer[i]=subset(cl_clas$id_mer,cl_clas$id == dat_cum$ngroup[i])[1] 
    id_rel[i]=subset(cl_clas$id_rel,cl_clas$id == dat_cum$ngroup[i])[1] 
  }
  dat_cum<-data.frame(dat_cum,id_mer,id_rel)
  write.table(dat_cum,file=name.groups,row.names=F,quote=F)

  dat_gal<-read.table(name.gal,header=T)
  dat_cum<-subset(dat_cum,dat_cum$id_mer==1)

  for(i in 1:length(dat_cum$ngroup)){
    sub<-subset(dat_gal,dat_gal$ngroup == dat_cum$ngroup[i])
    sub_aux<-subset(aux_id,aux_id$id == dat_cum$ngroup[i])
    if(i == 1){
      aux<-sub
      id<-sub_aux$class
    } else {
      aux<-rbind(aux,sub)
      id<-c(id,sub_aux$class)
    }
  }
  dat_gal<-data.frame(aux,id)
  write.table(dat_gal,file=name.gal,row.names=F,quote=F)
}


#Dynamical Classification
if (is.null(opt$file) == FALSE){
  source(opt$file)
  
  if(ncluster != 0){
    library('foreach')
    library('doParallel')
    library('doSNOW')

    merclust.par(dat=dat,cum=cum,gal=gal,rank=rank,relaxed=relaxed,name.groups=name.groups,name.gal=name.gal,rank.name=rank.name,relaxed.name=relaxed.name,folder=folder,nrank=nrank,name_trainset_cum=name_trainset_cum,name_trainset_gal=name_trainset_gal,ntotal=0,ncluster = ncluster)
  } else {
    merclust(dat=dat,cum=cum,gal=gal,rank=rank,relaxed=relaxed,name.groups=name.groups,name.gal=name.gal,rank.name=rank.name,relaxed.name=relaxed.name,folder=folder,nrank=nrank,name_trainset_cum=name_trainset_cum,name_trainset_gal=name_trainset_gal,ntotal=0)
  }

  beep(8)
}
