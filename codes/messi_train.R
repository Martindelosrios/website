#---------------------------INPUT FILES-----------------------------------

#The name of the file where the galaxies information is stored.
#It must have the same structure as 'prueba_train.dat'
input.file.name='prueba_train.dat'
#The name of the file where the cluster classification is stored.
#It must have the same structure as 'id_cum.dat'
id.file.name='id_cum.dat'

#---------------------------OUTPUT FILES-----------------------------------

#The name of the output file where the clusters properties will be stored. 
#The default is 'trainset_cum.dat'
name.groups='trainset_cum.dat'
#The name of the output file where the galxies properties will be stored. 
#The default is 'trainset_gal.dat'
name.gal='trainset_gal.dat'
#The name of the folder where all the information will be stored 
#The default is the working directory
folder=getwd()

 

dat<-read.table(input.file.name,header=T)
aux_id<-subset(dat,select=c('id','class'))
dat<-subset(dat,select=c('ra','dec','z','mag','color','id'))
cl_clas<-read.table(id.file.name,header=T)
