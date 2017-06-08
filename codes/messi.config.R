#---------------------------INPUT FILES-----------------------------------

#The name of the file where the galaxies information is stored, the default is 'prueba.dat'.
#It must have the same structure as 'prueba.dat'
input.file.name='galaxias_relax.dat'

#---------------------------TRAINSET FILES-----------------------------------

#The name of the trainset file for the clusters. The default is 'trainset_cum.dat'
#If you want to change the trainset for the clusters you must have
#  the same structure as 'trainset_cum.dat'
name_trainset_cum='trainset_cum.dat'
#The name of the trainset file for galaxies. The default is 'trainset_cum.dat'
#If you want to change the trainset for the galaxies you must have
#  the same structure as 'trainset_galaxies.dat'
name_trainset_gal='trainset_gal.dat'


#---------------------------OUTPUT FILES-----------------------------------

#The name of the folder where all the information will be stored 
#The default is 'prueba_messi'
folder='../sdss/analisis_sdss/messi_par'
#The name of the output file where the galaxies properties will be stored. 
#The default is 'estadisticos_galaxias.dat'
name.gal='estadisticos_galaxias.dat'
#The name of the output file where the clusters properties will be stored. 
#The default is 'estadisticos_grupos.dat'
name.groups='estadisticos_grupos.dat'
#The name of the output file where the merging-classification properties 
#  will be stored. 
#The default is 'ranking.dat'
rank.name='ranking.dat'
#The name of the output file where the relaxed-classification properties 
#  will be stored. 
#The default is 'ranking_relaxed.dat'
relaxed.name='ranking_relaxed.dat'
 
#-------------------------CONFIGURATION PARAMETERS--------------------------

nrank=30 #Number of random forest instances you want to do
cum=FALSE #change to FALSE if you have alredy estimate the clusters properties.
gal=FALSE #change to FALSE if you have alredy estimate the galaxy properties.
rank=TRUE #change to FALSE if you do not want to performed a merging classification
relaxed=FALSE #change to FALSE if you do not want to performed a relaxed classification
ncluster = 8 #Number of cluster to be used in the computation, if you set to 0, then the computation will not be parallel

dat<-read.table(input.file.name,header=T)

