options(scipen=999)
#Required libraries
library('nortest')
library('cosmoFns')
library('mclust')
library('e1071')
library('randomForest')
library('beepr')
#dressler(x,y,vel)  Calcula Delta2d/ngal para un cumulo.Coordenadas en grados decimales
#{{{
dressler<-function(x,y,vel){
	n=length(x)

        x=x*pi/180
        y=y*pi/180
#Calculo de la velocidad peculiar para ngal galaxias
	delta<-1:n
	for(j in 1:n){
                x0=x[j]
                y0=y[j]
		r<-distancia1(x,x0,y,y0)
		delta[j]<-ds1(vel,r)
			}
	DELTA=mean(delta)


	return(DELTA)

}
#}}}
#dressler2(x,y,vel)  Calcula Delta2d/ngal para un cumulo.Coordenadas en grados decimales
#{{{
dressler2<-function(x,y,vel){
	n=length(x)

        x=x*pi/180
        y=y*pi/180
#Calculo de la velocidad peculiar para ngal galaxias
	delta<-1:n
	for(j in 1:n){
                x0=x[j]
                y0=y[j]
		r<-distancia2(x,x0,y,y0)
		delta[j]<-ds1(vel,r)
			}
	DELTA=mean(delta)


	return(DELTA)

}
#}}}
#dressler.iter(x,y,vel)  Dressler shectman iterativo.Coordenadas en grados decimales
#{{{
dressler.iter<-function(ra,dec,vel){
  ngal=length(ra)
  ngal0=ngal/2
  ra1=ra*pi/180
  dec1=dec*pi/180
  vel1=vel
  mat<-data.frame(ra1,dec1,vel1)
  niter=30
  corte=1000
  c=0
#Calculo de la velocidad peculiar para ngal galaxias
  for(i in 1:niter){ 
    if(corte>0){
      if(ngal>ngal0){
          c=c+1
          n=length(mat$ra1)
          delta<-1:n
          for(j in 1:n){
                  x0=mat$ra1[j]
                  y0=mat$dec1[j]
                  ra1=mat$ra1
                  dec1=mat$dec1
                  vel1=mat$vel1
          	r<-distancia2(ra1,x0,dec1,y0)
          	delta[j]<-ds1(vel1,r)
          		}
          mat<-data.frame(ra1,dec1,vel1,delta)
          #mat<-subset(mat,mat$delta>1.4)
          mat<-subset(mat,mat$delta>0.5*mean(mat$delta))
          corte=ngal-length(mat$ra1)
          ngal=length(mat$ra1)
      } else {
          c=-1
      }  
    }
  }
  return(c)
}
#}}}
#dressler_graf(x,y,vel)  Calcula delta para c/galaxia de un cumulo.Coordenadas en grados decimales
#{{{
dressler_graf<-function(x,y,vel){
	n=length(x)
    
        x=x*pi/180
        y=y*pi/180

#Calculo de la velocidad peculiar para ngal galaxias
	delta<-1:n
	for(j in 1:n){
                x0=x[j]
                y0=y[j]
		r<-distancia1(x,x0,y,y0)
		delta[j]<-ds1(vel,r)
			}
	return(delta)

}
#}}}
#dressler_graf2(x,y,vel)  Calcula delta para c/galaxia de un cumulo.Coordenadas en grados decimales
#{{{
dressler_graf2<-function(x,y,vel){
	n=length(x)

        x=x*pi/180
        y=y*pi/180

#Calculo de la velocidad peculiar para ngal galaxias
	delta<-1:n
	for(j in 1:n){
                x0=x[j]
                y0=y[j]
		r<-distancia2(x,x0,y,y0)
		delta[j]<-ds1(vel,r)
			}
	return(delta)

}
#}}}
#distancia1(x,x0,y,y0) Calcula las distancias de cada galaxia a la galxia (x0,y0) y devuelve el vector con los indices
#{{{
distancia1<-function(x,x0,y,y0){
	n=length(x)
	r<-1:n
        nr=10
	for(i in 1:n){
		r[i]=sqrt((((x0-x[i])**2)*cos(y0)*cos(y0))+((y[i]-y0)**2))
			}
	r<-sort(r,index.return=TRUE)
        ri<-r$ix[2:(nr+1)]
	return(ri)
}
#}}}
#distancia2(x,x0,y,y0) Calcula las distancias de cada galaxia a la galxia (x0,y0) y devuelve el vector con los indices
#{{{
distancia2<-function(x,x0,y,y0){
	n=length(x)
	r<-1:n
        nr=floor(sqrt(n))
	for(i in 1:n){
		r[i]=sqrt((((x0-x[i])**2)*cos(y0)*cos(y0))+((y[i]-y0)**2))
			}
	r<-sort(r,index.return=TRUE)
        ri<-r$ix[2:(nr+1)]
	return(ri)
}
#}}}
#ds1(vel,r) calcula el delta de la galaxia cuyos indices de vecinos es r
#{{{
ds1<-function(vr,r){
	n=length(vr)
        vel_loc<-1:length(r)
#Calculo de velocidad media del cumulo
	velmed=mean(vr)
#calculo de la dispersion de velocidades del cumulo
	disp=0
	for(j in 1:n){
			disp=((velmed-vr[j])**2)+disp}
	disp=sqrt(disp/(n+1))
#Calculo de velocidad local media
	for(i in 1:length(r)){
             vel_loc[i]=vr[r[i]]
		}

	velloc=mean(vel_loc)
#Calculo de dispersion local de velocidades
	disp_loc=0
	for(i in 1:length(r)){
	        disp_loc=((vel_loc[i]-velloc)**2)+disp_loc		
		}
		disp_loc=sqrt(disp_loc/(length(r)+1))
#calculo del estimador delta

	delta=sqrt(((length(r)+1)/(disp**2))*(((velloc-velmed)**2)+((disp_loc-disp)**2)))
	
	return(delta)
}
#}}}
#mont(x,y,vel) calcula el pval2d para un cumulo.Coordenadas en grados decimales
#{{{
mont<-function(x,y,vel){
  n=length(vel)
  nran=0
  nmont=100
  del<-dressler(x,y,vel)

 for(i in 1:nmont){
  #print(i)
  vel<-sample(vel,size=n,replace=FALSE)
  del_mon<-dressler(x,y,vel)
  if(del_mon > del){
    nran=nran+1
    }
  }
 nran=nran/nmont
 return(nran)

}
#}}}
#mont_par(x,y,vel) calcula el pval2d para un cumulo de manera paralela.Coordenadas en grados decimales
#{{{
mont_par<-function(x,y,vel,ncluster){
  cl <- makeCluster(ncluster)
  registerDoParallel(cl)

  n=length(vel)
  nran=0
  nmont=30
  del<-dressler(x,y,vel)

  deltas_mont<-foreach(i=1:nmont,.combine='c',.export=c('dressler','distancia1','ds1')) %dopar% {
    vel<-sample(vel,size=n,replace=FALSE)
    del_mon<-dressler(x,y,vel)
    return(del_mon)
  }
 
 nran=length(which(deltas_mont > del))
 nran=nran/nmont

 stopCluster(cl)
 return(nran)

}
#}}}
#mont2(x,y,vel) calcula el pval2d para un cumulo.Coordenadas en grados decimales
#{{{
mont2<-function(x,y,vel){
  n=length(vel)
  nran=0
  nmont=100
  del<-dressler(x,y,vel)

 for(i in 1:nmont){
  #print(i)
  vel<-sample(vel,size=n,replace=FALSE)
  del_mon<-dressler2(x,y,vel)
  if(del_mon > del){
    nran=nran+1
    }
  }
 nran=nran/nmont
 return(nran)

}
#}}}
#gap(x)
#{{{
gap<-function(x){
#Esta funcion calcula el gap necesario para calcular la dispersion de velocidades de un cumulo.
#Datos de entrada: x= vector con las velocidades radiales de las galaxias.
#Datos de salida: gap. (Ver )	
	n=length(x)-1
	w<-1:2
	g<-1:2
	x<-sort(x,decreasing=FALSE)
	for(i in 1:n){
		w[i]=i*(n+1-i)
		g[i]=x[i+1]-x[i]
		}
	mul<-g*w
	mul<-sum(mul)
	return(mul)
}

#}}}
#dist(x,y,vel)
#{{{
dist<-function(x,y,vel){
#Esta funcion calcula las distancias proyectadas entre las galaxias. Necesaria para calcular el rvir
#Dato de entrada: (x,y,vel) vectores con las coordenadas de las galaxias en grados.
#Datos de salida: sum(1/r): Donde r es un vector que contiene las distancias proyectadas entre las galaxias en Mpc.
        x=x*pi/180
        y=y*pi/180
        id<-sort(vel,decreasing=FALSE,index.return=TRUE)$ix
        vel=vel[id]
        x=x[id]
        y=y[id]
	n=length(x)
	r<-1:2
	d=0
#pb <- txtProgressBar(title = "progress bar", min = 0,max = n, width = 82)
       for(i in 1:length(x)){
#setTxtProgressBar(pb, i, label=paste( round(i/n*100, 0),"% done"))
		if(i<length(x)){
		for(j in (i+1):length(x)){
                    if(x[i] != x[j]){
                    if(y[i] != y[j]){
                        d=d+1
			#r[d]=(sqrt((y[i]-y[j])**2+(cos(y[i])*(x[i]-x[j]))**2))*(vel[i]/100)/(1*(1+(vel[i]/300000)))
			r[d]=(sqrt((y[i]-y[j])**2+(cos(y[i])*(x[i]-x[j]))**2))
                        r[d]=r[d]*D.A(mean(vel/300000))
			}
                   }}
		}
		}
#close(pb)
	r=1/r
	r<-subset(r,r<10000000)
	r=sum(r)
	return(r)
	}
#}}}

#Esta funcion calcula el valor mas repetido de un vector.
#Dato de entrada: vector a estudiar.
#Dato de salida: Cantidad de veces que se repite el valor mas repetido.
#{{{
how1<-function(x){
	n=length(x)
	c=0
        l<-c(-1,-1)
	for(i in 1:n){
		if(length(x)>0){
		sub<-subset(x,x==x[1])
		c=c+1
                l[c]=length(sub)
		x<-subset(x,x!=x[1])
			}
		}
        l<-sort(l,decreasing=TRUE)
        c<-l[1]
	return(c)
}
#}}}

#mixt1(n,x,y,vel,delta) Realiza una mixtura de gaussianas pesadas por delta . n=numero de gaussianas permitidas
#{{{
mixt1<-function(n,x,y,vel,delta,delta_min){
  mat<-data.frame(x,y,vel,delta)
  mat<-subset(mat,mat$delta>delta_min)
  
  conta=0
  xnew<-1:length(mat$x)
  ynew<-1:length(mat$y)
  velnew<-1:length(mat$vel)
  for(k in 1:length(mat$x)){
  	dmin=min(mat$delta)
  	dmax=max(mat$delta)
  	npoints<-floor(((mat$delta[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-mat$x[k]
  		ynew[conta]<-mat$y[k]
  		velnew[conta]<-vel[k]
  		}
  	}

  mat<-data.frame(xnew,ynew) 
  Mclust(mat,G=1:n)->f
  fg<-f$classification
  mat<-data.frame(xnew,ynew,velnew,fg) 
  plot(x,y,pch=20,cex=2.2)
  #points(x0,y0,pch=20,col='red',cex=2.2)
  #points(x1,y1,pch=20,col='blue',cex=2.2)
  g1<-subset(mat,mat$fg==1)
  g2<-subset(mat,mat$fg==2)
  points(g1$xnew,g1$ynew,col='yellow',cex=0.5,pch=4)
  points(g2$xnew,g2$ynew,col='green',cex=0.5,pch=4)
   return(f)
}
#}}}
#mixt2(n,x,y,vel,delta) Realiza una mixtura de gaussianas. n=numero de gaussianas permitidas
#{{{
mixt2<-function(n,x,y,vel,delta,delta_min){
  mat<-data.frame(x,y,vel,delta)
  mat<-subset(mat,mat$delta>delta_min)
  mat<-data.frame(mat$x,mat$y) 
  Mclust(mat,G=1:n)->f
  fg<-f$classification
  mat1<-data.frame(mat$mat.x,mat$mat.y,fg) 
  plot(x,y,pch=20,cex=2.2)
  g1<-subset(mat1,mat1$fg==1)
  g2<-subset(mat1,mat1$fg==2)
  points(g1$mat.mat.x,g1$mat.mat.y,col='yellow',cex=0.5,pch=4)
  points(g2$mat.mat.x,g2$mat.mat.y,col='green',cex=0.5,pch=4)
  return(f)
}
#}}}

#mixt.mas(haloID,x,y,vel,delta,peso,delta_min, n_G=2) Realiza una mixtura de gaussianas pesadas por delta . n=numero de gaussianas permitidas. Sin grafico
#{{{
mixt.mas<-function(haloID,x,y,vel,delta,peso,delta_min,n_G=2){
  mat<-data.frame(x,y)
  Mclust(mat,G=1)->f
  sigma0_ra=f$parameters$variance$sigma[1,1,1]
  sigma0_dec=f$parameters$variance$sigma[2,2,1]
  ngal=length(x) 
  mat<-data.frame(x,y,vel,delta,peso,haloID)
  mat<-subset(mat,mat$delta>delta_min)

  id1=-99
  id2=-99
  l1=-99
  l2=-99
  l3=-99
  l4=-99
  rvir1=-99
  rvir2=-99
  rvir3=-99
  rvir4=-99
  dvel1=-99
  dvel2=-99
  dvel3=-99
  dvel4=-99
  mas1=-99
  mas2=-99
  mas3=-99
  mas4=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  racen1=-99
  deccen1=-99
  velcen1=-99
  racen2=-99
  deccen2=-99
  velcen2=-99
  racen3=-99
  deccen3=-99
  velcen3=-99
  racen4=-99
  deccen4=-99
  velcen4=-99
  nclus=-99
  nclus5=-99
  nmat=length(mat$delta)
  if(nmat>10){
  x1=mat$x
  y1=mat$y
  vel1=mat$vel
  haloid1=mat$haloID
#Mclust pesado por peso 
#{{{
  conta=0
  xnew<-1:2
  ynew<-1:2
  galid<-1:2
  for(k in 1:length(mat$x)){
  	dmin=min(mat$peso)
  	dmax=max(mat$peso)
  	npoints<-floor(((mat$peso[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-mat$x[k]
  		ynew[conta]<-mat$y[k]
                galid[conta]<-k
  		}
  	}

  mat<-data.frame(xnew,ynew) 
  Mclust(mat,G=1:n_G)->f
  ocu_rel=-999
  tot=-999
  if(f$G>1){
    nclus=f$G
    fg<-f$classification
    mat<-data.frame(xnew,ynew,fg,galid) 
    Gid<-1:length(x1)
    for(k in 1:length(x1)){
      s<-subset(mat,mat$galid==k)
      Gid[k]=s$fg[1]
      }
#}}} 
#Ocupacion relativa
#{{{
    mat<-data.frame(x1,y1,vel1,haloid1,Gid)
    name<-paste('galaxias_en_mock/',toString(backh),'_galaxies.dat',sep='')
    write.table(mat,file=name,row.names=FALSE)
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$Gid==j)
      l[c]=length(g$Gid)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    nclus5=length(which(l$x>5))
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    if(l1>4){
      if(l2>4){
        g<-subset(mat,mat$Gid==l$ix[1])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir1=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel1=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
        racen1=mean(g$x)*pi/180
        deccen1=mean(g$y)*pi/180
        velcen1=mean(g$vel)
        #mas1<-format(mas1,scientific=TRUE)
        

        g<-subset(mat,mat$Gid==l$ix[2])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir2=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel2=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
        racen2=mean(g$x)*pi/180 
        deccen2=mean(g$y)*pi/180
        velcen2=mean(g$vel)
        #mas2<-format(mas2,scientific=TRUE)

        if(length(subset(mat$Gid,mat$Gid==l$ix[3]))>4){
          g<-subset(mat,mat$Gid==l$ix[3])
          dd<-gap(g$vel)
          r<-dist(g$x,g$y,g$vel)
          rvir3=(pi*(length(g$x)-1)*length(g$x))/(2*r)
          dvel3=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
          mas3=(5*(dvel3**2)*rvir3)/(4.314465e-09)
          l3=l$x[3]
          racen3=mean(g$x)*pi/180 
          deccen3=mean(g$y)*pi/180
          velcen3=mean(g$vel)
        }
        if(length(subset(mat$Gid,mat$Gid==l$ix[4]))>4){
          g<-subset(mat,mat$Gid==l$ix[4])
          dd<-gap(g$vel)
          r<-dist(g$x,g$y,g$vel)
          rvir4=(pi*(length(g$x)-1)*length(g$x))/(2*r)
          dvel4=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
          mas4=(5*(dvel4**2)*rvir4)/(4.314465e-09)
          l4=l$x[4]
          racen4=mean(g$x)*pi/180 
          deccen4=mean(g$y)*pi/180
          velcen4=mean(g$vel)
        }

      }
    }
 
    #masas<-c(mas1,mas2,mas3,mas4)
    #radios_vir<-c(rvir1,rvir2,rvir3,rvir4)
    #disp_vel<-c(dvel1,dvel2,dvel3,dvel4)
    #ra_centros<-c(racen1,racen2,racen3,racen4)
    #dec_centros<-c(deccen1,deccen2,deccen3,deccen4)
    #vel_centros<-c(velcen1,velcen2,velcen3,velcen4)
    #masas<-sort(masas,index.return=TRUE,decreasing=T)
    #l1=l$x[masas$ix[1]]
    #l2=l$x[masas$ix[2]]

    masas<-c(mas1,mas2)
    radios_vir<-c(rvir1,rvir2)
    disp_vel<-c(dvel1,dvel2)
    ra_centros<-c(racen1,racen2)
    dec_centros<-c(deccen1,deccen2)
    vel_centros<-c(velcen1,velcen2)
    masas<-sort(masas,index.return=TRUE,decreasing=T)
    l1=l$x[masas$ix[1]]
    l2=l$x[masas$ix[2]]

    deccen=(dec_centros[masas$ix[1]]+dec_centros[masas$ix[2]])/2
    racen=(ra_centros[masas$ix[1]]+ra_centros[masas$ix[2]])/2
    velcen=(vel_centros[masas$ix[1]]+vel_centros[masas$ix[2]])/2
    dist12=(sqrt((dec_centros[masas$ix[1]]-dec_centros[masas$ix[2]])**2+(cos(deccen)*(ra_centros[masas$ix[1]]-ra_centros[masas$ix[2]]))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(radios_vir[masas$ix[1]]+radios_vir[masas$ix[2]])
    par2=abs(vel_centros[masas$ix[1]]-vel_centros[masas$ix[2]])/(disp_vel[masas$ix[1]]+disp_vel[masas$ix[2]])
#}}}
#Tamaño relativo
#{{{
sigma_ra<-1:f$G
sigma_dec<-1:f$G
for(i in 1:f$G){
  sigma_ra[i]=f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
#Asosiacion con subhalo
#{{{
#g<-subset(mat,mat$Gid==masas$ix[1])
g<-subset(mat,mat$Gid==l$ix[1])
id1=how(g$haloid1)
#g<-subset(mat,mat$Gid==masas$ix[2])
g<-subset(mat,mat$Gid==l$ix[2])
id2=how(g$haloid1)
#}}}
  }
  }
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2,nclus5,nclus)
  #v<-c(l1,l2,l3,mas3,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2,nclus5,nclus)
  return(v)
}
#}}}

#mixt(haloID,x,y,vel,delta,peso,delta_min,n_G=2) Realiza una mixtura de gaussianas pesadas por delta . n=numero de gaussianas permitidas. Sin grafico
#{{{
mixt<-function(haloID,x,y,vel,delta,peso,delta_min,n_G=2){
  mat<-data.frame(x,y)
  Mclust(mat,G=1)->f
  sigma0_ra=f$parameters$variance$sigma[1,1,1]
  sigma0_dec=f$parameters$variance$sigma[2,2,1]
  ngal=length(x) 
  mat<-data.frame(x,y,vel,delta,peso,haloID)
  mat<-subset(mat,mat$delta>delta_min)
  id1=-9
  id2=-99
  id3=-999
  id4=-9999
  l1=-99
  l2=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  racen1=-99
  deccen1=-99
  velcen1=-99
  racen2=-99
  deccen2=-99
  velcen2=-99
  nmat=length(mat$delta)
  nclus=-99
  nclus5=-99
  if(nmat>10){
  x1=mat$x
  y1=mat$y
  vel1=mat$vel
  haloid1=mat$haloID
  peso1=mat$peso
#Mclust pesado por peso 
#{{{
  conta=0
  xnew<-1:2
  ynew<-1:2
  galid<-1:2
  for(k in 1:length(mat$x)){
  	dmin=min(mat$peso)
  	dmax=max(mat$peso)
  	npoints<-floor(((mat$peso[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-mat$x[k]
  		ynew[conta]<-mat$y[k]
                galid[conta]<-k
  		}
  	}

  mat<-data.frame(xnew,ynew) 
  Mclust(mat,G=1:n_G)->f
  ocu_rel=-999
  tot=-999
  if(f$G>1){
    nclus=f$G
    fg<-f$classification
    mat<-data.frame(xnew,ynew,fg,galid) 
    Gid<-1:length(x1)
    for(k in 1:length(x1)){
      s<-subset(mat,mat$galid==k)
      Gid[k]=s$fg[1]
      }
#}}} 
#Ocupacion relativa
#{{{
    mat<-data.frame(x1,y1,vel1,haloid1,Gid,peso1)
    c=0
    l<-1:2
    pesos_medios<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$Gid==j)
      l[c]=length(g$Gid)
      pesos_medios[c]=mean(g$peso1)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    pesos_medios<-sort(pesos_medios,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    nclus5=length(which(l$x>5))
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    if(l1>4){
      if(l2>4){
        g<-subset(mat,mat$Gid==l$ix[1])
        #g<-subset(mat,mat$Gid==pesos_medios$ix[1])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir1=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel1=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
        racen1=mean(g$x)*pi/180
        deccen1=mean(g$y)*pi/180
        velcen1=mean(g$vel)
        #mas1<-format(mas1,scientific=TRUE)
        

        g<-subset(mat,mat$Gid==l$ix[2])
        #g<-subset(mat,mat$Gid==pesos_medios$ix[2])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir2=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel2=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
        racen2=mean(g$x)*pi/180 
        deccen2=mean(g$y)*pi/180
        velcen2=mean(g$vel)
        #mas2<-format(mas2,scientific=TRUE)
      }
    }

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12#/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)#/(dvel1+dvel2)
#}}}
#Tamaño relativo
#{{{
sigma_ra<-1:f$G
sigma_dec<-1:f$G
for(i in 1:f$G){
  sigma_ra[i]=f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
#Asosiacion con subhalo
#{{{
#g<-subset(mat,mat$Gid==pesos_medios$ix[1])
g<-subset(mat,mat$Gid==l$ix[1])
id1=how(g$haloid1)
#g<-subset(mat,mat$Gid==pesos_medios$ix[2])
g<-subset(mat,mat$Gid==l$ix[2])
id2=how(g$haloid1)
if(length(subset(mat$Gid,mat$Gid==l$ix[3])>0)){
  g<-subset(mat,mat$Gid==l$ix[3])
  id3=how(g$haloid1)
}
if(length(subset(mat$Gid,mat$Gid==l$ix[4])>0)){
  g<-subset(mat,mat$Gid==l$ix[4])
  id4=how(g$haloid1)
}
#}}}
  }
  }
  idis<-c(id1,id2,id3,id4)
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2,nclus5,nclus)
  salida<-list(v,idis)
  return(salida)
}
#}}}

#library('Rmixmod')
#mixt.mix3d(haloID,x,y,vel,delta,peso,delta_min,n_G=4) Realiza una mixtura de gaussianas pesadas por delta . n=numero de gaussianas permitidas. Sin grafico
#{{{
mixt.mix3d<-function(haloID,x,y,vel,delta,peso,delta_min,crit='BIC',n_G=4){
  sigma0_ra=-9  
  sigma0_dec=-9 
  ngal=length(x) 
  mat<-data.frame(x,y,vel,delta,peso,haloID)
  mat<-subset(mat,mat$delta>delta_min)
  id1=-9
  id2=-99
  id3=-999
  id4=-9999
  l1=-99
  l2=-99
  l3=-99
  l4=-99
  rvir1=-99
  rvir2=-99
  rvir3=-99
  rvir4=-99
  dvel1=-99
  dvel2=-99
  dvel3=-99
  dvel4=-99
  mas1=-99
  mas2=-99
  mas3=-99
  mas4=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  racen1=-99
  deccen1=-99
  velcen1=-99
  racen2=-99
  deccen2=-99
  velcen2=-99
  racen3=-99
  deccen3=-99
  velcen3=-99
  racen4=-99
  deccen4=-99
  velcen4=-99
  nclus=-99
  nclus5=-99
  nmat=length(mat$delta)
  if(nmat>10){
  peso1=mat$peso
#mixmod pesado por peso 
#{{{
  conta=0
  xnew<-1:2
  ynew<-1:2
  velnew<-1:2
  xnew0=mean(mat$x)
  ynew0=mean(mat$y)
  velnew0=mean(mat$vel)
  galid<-1:2
  for(k in 1:length(mat$x)){
  	dmin=min(mat$peso)
  	dmax=max(mat$peso)
  	npoints<-floor(((mat$peso[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-(mat$x[k]-xnew0)*cos(ynew0) 
  		ynew[conta]<-(mat$y[k]-ynew0)
  		velnew[conta]<-(mat$vel[k]-velnew0)
                galid[conta]<-k
  		}
  	}

  mat_aux<-data.frame(xnew,ynew,velnew) 
  mat_aux<-data.frame(sapply(mat_aux,FUN=scale))

  mixmodCluster(mat_aux,nbCluster=1:n_G,criterion=crit)->f
  ocu_rel=-999
  tot=-999
  if(f@bestResult@nbCluster>1){
    nclus=f@bestResult@nbCluster
    fg<-f@bestResult@partition
    mat_aux<-data.frame(mat_aux,fg,galid) 
    Gid<-1:length(peso)
    x1<-1:length(peso)
    y1<-1:length(peso)
    vel1<-1:length(peso)
    haloid1<-1:length(peso)
    peso1<-1:length(peso)
    for(k in 1:length(peso)){
      s<-subset(mat_aux,mat_aux$galid==k)
      Gid[k]=s$fg[1]
      x1[k]=mat$x[k]
      y1[k]=mat$y[k]
      vel1[k]=mat$vel[k]
      haloid1[k]=mat$haloID[k]
      peso1[k]=mat$peso[k]
    }
#}}} 
#Ocupacion relativa
#{{{
    mat<-data.frame(x1,y1,vel1,haloid1,Gid,peso1)
    name<-paste('galaxias_en_mock/',toString(backh),'_galaxies.dat',sep='')
    write.table(mat,file=name,row.names=FALSE)
    c=0
    l<-1:2
    pesos_medios<-1:2
    for(j in 1:f@bestResult@nbCluster){
      c=c+1
      g<-subset(mat,mat$Gid==j)
      l[c]=length(g$Gid)
      pesos_medios[c]=mean(g$peso1)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    pesos_medios<-sort(pesos_medios,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    nclus5=length(which(l$x>5))
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    if(l1>4){
      if(l2>4){
        g<-subset(mat,mat$Gid==l$ix[1])
        #g<-subset(mat,mat$Gid==pesos_medios$ix[1])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir1=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel1=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
        racen1=mean(g$x)*pi/180
        deccen1=mean(g$y)*pi/180
        velcen1=mean(g$vel)
        #mas1<-format(mas1,scientific=TRUE)
        

        g<-subset(mat,mat$Gid==l$ix[2])
        #g<-subset(mat,mat$Gid==pesos_medios$ix[2])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir2=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel2=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
        racen2=mean(g$x)*pi/180 
        deccen2=mean(g$y)*pi/180
        velcen2=mean(g$vel)
        #mas2<-format(mas2,scientific=TRUE)

        if(length(subset(mat$Gid,mat$Gid==l$ix[3]))>4){
          g<-subset(mat,mat$Gid==l$ix[3])
          dd<-gap(g$vel)
          r<-dist(g$x,g$y,g$vel)
          rvir3=(pi*(length(g$x)-1)*length(g$x))/(2*r)
          dvel3=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
          mas3=(5*(dvel3**2)*rvir3)/(4.314465e-09)
          l3=l$x[3]
          racen3=mean(g$x)*pi/180 
          deccen3=mean(g$y)*pi/180
          velcen3=mean(g$vel)
        }
        if(length(subset(mat$Gid,mat$Gid==l$ix[4]))>4){
          g<-subset(mat,mat$Gid==l$ix[4])
          dd<-gap(g$vel)
          r<-dist(g$x,g$y,g$vel)
          rvir4=(pi*(length(g$x)-1)*length(g$x))/(2*r)
          dvel4=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
          mas4=(5*(dvel4**2)*rvir4)/(4.314465e-09)
          l4=l$x[4]
          racen4=mean(g$x)*pi/180 
          deccen4=mean(g$y)*pi/180
          velcen4=mean(g$vel)
        }

      }
    }
 
    #masas<-c(mas1,mas2,mas3,mas4)
    #radios_vir<-c(rvir1,rvir2,rvir3,rvir4)
    #disp_vel<-c(dvel1,dvel2,dvel3,dvel4)
    #ra_centros<-c(racen1,racen2,racen3,racen4)
    #dec_centros<-c(deccen1,deccen2,deccen3,deccen4)
    #vel_centros<-c(velcen1,velcen2,velcen3,velcen4)
    #masas<-sort(masas,index.return=TRUE,decreasing=T)
    #l1=l$x[masas$ix[1]]
    #l2=l$x[masas$ix[2]]

    masas<-c(mas1,mas2)
    radios_vir<-c(rvir1,rvir2)
    disp_vel<-c(dvel1,dvel2)
    ra_centros<-c(racen1,racen2)
    dec_centros<-c(deccen1,deccen2)
    vel_centros<-c(velcen1,velcen2)
    masas<-sort(masas,index.return=TRUE,decreasing=T)
    l1=l$x[1]
    l2=l$x[2]

    deccen=(dec_centros[masas$ix[1]]+dec_centros[masas$ix[2]])/2
    velcen=(vel_centros[masas$ix[1]]+vel_centros[masas$ix[2]])/2
    dist12=(sqrt((dec_centros[masas$ix[1]]-dec_centros[masas$ix[2]])**2+(cos(deccen)*(ra_centros[masas$ix[1]]-ra_centros[masas$ix[2]]))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(radios_vir[masas$ix[1]]+radios_vir[masas$ix[2]])
    par2=abs(vel_centros[masas$ix[1]]-vel_centros[masas$ix[2]])/(disp_vel[masas$ix[1]]+disp_vel[masas$ix[2]])
#}}}
#Tamaño relativo
#{{{
sigma_ra<-1:f@bestResult@nbCluster
sigma_dec<-1:f@bestResult@nbCluster
for(i in 1:f@bestResult@nbCluster){
  sigma_ra[i]=99#f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=99#f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
#Asosiacion con subhalo
#{{{
#g<-subset(mat,mat$Gid==pesos_medios$ix[1])
g<-subset(mat,mat$Gid==l$ix[1])
id1=how(g$haloid1)
#g<-subset(mat,mat$Gid==pesos_medios$ix[2])
g<-subset(mat,mat$Gid==l$ix[2])
id2=how(g$haloid1)
if(length(subset(mat$Gid,mat$Gid==l$ix[3])>0)){
  g<-subset(mat,mat$Gid==l$ix[3])
  id3=how(g$haloid1)
} 
if(length(subset(mat$Gid,mat$Gid==l$ix[4])>0)){
  g<-subset(mat,mat$Gid==l$ix[4])
  id4=how(g$haloid1)
}
#}}}
  }
  }
  idis<-c(id1,id2,id3,id4)
  v<-c(l1,l2,l3,mas3,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2,nclus5,nclus)
  salida<-list(v,idis)
  return(salida)
}
#}}}
#mixt.mix(haloID,x,y,vel,delta,peso,delta_min) Realiza una mixtura de gaussianas pesadas por delta . n=numero de gaussianas permitidas. Sin grafico
#{{{
mixt.mix<-function(haloID,x,y,vel,delta,peso,delta_min,crit='BIC'){
  n_G=4
  mat<-data.frame(x,y)
  Mclust(mat,G=1)->f
  sigma0_ra=f$parameters$variance$sigma[1,1,1]
  sigma0_dec=f$parameters$variance$sigma[2,2,1]
  ngal=length(x) 
  mat<-data.frame(x,y,vel,delta,peso,haloID)
  mat<-subset(mat,mat$delta>delta_min)
  id1=-9
  id2=-99
  id3=-999
  id4=-9999
  l1=-99
  l2=-99
  l3=-99
  l4=-99
  rvir1=-99
  rvir2=-99
  rvir3=-99
  rvir4=-99
  dvel1=-99
  dvel2=-99
  dvel3=-99
  dvel4=-99
  mas1=-99
  mas2=-99
  mas3=-99
  mas4=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  racen1=-99
  deccen1=-99
  velcen1=-99
  racen2=-99
  deccen2=-99
  velcen2=-99
  racen3=-99
  deccen3=-99
  velcen3=-99
  racen4=-99
  deccen4=-99
  velcen4=-99
  nclus=-99
  nclus5=-99
  nmat=length(mat$delta)
  if(nmat>10){
  x1=mat$x
  y1=mat$y
  vel1=mat$vel
  peso1=mat$peso
  haloid1=mat$haloID
#mixmod pesado por peso 
#{{{
  conta=0
  xnew<-1:2
  ynew<-1:2
  galid<-1:2
  for(k in 1:length(mat$x)){
  	dmin=min(mat$peso)
  	dmax=max(mat$peso)
  	npoints<-floor(((mat$peso[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-mat$x[k]
  		ynew[conta]<-mat$y[k]
                galid[conta]<-k
  		}
  	}

  mat<-data.frame(xnew,ynew) 
  mixmodCluster(mat,nbCluster=1:n_G,criterion=crit)->f
  ocu_rel=-999
  tot=-999
  if(f@bestResult@nbCluster>1){
    nclus=f@bestResult@nbCluster
    fg<-f@bestResult@partition
    mat<-data.frame(xnew,ynew,fg,galid) 
    Gid<-1:length(x1)
    for(k in 1:length(x1)){
      s<-subset(mat,mat$galid==k)
      Gid[k]=s$fg[1]
      }
#}}} 
#Ocupacion relativa
#{{{
    mat<-data.frame(x1,y1,vel1,haloid1,Gid,peso1)
    name<-paste('galaxias_en_mock/',toString(backh),'_galaxies.dat',sep='')
    write.table(mat,file=name,row.names=FALSE)
    c=0
    l<-1:2
    pesos_medios<-1:2
    for(j in 1:f@bestResult@nbCluster){
      c=c+1
      g<-subset(mat,mat$Gid==j)
      l[c]=length(g$Gid)
      pesos_medios[c]=mean(g$peso1)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    pesos_medios<-sort(pesos_medios,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    nclus5=length(which(l$x>5))
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    if(l1>4){
      if(l2>4){
        g<-subset(mat,mat$Gid==l$ix[1])
        #g<-subset(mat,mat$Gid==pesos_medios$ix[1])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir1=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel1=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
        racen1=mean(g$x)*pi/180
        deccen1=mean(g$y)*pi/180
        velcen1=mean(g$vel)
        #mas1<-format(mas1,scientific=TRUE)
        

        g<-subset(mat,mat$Gid==l$ix[2])
        #g<-subset(mat,mat$Gid==pesos_medios$ix[2])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir2=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel2=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
        racen2=mean(g$x)*pi/180 
        deccen2=mean(g$y)*pi/180
        velcen2=mean(g$vel)
        #mas2<-format(mas2,scientific=TRUE)

        if(length(subset(mat$Gid,mat$Gid==l$ix[3]))>4){
          g<-subset(mat,mat$Gid==l$ix[3])
          dd<-gap(g$vel)
          r<-dist(g$x,g$y,g$vel)
          rvir3=(pi*(length(g$x)-1)*length(g$x))/(2*r)
          dvel3=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
          mas3=(5*(dvel3**2)*rvir3)/(4.314465e-09)
          l3=l$x[3]
          racen3=mean(g$x)*pi/180 
          deccen3=mean(g$y)*pi/180
          velcen3=mean(g$vel)
        }
        if(length(subset(mat$Gid,mat$Gid==l$ix[4]))>4){
          g<-subset(mat,mat$Gid==l$ix[4])
          dd<-gap(g$vel)
          r<-dist(g$x,g$y,g$vel)
          rvir4=(pi*(length(g$x)-1)*length(g$x))/(2*r)
          dvel4=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
          mas4=(5*(dvel4**2)*rvir4)/(4.314465e-09)
          l4=l$x[4]
          racen4=mean(g$x)*pi/180 
          deccen4=mean(g$y)*pi/180
          velcen4=mean(g$vel)
        }

      }
    }
 
    #masas<-c(mas1,mas2,mas3,mas4)
    #radios_vir<-c(rvir1,rvir2,rvir3,rvir4)
    #disp_vel<-c(dvel1,dvel2,dvel3,dvel4)
    #ra_centros<-c(racen1,racen2,racen3,racen4)
    #dec_centros<-c(deccen1,deccen2,deccen3,deccen4)
    #vel_centros<-c(velcen1,velcen2,velcen3,velcen4)
    #masas<-sort(masas,index.return=TRUE,decreasing=T)
    #l1=l$x[masas$ix[1]]
    #l2=l$x[masas$ix[2]]

    masas<-c(mas1,mas2)
    radios_vir<-c(rvir1,rvir2)
    disp_vel<-c(dvel1,dvel2)
    ra_centros<-c(racen1,racen2)
    dec_centros<-c(deccen1,deccen2)
    vel_centros<-c(velcen1,velcen2)
    masas<-sort(masas,index.return=TRUE,decreasing=T)
    l1=l$x[1]
    l2=l$x[2]

    deccen=(dec_centros[masas$ix[1]]+dec_centros[masas$ix[2]])/2
    velcen=(vel_centros[masas$ix[1]]+vel_centros[masas$ix[2]])/2
    dist12=(sqrt((dec_centros[masas$ix[1]]-dec_centros[masas$ix[2]])**2+(cos(deccen)*(ra_centros[masas$ix[1]]-ra_centros[masas$ix[2]]))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(radios_vir[masas$ix[1]]+radios_vir[masas$ix[2]])
    par2=abs(vel_centros[masas$ix[1]]-vel_centros[masas$ix[2]])/(disp_vel[masas$ix[1]]+disp_vel[masas$ix[2]])
#}}}
#Tamaño relativo
#{{{
sigma_ra<-1:f@bestResult@nbCluster
sigma_dec<-1:f@bestResult@nbCluster
for(i in 1:f@bestResult@nbCluster){
  sigma_ra[i]=99#f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=99#f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
#Asosiacion con subhalo
#{{{
#g<-subset(mat,mat$Gid==pesos_medios$ix[1])
g<-subset(mat,mat$Gid==l$ix[1])
id1=how(g$haloid1)
#g<-subset(mat,mat$Gid==pesos_medios$ix[2])
g<-subset(mat,mat$Gid==l$ix[2])
id2=how(g$haloid1)
if(length(subset(mat$Gid,mat$Gid==l$ix[3])>0)){
  g<-subset(mat,mat$Gid==l$ix[3])
  id3=how(g$haloid1)
} 
if(length(subset(mat$Gid,mat$Gid==l$ix[4])>0)){
  g<-subset(mat,mat$Gid==l$ix[4])
  id4=how(g$haloid1)
}
#}}}
  }
  }
  idis<-c(id1,id2,id3,id4)
  v<-c(l1,l2,l3,mas3,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2,nclus5,nclus)
  salida<-list(v,idis)
  return(salida)
}
#}}}
#mixt.real(x,y,vel,delta,peso,delta_min) Realiza una mixtura de gaussianas pesadas por delta . n=numero de gaussianas permitidas. Sin grafico
#{{{
mixt.real<-function(x,y,vel,mag,color,delta,peso,delta_min,grupo.id){
  n_G=2
  mat<-data.frame(x,y)
  Mclust(mat,G=1)->f
  sigma0_ra=f$parameters$variance$sigma[1,1,1]
  sigma0_dec=f$parameters$variance$sigma[2,2,1]
  ngal=length(x) 
  mat<-data.frame(x,y,vel,mag,color,delta,peso)
  #mat<-subset(mat,mat$delta>delta_min)
  id1=-99
  id2=-99
  l1=-99
  l2=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  racen1=-99
  deccen1=-99
  velcen1=-99
  racen2=-99
  deccen2=-99
  velcen2=-99
  nmat=length(mat$delta)
  if(nmat>10){
  x1=mat$x
  y1=mat$y
  vel1=mat$vel
  mag1=mat$mag
  color1=mat$color
#Mclust pesado por peso 
#{{{
  conta=0
  xnew<-1:2
  ynew<-1:2
  galid<-1:2
  for(k in 1:length(mat$x)){
  	dmin=min(mat$peso)
  	dmax=max(mat$peso)
  	npoints<-floor(((mat$peso[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-mat$x[k]
  		ynew[conta]<-mat$y[k]
                galid[conta]<-k
  		}
  	}

  mat<-data.frame(xnew,ynew) 
  Mclust(mat,G=1:n_G)->f
  ocu_rel=-999
  tot=-999
  if(f$G>1){
    fg<-f$classification
    mat<-data.frame(xnew,ynew,fg,galid) 
    Gid<-1:length(x1)
    for(k in 1:length(x1)){
      s<-subset(mat,mat$galid==k)
      Gid[k]=s$fg[1]
      }
#}}} 
#Ocupacion relativa
#{{{
    mat<-data.frame(x1,y1,vel1,mag1,color1,Gid)
    name<-paste(toString(grupo.id),'_galaxies.dat',sep='')

    if(file.exists(name)==FALSE){
     vec<-c('ra','dec','vel','mag','color','id')
     write.table(mat,file=name,row.names=FALSE,col.names=vec)
    }
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$Gid==j)
      l[c]=length(g$Gid)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    if(l1>4){
      if(l2>4){
        g<-subset(mat,mat$Gid==l$ix[1])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir1=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel1=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
        racen1=mean(g$x)*pi/180
        deccen1=mean(g$y)*pi/180
        velcen1=mean(g$vel)
        #mas1<-format(mas1,scientific=TRUE)
        

        g<-subset(mat,mat$Gid==l$ix[2])
        dd<-gap(g$vel)
        r<-dist(g$x,g$y,g$vel)
        rvir2=(pi*(length(g$x)-1)*length(g$x))/(2*r)
        dvel2=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
        mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
        racen2=mean(g$x)*pi/180 
        deccen2=mean(g$y)*pi/180
        velcen2=mean(g$vel)
        #mas2<-format(mas2,scientific=TRUE)
      }
    }

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)/(dvel1+dvel2)
#}}}
#Tamaño relativo
#{{{
sigma_ra<-1:f$G
sigma_dec<-1:f$G
for(i in 1:f$G){
  sigma_ra[i]=f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
  }
  }
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,racen1,racen2,deccen1,deccen2,velcen1,velcen2)
  return(v)
}
#}}}
#re.escaleo(x,y,vel)
#{{{
re.escaleo<-function(x,y,vel){
  
}
#}}}

#mixt.clas(haloID,x,y,vel,delta,delta_min) Realiza una mixtura de gaussianas sin pesar .Sin grafico
#{{{
mixt.clas<-function(haloID,x,y,vel,delta,delta_min){
  n_G=2
  mat<-data.frame(x,y)
  Mclust(mat,G=1)->f
  sigma0_ra=f$parameters$variance$sigma[1,1,1]
  sigma0_dec=f$parameters$variance$sigma[2,2,1]
  ngal=length(x) 
  mat<-data.frame(x,y,delta,vel,haloID)
  mat<-subset(mat,mat$delta>delta_min)
  id1=-99
  id2=-99
  l1=-99
  l2=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  nmat=length(mat$delta)
  if(nmat>10){
#Mclust 
#{{{
  mat1<-data.frame(mat$x,mat$y)
  Mclust(mat1,G=1:n_G)->f
  ocu_rel=-999
  tot=-999
  fg<-f$classification
  mat<-data.frame(mat1,mat$vel,fg,mat$haloID) 
#}}} 


if(f$G>1){
#Ocupacion relativa
#{{{
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$fg==j)
      l[c]=length(g$fg)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    g<-subset(mat,mat$fg==l$ix[1])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir1=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel1=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
    racen1=mean(g$mat.x)*pi/180
    deccen1=mean(g$mat.y)*pi/180
    velcen1=mean(g$mat.vel)
    #mas1<-format(mas1,scientific=TRUE)

    g<-subset(mat,mat$fg==l$ix[2])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir2=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel2=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
    racen2=mean(g$mat.x)*pi/180
    deccen2=mean(g$mat.y)*pi/180
    velcen2=mean(g$mat.vel)
    #mas2<-format(mas2,scientific=TRUE)

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)/(dvel1+dvel2)
  
#}}}
#Tamaño relativo
#{{{
sigma_ra<-1:f$G
sigma_dec<-1:f$G
for(i in 1:f$G){
  sigma_ra[i]=f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
#Asosiacion con subhalo
#{{{
g<-subset(mat,mat$fg==l$ix[1])
id1=how(g$mat.haloID)
g<-subset(mat,mat$fg==l$ix[2])
id2=how(g$mat.haloID)
#}}}
}
}
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2)
  return(v)
}
#}}}


#mixt.clas_graf(x,y,vel,delta,delta_min) Realiza una mixtura de gaussianas sin pesar .Con grafico
#{{{
mixt.clas_graf<-function(x,y,vel,delta,delta_min){
  ngal=length(x) 
  plot(x,y,pch=20,cex=0.5,xlab='Alfa',ylab='Delta')
  points(x,y,cex=exp(delta)/5)
  mat<-data.frame(x,y,delta,vel)
  mat<-subset(mat,mat$delta>delta_min)
  l1=-99
  l2=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  par1=-99
  par2=-99
  tot=-99
  nmat=length(mat$delta)
  if(nmat>10){
#Mclust 
#{{{
  mat1<-data.frame(mat$x,mat$y)
  Mclust(mat1)->f
  ocu_rel=-999
  tot=-999
  fg<-f$classification
  mat<-data.frame(mat1,mat$vel,fg) 
#}}} 


if(f$G>1){
#Ocupacion relativa
#{{{
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$fg==j)
      points(g$mat.x,g$mat.y,pch=20,col=(j+1))
      l[c]=length(g$fg)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    g<-subset(mat,mat$fg==l$ix[1])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir1=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel1=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
    racen1=mean(g$mat.x)*pi/180
    deccen1=mean(g$mat.y)*pi/180
    velcen1=mean(g$mat.vel)
    #mas1<-format(mas1,scientific=TRUE)

    g<-subset(mat,mat$fg==l$ix[2])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir2=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel2=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
    racen2=mean(g$mat.x)*pi/180
    deccen2=mean(g$mat.y)*pi/180
    velcen2=mean(g$mat.vel)
    #mas2<-format(mas2,scientific=TRUE)

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)/(dvel1+dvel2)
  
#}}}
}
}
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot)
  return(v)
}
#}}}

#mixt.clas_3d(haloID,x,y,vel,delta,delta_min) Realiza una mixtura de gaussianas sin pesar .Sin grafico
#{{{
mixt.clas_3d<-function(haloID,x,y,vel,delta,delta_min){
  n_G=2
  ngal=length(x) 
  mat<-data.frame(x,y)
  Mclust(mat,G=1)->f
  sigma0_ra=f$parameters$variance$sigma[1,1,1]
  sigma0_dec=f$parameters$variance$sigma[2,2,1]
  mat<-data.frame(x,y,delta,vel,haloID)
  mat<-subset(mat,mat$delta>delta_min)
  id1=-99
  id2=-99
  l1=-99
  l2=-99
  rvir1=-99
  rvir2=-99
  dvel1=-99
  dvel2=-99
  mas1=-99
  mas2=-99
  par1=-99
  par2=-99
  tot=-99
  sigma1_ra=-99
  sigma1_dec=-99 
  sigma2_ra=-99
  sigma2_dec=-99
  nmat=length(mat$delta)
  if(nmat>10){
#Mclust 
#{{{
  mat1<-data.frame(mat$x,mat$y,mat$vel)
  Mclust(mat1,G=1:n_G)->f
  ocu_rel=-999
  tot=-999
  fg<-f$classification
  mat<-data.frame(mat1,mat$vel,fg,mat$haloID) 
#}}} 


if(f$G>1){
#Ocupacion relativa
#{{{
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$fg==j)
      l[c]=length(g$fg)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    tot=(l1+l2)/ngal
#}}}
#Masa relativa, rvir, y dvel
#{{{
    g<-subset(mat,mat$fg==l$ix[1])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir1=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel1=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
    racen1=mean(g$mat.x)*pi/180
    deccen1=mean(g$mat.y)*pi/180
    velcen1=mean(g$mat.vel)
    #mas1<-format(mas1,scientific=TRUE)

    g<-subset(mat,mat$fg==l$ix[2])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir2=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel2=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
    racen2=mean(g$mat.x)*pi/180
    deccen2=mean(g$mat.y)*pi/180
    velcen2=mean(g$mat.vel)
    #mas2<-format(mas2,scientific=TRUE)

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)/(dvel1+dvel2)
  
#}}}

#Tamaño relativo
#{{{
sigma_ra<-1:f$G
sigma_dec<-1:f$G
for(i in 1:f$G){
  sigma_ra[i]=f$parameters$variance$sigma[1,1,i]
  sigma_dec[i]=f$parameters$variance$sigma[2,2,i]
 }
sigma1_ra=sigma_ra[l$ix[1]]/sigma0_ra
sigma1_dec=sigma_dec[l$ix[1]]/sigma0_dec
sigma2_ra=sigma_ra[l$ix[2]]/sigma0_ra
sigma2_dec=sigma_dec[l$ix[2]]/sigma0_dec
#}}}
#Asosiacion con subhalo
#{{{
g<-subset(mat,mat$fg==l$ix[1])
id1=how(g$mat.haloID)
g<-subset(mat,mat$fg==l$ix[2])
id2=how(g$mat.haloID)
#}}}
}
}
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot,sigma1_ra,sigma1_dec,sigma2_ra,sigma2_dec,id1,id2,racen1,racen2,deccen1,deccen2,velcen1,velcen2)
  return(v)
}
#}}}

#mixt_2(x,y,vel,delta,delta_min) Realiza una mixtura de gaussianas pesadas por delta .  Sin grafico
#{{{
mixt_2<-function(x,y,vel,delta,peso,delta_min){
  ngal=length(x) 
  mat<-data.frame(x,y,vel,delta,peso)
  mat<-subset(mat,mat$delta>delta_min)
  l1=-99
  l2=-99
rvir1=-99
rvir2=-99
dvel1=-99
dvel2=-99
mas1=-99
mas2=-99
  par1=-99
  par2=-99
  tot=-99
  nmat=length(mat$delta)
  if(nmat>10){
  x1=mat$x
  y1=mat$y
  vel1=mat$vel
#Mclust pesado por peso 
#{{{
  conta=0
  xnew<-1:2
  ynew<-1:2
  galid<-1:2
  for(k in 1:length(mat$x)){
  	dmin=min(mat$peso)
  	dmax=max(mat$peso)
  	npoints<-floor(((mat$peso[k]-dmin)/(dmax-dmin))*49+1)
  	for(kk in 1:npoints){
  		conta=conta+1
  		xnew[conta]<-mat$x[k]
  		ynew[conta]<-mat$y[k]
                galid[conta]<-k
  		}
  	}

  mat<-data.frame(xnew,ynew) 
  Mclust(mat,G=1:2)->f
  ocu_rel=-999
  tot=-999
  if(f$G>1){
    fg<-f$classification
    mat<-data.frame(xnew,ynew,fg,galid) 
    Gid<-1:length(x1)
    for(k in 1:length(x1)){
      s<-subset(mat,mat$galid==k)
      Gid[k]=s$fg[1]
      }
    }
#}}} 
#Ocupacion relativa
#{{{
    mat<-data.frame(x1,y1,vel1,Gid)
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$Gid==j)
      l[c]=length(g$Gid)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    tot=(l1+l2)/sum(l$x)
#}}}
#Masa relativa, rvir, y dvel
#{{{
    g<-subset(mat,mat$Gid==l$ix[1])
    dd<-gap(g$vel)
    r<-dist(g$x,g$y,g$vel)
    rvir1=(pi*(length(g$x)-1)*length(g$x))/(2*r)
    dvel1=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
    mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
    racen1=mean(g$x)*pi/180
    deccen1=mean(g$y)*pi/180
    velcen1=mean(g$vel)
    #mas1<-format(mas1,scientific=TRUE)

    g<-subset(mat,mat$Gid==l$ix[2])
    dd<-gap(g$vel)
    r<-dist(g$x,g$y,g$vel)
    rvir2=(pi*(length(g$x)-1)*length(g$x))/(2*r)
    dvel2=(sqrt(pi)*dd)/(length(g$x)*(length(g$x)-1))
    mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
    racen2=mean(g$x)*pi/180
    deccen2=mean(g$y)*pi/180
    velcen2=mean(g$vel)
    #mas2<-format(mas2,scientific=TRUE)

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)/(dvel1+dvel2)
  
#}}}
}
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot)
  return(v)
}
#}}}

#mixt.clas_2(x,y,vel,delta,delta_min) Realiza una mixtura de gaussianas sin pesar .Sin grafico
#{{{
mixt.clas_2<-function(x,y,vel,delta,delta_min){
  ngal=length(x) 
  mat<-data.frame(x,y,delta,vel)
  mat<-subset(mat,mat$delta>delta_min)
  l1=-99
  l2=-99
rvir1=-99
rvir2=-99
dvel1=-99
dvel2=-99
mas1=-99
mas2=-99
  tot=-99
  par1=-99
  par2=-99
  nmat=length(mat$delta)
  if(nmat>10){

#Mclust 
#{{{
  mat1<-data.frame(mat$x,mat$y)
  Mclust(mat1,G=1:2)->f
  ocu_rel=-999
  tot=-999
  fg<-f$classification
  mat<-data.frame(mat1,mat$vel,fg) 
#}}} 

if(f$G>1){
#Ocupacion relativa
#{{{
    c=0
    l<-1:2
    for(j in 1:f$G){
      c=c+1
      g<-subset(mat,mat$fg==j)
      l[c]=length(g$fg)
    }
    l<-sort(l,decreasing=TRUE,index.return=TRUE)
    l1=l$x[1]
    l2=l$x[2]
    tot=(l1+l2)/sum(l$x)
#}}}
#Masa relativa, rvir, y dvel
#{{{
    g<-subset(mat,mat$fg==l$ix[1])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir1=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel1=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas1=(5*(dvel1**2)*rvir1)/(4.314465e-09)
    racen1=mean(g$mat.x)*pi/180
    deccen1=mean(g$mat.y)*pi/180
    velcen1=mean(g$mat.vel)
    #mas1<-format(mas1,scientific=TRUE)

    g<-subset(mat,mat$fg==l$ix[2])
    dd<-gap(g$mat.vel)
    r<-dist(g$mat.x,g$mat.y,g$mat.vel)
    rvir2=(pi*(length(g$mat.x)-1)*length(g$mat.x))/(2*r)
    dvel2=(sqrt(pi)*dd)/(length(g$mat.x)*(length(g$mat.x)-1))
    mas2=(5*(dvel2**2)*rvir2)/(4.314465e-09)
    racen2=mean(g$mat.x)*pi/180
    deccen2=mean(g$mat.y)*pi/180
    velcen2=mean(g$mat.vel)
    #mas2<-format(mas2,scientific=TRUE)

    deccen=(deccen1+deccen2)/2
    velcen=(velcen1+velcen2)/2
    dist12=(sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000)))
    par1=dist12/(rvir1+rvir2)
    par2=abs(velcen1-velcen2)/(dvel1+dvel2)
  
#}}}
}
}
  v<-c(l1,l2,rvir1,rvir2,dvel1,dvel2,mas1,mas2,par1,par2,tot)
  return(v)
}
#}}}

#subs(halid,gid,masas,radios) Calcula la contaminacion, pureza y completitud de un subhalo.
#{{{
subs<-function(halid,gid,masas,radios){
  dat<-data.frame(halid,gid,masas,radios)
  dat<-subset(dat,dat$gid != 0)


  how(dat$gid)->n1

  g1<-subset(dat,dat$gid==n1)  
  ltot=length(g1$gid)
  ltot1=ltot #numero de galaxias del grupo identificado

  how(g1$halid)->sh1 #subhaloid del halo asosiado
  s1<-subset(halid,halid==sh1)
  l0=length(s1) #numero de galaxias del subhalo asosiado
  g1<-subset(g1,g1$halid==sh1)  
  masa1<-g1$masas[1]
  radio1<-g1$radios[1]
  l1=length(g1$halid) #numero de galaxias del grupo que pertenecen al subhalo asosiado

  cont1=(ltot1-l1)/l0 
  pur1=l0/ltot1
  compl1=l1/l0
  dat<-subset(dat,dat$gid != n1)

  how(dat$gid)->n2

  g2<-subset(dat,dat$gid==n2)  
  ltot=length(g2$gid)
  ltot2=ltot

  how(g2$halid)->sh2
  s2<-subset(halid,halid==sh2)
  l02=length(s2)
  g2<-subset(g2,g2$halid==sh2)  
  masa2<-g2$masas[1]
  radio2<-g2$radios[1]
  l2=length(g2$halid)

  cont2=(ltot2-l2)/l02 
  pur2=l0/ltot2
  compl2=l2/l02
 
  vec<-c(ltot1,cont1,pur1,compl1,masa1,radio1,ltot2,cont2,pur2,compl2,masa2,radio2)
  return(vec)

}
#}}}
#library("NORMT3") #Libreria que contiene a la funcion erf
#AD(x)
#Esta funcion calcula el estadistico A^2 de una distribucion
#Datos de entrada: x=vector de datos 
#Datos de salida: Alfa= nivel de significancia de la hipotesis nula. Alfa chico implica una distribucion no gaussiana.
#{{{
AD<-function(x){

   x<-sort(x,decreasing=FALSE)
   n=length(x)
   mean=mean(x)
   std.dev(x)->disp
   aux=0
   for(i in 1:n){
     x1=x[i]
     aux1<-CDF(x1,mean,disp)	
     x1=x[(n+1)-i]
     aux2<-CDF(x1,mean,disp)	
     aux=((2*i-1)*(log(aux1)+log(1-aux2)))+aux	
     }
   A2=-n-(1/n)*aux
   A2=A2*(1+4/n-25/n) 
   a=3.6789468
   b=0.1749916
   alfa=a*exp(-A2/b)
   return(alfa)
   #return(A2)

}
#}}}
#CDF(x,mean,std)
#Esta funcion calcula la funcion distribucion acumulada de una gaussiana con mu=mean y sigma=std. Necesaria para la funcion AD()
#Datos de salida: CDF(x)
#{{{
CDF<-function(x0,mean,std){
   x0=(x0-mean)/(std*sqrt(2))
   cdf=1+Re(erf(x0))
   cdf=cdf/2
   return(cdf)
}
#}}}
#std.dev(x)
#Datos de entrada: x=vector de datos
#Esta funcion calcula la desviacion estandar de un conjunto de puntos
#Datos de salida: std
#{{{
std.dev<-function(x){
   n=length(x)
   mean<-mean(x)
   std=0
   for(i in 1:n){
      x0=x[i]
      std=std+((x0-mean)**2)
      }
   std=std/(n+1)
   std=sqrt(std)
   return(std)

}
#}}}
#sph.test(x,y)
#Datos de entrada: x,y= coordenadas
#Datos de salida: gap.sph (gap.sph=1 esferico, gap.sph >> 1 no esferico)
#{{{
sph.test<-function(x,y){
   x0=mean(x)
   y0=mean(y)
   x=x-x0
   y=y-y0
   dat<-data.frame(x,y)

   dev.std<-1:10
   for(i in 1:10){
      tita=2*pi*(i-1)/10
      xnew=dat$x*cos(tita)-dat$y*sin(tita)
      ynew=dat$x*sin(tita)+dat$y*cos(tita)
      dev.std[i]=sd(xnew)
   }

   dev.std<-sort(dev.std,decreasing=T)
   gap.sph=dev.std[1]/dev.std[length(dev.std)]   
   return(gap.sph)
}
#}}}
#dens.test(x,y)
#{{{
dens.test<-function(x,y){
  x0=mean(x)
  y0=mean(y)
  x=x-x0
  y=y-y0
  dat<-data.frame(x,y)

  cuad1<-subset(dat,dat$x > 0 & dat$y > 0)
  cuad2<-subset(dat,dat$x > 0 & dat$y < 0)
  cuad3<-subset(dat,dat$x < 0 & dat$y < 0)
  cuad4<-subset(dat,dat$x < 0 & dat$y > 0)

  tita1=atan2(cuad1$y,cuad1$x) 
  tita2=(pi/2-atan2(cuad2$y,cuad2$x))
  tita3=(atan2(cuad3$y,cuad3$x)+2*pi)
  tita4=(atan2(cuad4$y,cuad4$x)+pi)

  cuad1<-data.frame(cuad1,tita1)
  colnames(cuad1)<-c('x','y','tita')
  cuad2<-data.frame(cuad2,tita2)
  colnames(cuad2)<-c('x','y','tita')
  cuad3<-data.frame(cuad3,tita3)
  colnames(cuad3)<-c('x','y','tita')
  cuad4<-data.frame(cuad4,tita4)
  colnames(cuad4)<-c('x','y','tita')

  dat<-rbind(cuad1,cuad2,cuad3,cuad4)

  dens<-1:2
  for(i in 1:10){
    tita.min=2*pi*(i-1)/10
    tita.max=2*pi*i/10
    sub<-subset(dat,dat$tita > tita.min & dat$tita < tita.max)
    dens[i]=length(sub$tita)
  }

  dens<-sort(dens,decreasing=T)
 
  if(dens[length(dens)] > 0){ 
    gap<-dens[1]/dens[length(dens)]
  } else {
    gap=99999
  }
  return(gap)
}
#}}}
#b.skew(data) bootstrap Skewness

 #{{{
b.skew <- function(data){
#bootstrap para la skewness, devuelbe el error standard
    num=length(data)/2
    if(num>30)num=30
    resamples <- lapply(1:num, function(i) sample(data, replace=T))
    r.skew <- sapply(resamples, skewness)
    std.err <- sqrt(var(r.skew))
    #list(std.err=std.err, resamples=resamples, skew=r.skew)
    return(std.err)
} #}}}*/
#how(x)
#Esta funcion calcula cual es el valor que mas se repite en un vector.
#Dato de entrada: Vector a estudiar.
#Dato de salida: Valor mas repetido.
#{{{
how<-function(x){
	n=length(x)
	r<-c(0,0)
	y<-c(0,0)
	c=0
	for(i in 1:n){
		if(length(x)>0){
		sub<-subset(x,x==x[1])
		c=c+1
		r[c]=length(sub)
		y[c]=x[1]
		x<-subset(x,x!=x[1])
			}
		}
	r<-sort(r,decreasing=TRUE,index.return=TRUE)
	rr<-y[r$ix[1]]
	return(rr)
}
#}}}
#how2(x)
#Esta funcion calcula cuales son los valores que mas se repiten en un vector.
#Dato de entrada: Vector a estudiar.
#Dato de salida: data frame con los Valores mas repetidos y la cantidad de veces que se repitieron.
#{{{
how2<-function(x){
	n=length(x)
	r<-c(0,0)
	y<-c(0,0)
	c=0
	for(i in 1:n){
		if(length(x)>0){
		sub<-subset(x,x==x[1])
		c=c+1
		r[c]=length(sub)
		y[c]=x[1]
		x<-subset(x,x!=x[1])
			}
		}
	r<-sort(r,decreasing=TRUE,index.return=TRUE)
        
        y1<-1:length(r$x)
        for(i in 1:length(r$x)){
          y1[i]=y[r$ix[i]]
        }
        rr<-data.frame(y1,r$x)
	return(rr)
}
#}}}
#how1(x)
#Esta funcion calcula cual es el valor que mas se repite en un vector.
#Dato de entrada: Vector a estudiar.
#Dato de salida: Numero de veces que se repite el valor mas repetidos.
#{{{
how1<-function(x){
	n=length(x)
	r<-c(0,0)
	y<-c(0,0)
	c=0
	for(i in 1:n){
		if(length(x)>0){
		sub<-subset(x,x==x[1])
		c=c+1
		r[c]=length(sub)
		y[c]=x[1]
		x<-subset(x,x!=x[1])
			}
		}
	r<-sort(r,decreasing=TRUE,index.return=TRUE)
	rr<-r$x[1]
	return(rr)
}
#}}}
#plot.elipse(x0,a,b,lty=1,col='black',lwd=1)
#{{{
plot.elipse<-function(x0,a,b,lty=1,col='black',lwd=1){
  theta <- seq(0, 2 * pi, length=(1000))
  x <- x0[1] + a * cos(theta)
  y <- x0[2] + b * sin(theta)
  points(x, y, type = "l",lty=lty,col=col,lwd=lwd) 
}
#}}}
#part(x,y,nstep) x:deltas de los positivos. y: deltas de los falsos.
#library('pracma')
#{{{
part<-function(x,y,nstep){
	v<-c(x,y)
	minval<-min(v)
	maxval<-max(v)
	dx=(maxval-minval)/(nstep-1)
	nx=length(x)
	ny=length(y)
		
	tpr<-1:nstep
	fpr<-1:nstep
	l<-1:nstep
	for(i in 1:nstep){
		lim=minval+(i-1)*dx	

		x<-subset(x,x>lim)
		y<-subset(y,y>lim)

		nx1=0
		ny1=0

		if(length(x) > 0){
		nx1=length(x)}
		if(length(y)>0){
		ny1=length(y)}

		tpr[i]=nx1/nx
		fpr[i]=ny1/ny
		l[i]=lim
	}	
	v<-data.frame(fpr,tpr,l)


        tprmax=max(tpr-fpr)

        delmax=which.max(tpr-fpr)
        fprmax=fpr[delmax]

        delmax=minval+(delmax-1)*dx

        area<-abs(trapz(fpr,tpr))

        leg0<-paste('max(TPR-FPR)=',format(tprmax,digits=2),sep='')
        leg1<-paste('FPR=',format(fprmax,digits=3),sep='')
        leg2<-paste('Delta/ngal=',format(delmax,digits=3),sep='')
        leg3<-paste('Area=',format(area,digits=3),sep='')
        leg<-c(leg0,leg1,leg2,leg3)
#GRAFICOS	
        par(lwd=3)
        par(cex=1.3)
	plot(v$fpr,v$tpr,cex=0.5,pch=19,typ="b",xlab="FPR",ylab="TPR")
	x<-c(0,1)
	lines(x,x)
        legend('bottomright',leg)
#/////////////////////////////////////////////////////////////////

	#return(v)
}
#}}}
#part2(x,y,nstep) x:pval de los positivos. y: pval de los falsos.
#{{{

part2<-function(x,y,nstep){
	v<-c(x,y)
	minval<-min(v)
	maxval<-max(v)
	dx=(maxval-minval)/(nstep-1)
	nx=length(x)
	ny=length(y)
		
	tpr<-1:nstep
	fpr<-1:nstep
	l<-1:nstep
	for(i in 1:nstep){
		lim=maxval-(i-1)*dx	

		x<-subset(x,x<=lim)
		y<-subset(y,y<=lim)

		nx1=0
		ny1=0

		if(length(x) > 0){
		nx1=length(x)}
		if(length(y)>0){
		ny1=length(y)}

		tpr[i]=nx1/nx
		fpr[i]=ny1/ny
		l[i]=lim
	}	
	v<-data.frame(fpr,tpr,l)

#GRAFICOS	
	plot(v$fpr,v$tpr,cex=0.5,pch=19,typ="b",xlab="FPR",ylab="TPR")
	x<-c(0,1)
	lines(x,x)
#/////////////////////////////////////////////////////////////////

	return(v)
}



#}}}
#dist_ang(x,y) calcula la distancia entre 2 posiciones x e y.
#{{{
dist_ang<-function(x,y){

}
#}}} 
#library('randomForestSRC')
#merclust(dat,cum=TRUE,gal=TRUE,rank=TRUE,sep=FALSE,relaxed=FALSE,nrank=100)
#dat debe ser un data frame con: id (numero),ra[° decimales],dec[° decimales],z (redshift),mag(aparente),color
#{{{
merclust<-function(dat,cum=TRUE,gal=TRUE,rank=TRUE,relaxed=TRUE,name.groups='estadisticos_grupos.dat',name.gal='estadisticos_galaxias.dat',rank.name='ranking.dat',relaxed.name='ranking_relaxed.dat',folder='folder',nrank=100,name_trainset_cum='trainset_cum.dat',name_trainset_gal='trainset_gal.dat',ntotal=0){
#Carpeta para guardar archivos
if(file.exists(folder)==FALSE){
  system(paste('mkdir',folder,sep=' '))
}
mc<-paste(folder,'/merging_clusters',sep='') 
if(file.exists(mc)==FALSE){
  system(paste('mkdir',mc,sep=' '))
}
name.groups<-paste(folder,'/',name.groups,sep='')
name.gal<-paste(folder,'/',name.gal,sep='')
rank.name<-paste(folder,'/',rank.name,sep='')
relaxed.name<-paste(folder,'/',relaxed.name,sep='')

#------------------------------------------------------------------------------
t1<-proc.time()
dat_gal=dat
if(cum == TRUE){
print('Empezando el calculo de estadisticos de los cumulos')
#Estadisticos cumulos
#{{{
dat=dat_gal
if(ntotal==0){
  ntotal=length(dat$ra)
}
ra_cum<-1:2
dec_cum<-1:2
z_cum<-1:2
delta<-1:2
delta_2<-1:2
ngroup<-1:2
ngal<-1:2
pval<-1:2
pval_2<-1:2
ind<-1:2
p_sw<-1:2
p_sf<-1:2
p_ad<-1:2
p_cvm<-1:2
p_lillie<-1:2
p_pearson<-1:2
color<-1:2
mag<-1:2
sph<-1:2
sph.dens<-1:2
gap.sph<-1:2
gamma<-1:2
c=0
#pb <- txtProgressBar(title = "progress bar", min = 0,max = ntotal, width = 82)
for(i in 1:ntotal){
#setTxtProgressBar(pb, i, label=paste( round(i/ntotal*100, 0),"% done"))
  if(length(dat$ra)>0){
print(i)
    groupid=dat$id[1]
    group<-subset(dat,dat$id==groupid)
    dat<-subset(dat,dat$id != groupid)

    ra=group$ra
    dec=group$dec
    vel=group$z*300000

    if(length(ra)>30){
      c=c+1    
      ngroup[c]=groupid
      ra_cum[c]=mean(ra)
      dec_cum[c]=mean(dec)
      z_cum[c]=mean(vel)/300000
      delta[c]=dressler(ra,dec,vel)
      delta_2[c]=dressler2(ra,dec,vel)
      ngal[c]=length(ra)
      pval[c]=mont(ra,dec,vel)
      pval_2[c]=pval[c]#mont2(ra,dec,vel) #No tiene sentido calcular y tarda mucho
      ind[c]=dressler.iter(ra,dec,vel) 
      shapiro.test(vel)->aux1
      p_sw[c]=aux1$p.val
      sf.test(vel)->aux1
      p_sf[c]=aux1$p.val
      ad.test(vel)->aux1
      p_ad[c]=aux1$p.val
      suppressWarnings(cvm.test(vel))->aux1
      p_cvm[c]=aux1$p.val
      lillie.test(vel)->aux1
      p_lillie[c]=aux1$p.val
      pearson.test(vel)->aux1
      p_pearson[c]=aux1$p.val
      gamma[c]= -99#relaxation(ra=ra,dec=dec,z=vel/300000,mag=group$mag)
      color[c]=mean(group$color)
      l<-sort(group$mag,decreasing=FALSE)
      mag[c]=l[1]
      gap.sph[c]=l[2]-l[1]
      sph[c]=sph.test(ra,dec)
      sph.dens[c]=dens.test(ra,dec)
    }

  }   
}

#close(pb)
grupos<-data.frame(ngroup,ra_cum,dec_cum,z_cum,delta,ngal,pval,ind,p_sw,p_sf,p_ad,p_cvm,p_lillie,p_pearson,delta_2,pval_2,color,mag,sph,sph.dens,gap.sph,gamma)
if(c==1){
  grupos<-grupos[1,]
}
vec<-c('ngroup','ra','dec','z','delta','ngal','pval','ind','p_sw','p_sf','p_ad','p_cvm','p_lillie','p_pearson','delta_2','pval_2','color','mag','sph','sph.dens','gap.sph','gamma')
dat_cum=grupos
write.table(grupos,file=name.groups,col.names=vec,row.names=FALSE)
#}}}
}
if(gal == TRUE){
print('Empezando el calculo de estadisticos de las galaxias')
#Estadisticos Galaxias
#{{{
dat_cum<-read.table(file=name.groups,header=TRUE)
dat=dat_gal
ntotal=length(dat_cum$ngroup)

c1=0
ra<-1:2
ngroup<-1:2
dec<-1:2
redshift<-1:2
mag_r<-1:2
g_r<-1:2
deltas_gal<-1:2
deltas2_gal<-1:2
del_cum<-1:2
pval_cum<-1:2
ngal<-1:2
sw_gal<-1:2 
sf_gal<-1:2
ad_gal<-1:2
cvm_gal<-1:2
lillie_gal<-1:2
pearson_gal<-1:2
sw_cum<-1:2  
sf_cum<-1:2 
ad_cum<-1:2 
cvm_cum<-1:2
lillie_cum<-1:2
pearson_cum<-1:2
col_cum<-1:2  
mag_cum<-1:2
ind_cum<-1:2

pb <- txtProgressBar(title = "progress bar", min = 0,max = ntotal, width = 82)
for(i in 1:ntotal){
setTxtProgressBar(pb, i, label=paste( round(i/ntotal*100, 0),"% done"))
  backh=dat_cum$ngroup[i]
  group<-subset(dat,dat$id==backh)
  ntodo=length(group$id)
  

  x=group$ra
  y=group$dec
  vel=group$z*300000

  delta<-dressler_graf(x,y,vel)
  delta2<-dressler_graf2(x,y,vel)

  p_sw<-1:ntodo
  p_sf<-1:ntodo
  p_ad<-1:ntodo
  p_cvm<-1:ntodo
  p_lillie<-1:ntodo
  p_pearson<-1:ntodo
  c=0
  for(j in 1:ntodo){
    x0=x[j]
    y0=y[j]
    r<-distancia1(x,x0,y,y0)
    vel_loc<-1:length(r)
    for(k in 1:length(r)){
      vel_loc[k]=vel[r[k]]
    }
  
    c=c+1 
    shapiro.test(vel_loc)->aux
    aux$p.value->p_sw[c]
    sf.test(vel_loc)->aux
    aux$p.value->p_sf[c]
    ad.test(vel_loc)->aux
    aux$p.value->p_ad[c]
    cvm.test(vel_loc)->aux
    aux$p.value->p_cvm[c]
    lillie.test(vel_loc)->aux
    aux$p.value->p_lillie[c]
    pearson.test(vel_loc)->aux
    aux$p.value->p_pearson[c]
  
  }

  for(j in 1:ntodo){
     c1=c1+1
     ngroup[c1]=backh
     ra[c1]=x[j]
     dec[c1]=y[j]
     redshift[c1]=vel[j]/300000
     mag_r[c1]=group$mag[j]
     g_r[c1]=group$color[j]
     deltas_gal[c1]=delta[j]
     deltas2_gal[c1]=delta2[j]
     del_cum[c1]=dat_cum$delta[i]
     pval_cum[c1]=dat_cum$pval[i]
     ngal[c1]=dat_cum$ngal[i]
     sw_gal[c1]=p_sw[j]
     sf_gal[c1]=p_sf[j]
     ad_gal[c1]=p_ad[j]
     cvm_gal[c1]=p_cvm[j]
     lillie_gal[c1]=p_lillie[j]
     pearson_gal[c1]=p_pearson[j]
     sw_cum[c1]=dat_cum$p_sw[i] 
     sf_cum[c1]=dat_cum$p_sf[i]
     ad_cum[c1]=dat_cum$p_ad[i]
     cvm_cum[c1]=dat_cum$p_cvm[i]
     lillie_cum[c1]=dat_cum$p_lillie[i]
     pearson_cum[c1]=dat_cum$p_pearson[i]
     col_cum[c1]=dat_cum$color[i] 
     mag_cum[c1]=dat_cum$mag[i]
     ind_cum[c1]=dat_cum$ind[i]

  }
}
close(pb)
vec<-c('ngroup','ra','dec','redshift','mag_r','g_r','deltas_gal','del_cum','pval_cum','ngal','sw_gal','sf_gal','ad_gal','cvm_gal','lillie_gal','pearson_gal','sw_cum','sf_cum','ad_cum','cvm_cum','lillie_cum','pearson_cum','col_cum','mag_cum','ind_cum','deltas2_gal')

galaxias<-data.frame(ngroup,ra,dec,redshift,mag_r,g_r,deltas_gal,del_cum,pval_cum,ngal,sw_gal,sf_gal,ad_gal,cvm_gal,lillie_gal,pearson_gal,sw_cum,sf_cum,ad_cum,cvm_cum,lillie_cum,pearson_cum,col_cum,mag_cum,ind_cum,deltas2_gal)
dat_gal=galaxias
write.table(galaxias,file=name.gal,col.names=vec,row.names=F,quote=F)
#}}}
}
if(rank== TRUE){
print('Empezando el rankiado de cumulos en merger')
#RANKING
#{{{
nrank=nrank

lim=0.3
lim.gal=0.4
#lim=0.2 #Version2
#lim.gal=0.3 #Version2

trainset<-read.table(file=name_trainset_cum,header=TRUE)
trainset_gal<-read.table(file=name_trainset_gal,header=TRUE)


dat_gals<-read.table(file=name.gal,header=TRUE)
if(length(dat_gals$rf.pred.gal) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
  dat_gals <- subset(dat_gals, select = -c(rf.pred.gal))
}
dat_cum<-read.table(file=name.groups,header=TRUE)
if(length(dat_cum$rf.pred.mer) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
  dat_cum <- subset(dat_cum, select = -c(rf.pred.mer))
}

cont.results=0
rf.pred.aux<-1:length(dat_cum$ngroup)
rf.pred.aux[]=0

pb <- txtProgressBar(title = "progress bar", min = 0,max = nrank, width = 82)
for(jj in 1:nrank){
setTxtProgressBar(pb, jj, label=paste( round(jj/nrank*100, 0),"% done"))
  
  rf.out<-suppressWarnings(randomForest(id_mer~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
  rf.pred.mer<-predict(rf.out,newdata=dat_cum)
  rf.pred.aux<-rf.pred.aux+rf.pred.mer
  
  if(length(dat_cum$rf.pred.mer) > 0){
    dat_cum$rf.pred.mer=rf.pred.mer
    dat<-dat_cum
  } else {
    dat<-data.frame(dat_cum,rf.pred.mer)
  }
  dat<-subset(dat,dat$rf.pred.mer>lim)
  
  
  rf.out.gal<-suppressWarnings(randomForest(id~deltas_gal*del_cum*g_r*sw_gal*col_cum,data=trainset_gal,importance=TRUE))
  rf.pred.gal<-predict(rf.out.gal,newdata=dat_gals)
  
  if(length(dat_gals$rf.pred.gal) > 0){
    dat_gals$rf.pred.gal=rf.pred.gal
    dat_gal<-dat_gals
  } else {
    dat_gal<-data.frame(dat_gals,rf.pred.gal)
  }
  dat_gal<-subset(dat_gal,dat_gal$rf.pred.gal>lim.gal)
  
  
  ncum=length(dat$delta)
  
  
  if(ncum > 0){ 
    aux<-1:(ncum*22)
    aux[1:(ncum*22)]=-99
    mat<-matrix(aux,ncol=22,nrow=ncum)
    for(i in 1:ncum){
      group<-subset(dat_gal,dat_gal$ngroup==dat$ngroup[i])
      if(length(group$ra)>10){
     
      ra=group$ra
      dec=group$dec
      vel=group$redshift*300000
      mag=group$mag_r
      color=group$g_r
      delta=group$deltas_gal
      peso=group$rf.pred.gal 
      
      delmin=0
      grupo.id<-paste(folder,'/merging_clusters/',toString(group$ngroup[1]),sep='')
      mixt.real(ra,dec,vel,mag,color,delta,peso,delmin,grupo.id)->v
      mat[i,1]=dat$ngroup[i]
      mat[i,2:22]=v
      }
    }    
    mat<-data.frame(mat)
    vec<-c('ngroup','l1','l2','rvir1','rvir2','dvel1','dvel2','mas1','mas2','par1','par2','tot','sigma1_ra','sigma1_dec','sigma2_ra','sigma2_dec','racen1','racen2','deccen1','deccen2','velcen1','velcen2')
    colnames(mat)<-vec
    mat<-subset(mat,mat$mas1>0)
    #mat<-subset(mat,mat$par1>0.22)
    #mat<-subset(mat,mat$par1>0.15)
   
    cont.results=cont.results+1
    if(cont.results==1){
      results=mat$ngroup
      results_mat=mat
    } else {
      results<-c(results,mat$ngroup)
      results_mat<-rbind(results_mat,mat)
    }
  } 
}
close(pb)


if(cont.results > 0){
  ntodo=length(results)
  name<-1:2
  largo<-1:2
  c=0
  for(i in 1:ntodo){
     if(length(results)>0){
        s1<-subset(results,results == results[1])
        c=c+1
        name[c]=s1[1]
        largo[c]=length(s1)
        results<-subset(results,results != results[1])
     }
  }
  
  largo=largo/nrank
  ngroup<-1:length(largo)
  l1<-1:length(largo)
  l2<-1:length(largo)
  l1_err<-1:length(largo)
  l2_err<-1:length(largo)
  rvir1<-1:length(largo)
  rvir2<-1:length(largo)
  rvir1_err<-1:length(largo)
  rvir2_err<-1:length(largo)
  dvel1<-1:length(largo)
  dvel2<-1:length(largo)
  dvel1_err<-1:length(largo)
  dvel2_err<-1:length(largo)
  m1<-1:length(largo)
  m2<-1:length(largo)
  par1<-1:length(largo)
  par2<-1:length(largo)
  par1_err<-1:length(largo)
  par2_err<-1:length(largo)
  tot<-1:length(largo)
  tot_err<-1:length(largo)
  sigma1_ra<-1:length(largo) 
  sigma2_ra<-1:length(largo) 
  sigma1_ra_err<-1:length(largo) 
  sigma2_ra_err<-1:length(largo) 
  sigma1_dec<-1:length(largo) 
  sigma2_dec<-1:length(largo) 
  sigma1_dec_err<-1:length(largo) 
  sigma2_dec_err<-1:length(largo) 
  m1_err<-1:length(largo)
  m2_err<-1:length(largo)
  ra1<-1:length(largo)
  ra2<-1:length(largo)
  ra1_err<-1:length(largo)
  ra2_err<-1:length(largo)
  dec1<-1:length(largo)
  dec2<-1:length(largo)
  dec1_err<-1:length(largo)
  dec2_err<-1:length(largo)
  z1<-1:length(largo)
  z2<-1:length(largo)
  z1_err<-1:length(largo)
  z2_err<-1:length(largo)
  
  for(i in 1:length(largo)){
    if(length(results_mat$ngroup)>0){
      mat_aux<-subset(results_mat,results_mat$ngroup==results_mat$ngroup[1])
      results_mat<-subset(results_mat,results_mat$ngroup!=results_mat$ngroup[1])
    
      ngroup[i]=mat_aux$ngroup[1]
      l1[i]=mean(mat_aux$l1) 
      l1_err[i]=sd(mat_aux$l1) 
      l2[i]=mean(mat_aux$l2) 
      l2_err[i]=sd(mat_aux$l2) 
      rvir1[i]=mean(mat_aux$rvir1) 
      rvir1_err[i]=sd(mat_aux$rvir1) 
      rvir2[i]=mean(mat_aux$rvir2) 
      rvir2_err[i]=sd(mat_aux$rvir2) 
      dvel1[i]=mean(mat_aux$dvel1) 
      dvel1_err[i]=sd(mat_aux$dvel1) 
      dvel2[i]=mean(mat_aux$dvel2) 
      dvel2_err[i]=sd(mat_aux$dvel2) 
      m1[i]=mean(mat_aux$mas1) 
      m1_err[i]=sd(mat_aux$mas1)
      m2[i]=mean(mat_aux$mas2) 
      m2_err[i]=sd(mat_aux$mas2)
      par1[i]=mean(mat_aux$par1) 
      par1_err[i]=sd(mat_aux$par1)
      par2[i]=mean(mat_aux$par2) 
      par2_err[i]=sd(mat_aux$par2)
      tot[i]=mean(mat_aux$tot)
      sigma1_ra[i]=mean(mat_aux$sigma1_ra) 
      sigma1_ra_err[i]=sd(mat_aux$sigma1_ra)
      sigma2_ra[i]=mean(mat_aux$sigma2_ra) 
      sigma2_ra_err[i]=sd(mat_aux$sigma2_ra)
      sigma1_dec[i]=mean(mat_aux$sigma1_dec) 
      sigma1_dec_err[i]=sd(mat_aux$sigma1_dec)
      sigma2_dec[i]=mean(mat_aux$sigma2_dec) 
      sigma2_dec_err[i]=sd(mat_aux$sigma2_dec)
      tot_err[i]=sd(mat_aux$tot)
      ra1[i]=mean(mat_aux$racen1)*180/pi
      ra1_err[i]=sd(mat_aux$racen1)*180/pi
      ra2[i]=mean(mat_aux$racen2)*180/pi
      ra2_err[i]=sd(mat_aux$racen2)*180/pi
      dec1[i]=mean(mat_aux$deccen1)*180/pi
      dec1_err[i]=sd(mat_aux$deccen1)*180/pi
      dec2[i]=mean(mat_aux$deccen2)*180/pi
      dec2_err[i]=sd(mat_aux$deccen2)*180/pi
      z1[i]=mean(mat_aux$velcen1)/300000
      z1_err[i]=sd(mat_aux$velcen1)/300000
      z2[i]=mean(mat_aux$velcen2)/300000
      z2_err[i]=sd(mat_aux$velcen2)/300000
    }
  }
  
  rank=largo
  rankin<-data.frame(ngroup,rank,m1,m1_err,ra1,ra1_err,dec1,dec1_err,z1,z1_err,m2,m2_err,ra2,ra2_err,dec2,dec2_err,z2,z2_err,l1,l1_err,l2,l2_err,tot,tot_err,rvir1,rvir1_err,rvir2,rvir2_err,dvel1,dvel1_err,dvel2,dvel2_err,par1,par1_err,par2,par2_err,sigma1_ra,sigma1_ra_err,sigma2_ra,sigma2_ra_err,sigma1_dec,sigma1_dec_err,sigma2_dec,sigma2_dec_err)

  if(c==1){
     rankin<-rankin[1,]
  }
} else {
  rankin<-'No hay cumulos en merger'
}
write.table(rankin,file=rank.name,row.names=FALSE)
rf.pred.mer<-rf.pred.aux/nrank
if(length(dat_cum$rf.pred.mer)>0){
  dat_cum$rf.pred.mer=rf.pred.mer
} else {
  dat_cum<-data.frame(dat_cum,rf.pred.mer)
}
write.table(dat_cum,file=name.groups,row.names=FALSE)
#}}}
}
if(relaxed==TRUE){
print('Empezando el rankiado de cumulos relajados')
#RANKING-RELAXED
#{{{
trainset<-read.table(file=name_trainset_cum,header=TRUE)
nrank=nrank

dat_cum<-read.table(file=name.groups,header=TRUE)
if(length(dat_cum$rf.pred.rel) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
  dat_cum <- subset(dat_cum, select = -c(rf.pred.rel))
}
lim=0.6

# Lectura de datos de cumulos
#if(mean(dat_cum$z) > 0.15){
#  trainset_rel<-read.table('/media/martin/store1/trabajos/mock/guo/mock_cumulos_snap56.dat',header=TRUE)
#} else {
#  trainset_rel<-read.table('/media/martin/store1/trabajos/mock/guo/mock_cumulos_snap63.dat',header=TRUE)
#}
cont.rel=0
rf.pred.aux<-1:length(dat_cum$ngroup)
rf.pred.aux[]=0

pb <- txtProgressBar(title = "progress bar", min = 0,max = nrank, width = 82)
for(jj in 1:nrank){
setTxtProgressBar(pb, jj, label=paste( round(jj/nrank*100, 0),"% done"))

  rf.out.rel<-suppressWarnings(randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
  #rf.out.rel<-randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie*sph*sph.dens*gap.sph,data=trainset,importance=TRUE)
  rf.pred.rel<-predict(rf.out.rel,newdata=dat_cum)
  rf.pred.aux<-rf.pred.aux+rf.pred.rel

  if(length(dat_cum$rf.pred.rel) > 0){
    dat_cum$rf.pred.rel=rf.pred.rel
    dat<-dat_cum
  } else {
    dat<-data.frame(dat_cum,rf.pred.rel)
  }
  dat<-subset(dat,dat$rf.pred.rel>lim)
 
  if(length(dat$ngroup)>0){
    cont.rel=cont.rel+1
    if(cont.rel==1){
      results=dat
    } else {
      results<-rbind(results,dat)
    }
  }
}
close(pb)

if(cont.rel>0){
  ntodo=length(dat_cum$ngroup)
  rank<-1:2
  for(i in 1:ntodo){
    aux<-subset(results,results$ngroup == dat_cum$ngroup[i])
    rank[i]=length(aux$ngroup)/nrank
  }
  rel_groups<-data.frame(dat_cum,rank)
  if(length(dat_cum$ngroup)==1){
    rel_groups<-rel_groups[1,]
  }
} else {
  rel_groups<-'No hay cumulos relajados'
}

write.table(rel_groups,file=relaxed.name,row.names=FALSE)
rf.pred.rel<-rf.pred.aux/nrank
if(length(dat_cum$rf.pred.rel)>0){
  dat_cum$rf.pred.rel=rf.pred.rel
} else {
  dat_cum<-data.frame(dat_cum,rf.pred.rel)
}
write.table(dat_cum,file=name.groups,row.names=FALSE)
#}}}
}
t2<-proc.time()
t<-(t2[3]-t1[3])
paste('El programa duro',toString(t/60),'minutos',sep=' ')->duracion
print(duracion)
}
#}}}
#mult.merger(group.id,method='distancia',...)
#{{{
#dat_gal<-read.table('estadisticos_galaxias.dat')
#mer<-read.table('ranking.dat')
mult.merger<-function(group.id,method='distancia',subfolder=0,dat_gal,mer,...){

   dat_gal<-subset(dat_gal,dat_gal$ngroup==group.id)
   dat_cum<-subset(mer,mer$ngroup==group.id)

   H0=73
   c=300000
#Metodo de separacion de las galaxias de las 2 subestructuras
if(method=='rvir'){
#Por radio virial
#{{{
   dist<-((dat_gal$ra*pi/180-dat_cum$ra1[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec1[1])**2)
   dist1=sqrt(dist)
   dist<-((dat_gal$ra*pi/180-dat_cum$ra2[1])**2)*cos(dat_cum$dec2[1])*cos(dat_cum$dec2[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec2[1])**2)
   dist2=sqrt(dist)
   aux<-data.frame(dat_gal,dist1,dist2)

   rvir1=dat_cum$rvir1[1]*H0/(c*dat_cum$z1[1])
   rvir2=dat_cum$rvir2[1]*H0/(c*dat_cum$z2[1])

   aux1<-subset(aux,aux$dist1 < rvir1)
   aux2<-subset(aux,aux$dist2 < rvir2)

   ra<-aux1$ra
   dec<-aux1$dec
   z<-aux1$redshift
   mag<-aux1$mag_r
   color<-aux1$g_r
   id<-1:length(aux1$ra)
   id[]=1
   dat1<-data.frame(ra,dec,z,mag,color,id)

   ra<-aux2$ra
   dec<-aux2$dec
   z<-aux2$redshift
   mag<-aux2$mag_r
   color<-aux2$g_r
   id<-1:length(aux2$ra)
   id[]=2

   dat2<-data.frame(ra,dec,z,mag,color,id)

#}}}
} else if(method=='distancia'){
#Por bisectriz o algo asi
#{{{
   dist<-((dat_gal$ra*pi/180-dat_cum$ra1[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec1[1])**2)
   dist1=sqrt(dist)
   dist<-((dat_gal$ra*pi/180-dat_cum$ra2[1])**2)*cos(dat_cum$dec2[1])*cos(dat_cum$dec2[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec2[1])**2)
   dist2=sqrt(dist)
   aux<-data.frame(dat_gal,dist1,dist2)

   dist_12<-((dat_cum$ra1[1]-dat_cum$ra2[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist_12=dist_12+((dat_cum$dec2[1]-dat_cum$dec1[1])**2)
   dist_12=sqrt(dist_12)/2

   aux1<-subset(aux,aux$dist1 < dist_12)
   aux2<-subset(aux,aux$dist2 < dist_12)

   ra<-aux1$ra
   dec<-aux1$dec
   z<-aux1$redshift
   mag<-aux1$mag_r
   color<-aux1$g_r
   id<-1:length(aux1$ra)
   id[]=1
   dat1<-data.frame(ra,dec,z,mag,color,id)

   ra<-aux2$ra
   dec<-aux2$dec
   z<-aux2$redshift
   mag<-aux2$mag_r
   color<-aux2$g_r
   id<-1:length(aux2$ra)
   id[]=2

   dat2<-data.frame(ra,dec,z,mag,color,id)

#}}}
} else if(method=='sustraccion'){
#Por sustraccion o algo asi
#{{{
   dist<-((dat_gal$ra*pi/180-dat_cum$ra1[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec1[1])**2)
   dist1=sqrt(dist)
   dist<-((dat_gal$ra*pi/180-dat_cum$ra2[1])**2)*cos(dat_cum$dec2[1])*cos(dat_cum$dec2[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec2[1])**2)
   dist2=sqrt(dist)
   aux<-data.frame(dat_gal,dist1,dist2)

   dist_12<-((dat_cum$ra1[1]-dat_cum$ra2[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist_12=dist_12+((dat_cum$dec2[1]-dat_cum$dec1[1])**2)
   dist_12=sqrt(dist_12)/2

   ind<-which(aux$dist1 < dist_12)
   aux1<-aux[-ind,]
   ind<-which(aux$dist2 < dist_12)
   aux2<-aux[-ind,]

   ra<-aux1$ra
   dec<-aux1$dec
   z<-aux1$redshift
   mag<-aux1$mag_r
   color<-aux1$g_r
   id<-1:length(aux1$ra)
   id[]=1
   dat1<-data.frame(ra,dec,z,mag,color,id)

   ra<-aux2$ra
   dec<-aux2$dec
   z<-aux2$redshift
   mag<-aux2$mag_r
   color<-aux2$g_r
   id<-1:length(aux2$ra)
   id[]=2

   dat2<-data.frame(ra,dec,z,mag,color,id)

#}}}
}
if(length(dat1$ra) > 30 & length(dat2$ra) > 30){
   dat<-rbind(dat1,dat2)
} else if(length(dat1$ra) > 30 & length(dat2$ra) < 30){
   dat<-dat1
} else if(length(dat1$ra) < 30 & length(dat2$ra) > 30){
   dat<-dat2
} 

flag=-1
if(length(dat$ra)>30){
#Analisis de las subestructuras por separado
   folder=paste(toString(group.id),'_subestructuras',sep='')
   if(subfolder != 0){
      folder<-paste(subfolder,'/',folder,sep='')
   }
   merclust(dat=dat,folder=folder,nrank=20)

   name<-paste(folder,'/ranking.dat',sep='')
   fl<-read.table(file=name,header=T)
   fl<-subset(fl,fl$rank>0.6)

   flag=0
   if(length(fl$rank)>0){
     if(length(fl$rank) == 2){
        flag=3
     } else {
        if(fl$ngroup[1]==1) {flag=1}
        if(fl$ngroup[1]==2) {flag=2}
     }
   }
}
   return(flag)
}
#}}}
#nose
#{{{
nose<-function(group.id,method='distancia'){
  dat_gal<-read.table('lemze_bien/estadisticos_galaxias.dat',header=T)
  mer<-read.table('ranking_v1.dat',header=T)

#------------------------------Graficos------------------------------------------
  dat_gal<-subset(dat_gal,dat_gal$ngroup==group.id)
  mer<-subset(mer,mer$ngroup==group.id)
  plot(dat_gal$ra,dat_gal$dec)
  points(mer$ra1[1]*180/pi,mer$dec1[1]*180/pi,pch=15,col='red')
  points(mer$ra2[1]*180/pi,mer$dec2[1]*180/pi,pch=15,col='blue')
#--------------------------------------------------------------------------------

  mult.merger(group.id=group.id,dat_gal=dat_gal,mer=mer,method=method)->flag
  folder=paste(toString(group.id),'_subestructuras',sep='')

  cont=1
  cont.fl3=0
  print(paste('flag= ',toString(flag),' in the ',toString(cont),' iteration',sep=''))
  while(flag > 0){
    if(flag==1){
       name<-paste(folder,'/estadisticos_galaxias.dat',sep='')
       dat_gal<-read.table(file=name,header=T)
       name<-paste(folder,'/ranking.dat',sep='')
       mer<-read.table(file=name,header=T)
       pl.sub(mer,pch=(cont+15))
       mult.merger(group.id=1,subfolder=folder,dat_gal=dat_gal,mer=mer,method=method)->flag
       folder=paste(folder,'/1_subestructuras',sep='')
    }
    if(flag==2){
       name<-paste(folder,'/estadisticos_galaxias.dat',sep='')
       dat_gal<-read.table(file=name,header=T)
       name<-paste(folder,'/ranking.dat',sep='')
       mer<-read.table(file=name,header=T)
       pl.sub(mer,pch=(cont+15))
       mult.merger(group.id=2,subfolder=folder,dat_gal=dat_gal,mer=mer,method=method)->flag
       folder=paste(folder,'/2_subestructuras',sep='')
    }
    if(flag==3){
       cont.fl3=cont.fl3+1
       name<-paste(folder,'/estadisticos_galaxias.dat',sep='')
       dat_gal<-read.table(file=name,header=T)
       name<-paste(folder,'/ranking.dat',sep='')
       mer<-read.table(file=name,header=T)
       pl.sub(mer,pch=(cont+15))
       mult.merger(group.id=1,subfolder=folder,dat_gal=dat_gal,mer=mer,method=method)->flag
       folder=paste(folder,'/1_subestructuras',sep='')
    }
    cont=cont+1
    print(paste('flag= ',toString(flag),' in the ',toString(cont),' iteration',sep=''))
  }

  print(paste('Hubo ',toString(cont.fl3),' flag = 3',sep=''))
  return(cont)
}
#}}}
#pl.sub<-function(mer,pch)
#{{{
pl.sub<-function(mer,pch){
  points(mer$ra1[1]*180/pi,mer$dec1[1]*180/pi,pch=pch,col='red')
  points(mer$ra2[1]*180/pi,mer$dec2[1]*180/pi,pch=pch,col='red')
  points(mer$ra1[2]*180/pi,mer$dec1[2]*180/pi,pch=pch,col='blue')
  points(mer$ra2[2]*180/pi,mer$dec2[2]*180/pi,pch=pch,col='blue')
}
#}}}
#train(dat,clas.mer,clas.rel,clas.gal,name='trainset_gal.dat',name.mer='trainset_cum.dat',name.rel='trainset_cum_rel.dat')
#{{{
train<-function(dat,clas.mer,clas.rel,clas.gal,name.gal='trainset_gal.dat',name.cum='trainset_cum.dat',folder=0){
#Carpeta para guardar archivos
if(folder!=0){
  system(paste('mkdir',folder,sep=' '))
  name.gal<-paste(folder,'/',name.gal,sep='')
  name.cum<-paste(folder,'/',name.cum,sep='')
}

#------------------------------------------------------------------------------
t1<-proc.time()
dat_gal=dat
if(cum == TRUE){
print('Empezando el calculo de estadisticos de los cumulos')
#Estadisticos cumulos
#{{{
dat=dat_gal
ntotal=length(dat$ra)

ra_cum<-1:2
dec_cum<-1:2
z_cum<-1:2
delta<-1:2
delta_2<-1:2
ngroup<-1:2
ngal<-1:2
pval<-1:2
pval_2<-1:2
ind<-1:2
p_sw<-1:2
p_sf<-1:2
p_ad<-1:2
p_cvm<-1:2
p_lillie<-1:2
p_pearson<-1:2
color<-1:2
mag<-1:2
sph<-1:2
sph.dens<-1:2
gap<-1:2
c=0
#pb <- txtProgressBar(title = "progress bar", min = 0,max = ntotal, width = 82)
for(i in 1:ntotal){
#setTxtProgressBar(pb, i, label=paste( round(i/ntotal*100, 0),"% done"))
  if(length(dat$ra)>0){
print(i)
    groupid=dat$id[1]
    group<-subset(dat,dat$id==groupid)
    dat<-subset(dat,dat$id != groupid)

    ra=group$ra
    dec=group$dec
    vel=group$z*300000

    if(length(ra)>30){
      c=c+1    
      ngroup[c]=groupid
      ra_cum[c]=mean(ra)
      dec_cum[c]=mean(dec)
      z_cum[c]=mean(vel)/300000
      delta[c]=dressler(ra,dec,vel)
      delta_2[c]=dressler2(ra,dec,vel)
      ngal[c]=length(ra)
      pval[c]=mont_par(ra,dec,vel)
      pval_2[c]=pval[c]#mont2(ra,dec,vel) #No tiene sentido calcular y tarda mucho
      ind[c]=dressler.iter(ra,dec,vel) 
      shapiro.test(vel)->aux1
      p_sw[c]=aux1$p.val
      sf.test(vel)->aux1
      p_sf[c]=aux1$p.val
      ad.test(vel)->aux1
      p_ad[c]=aux1$p.val
      cvm.test(vel)->aux1
      p_cvm[c]=aux1$p.val
      lillie.test(vel)->aux1
      p_lillie[c]=aux1$p.val
      pearson.test(vel)->aux1
      p_pearson[c]=aux1$p.val
      color[c]=mean(group$color)
      l<-sort(group$mag,decreasing=FALSE)
      mag[c]=l[1]
      gap[c]=l[2]-l[1]
      sph[c]=sph.test(ra,dec)
      sph.dens[c]=dens.test(ra,dec)
    }

  }   
}

#close(pb)
grupos<-data.frame(ngroup,ra_cum,dec_cum,z_cum,delta,ngal,pval,ind,p_sw,p_sf,p_ad,p_cvm,p_lillie,p_pearson,delta_2,pval_2,color,mag,sph,sph.dens,gap,clas.mer,clas.rel)
if(c==1){
  grupos<-grupos[1,]
}
vec<-c('ngroup','ra','dec','z','delta','ngal','pval','ind','p_sw','p_sf','p_ad','p_cvm','p_lillie','p_pearson','delta_2','pval_2','color','mag','sph','sph.dens','gap','id_mer','id_rel')
write.table(grupos,file=name.cum,col.names=vec,row.names=FALSE)
#}}}
}
if(gal == TRUE){
print('Empezando el calculo de estadisticos de las galaxias')
#Estadisticos Galaxias
#{{{
dat_cum<-read.table(file=name.groups,header=TRUE)
dat=dat_gal
ntotal=length(dat_cum$ngroup)

c1=0
ra<-1:2
ngroup<-1:2
dec<-1:2
redshift<-1:2
mag_r<-1:2
g_r<-1:2
deltas_gal<-1:2
deltas2_gal<-1:2
del_cum<-1:2
pval_cum<-1:2
ngal<-1:2
sw_gal<-1:2 
sf_gal<-1:2
ad_gal<-1:2
cvm_gal<-1:2
lillie_gal<-1:2
pearson_gal<-1:2
sw_cum<-1:2  
sf_cum<-1:2 
ad_cum<-1:2 
cvm_cum<-1:2
lillie_cum<-1:2
pearson_cum<-1:2
col_cum<-1:2  
mag_cum<-1:2
ind_cum<-1:2

pb <- txtProgressBar(title = "progress bar", min = 0,max = ntotal, width = 82)
for(i in 1:ntotal){
setTxtProgressBar(pb, i, label=paste( round(i/ntotal*100, 0),"% done"))
  backh=dat_cum$ngroup[i]
  group<-subset(dat,dat$id==backh)
  ntodo=length(group$id)
  

  x=group$ra
  y=group$dec
  vel=group$z*300000

  delta<-dressler_graf(x,y,vel)
  delta2<-dressler_graf2(x,y,vel)

  p_sw<-1:ntodo
  p_sf<-1:ntodo
  p_ad<-1:ntodo
  p_cvm<-1:ntodo
  p_lillie<-1:ntodo
  p_pearson<-1:ntodo
  c=0
  for(j in 1:ntodo){
    x0=x[j]
    y0=y[j]
    r<-distancia1(x,x0,y,y0)
    vel_loc<-1:length(r)
    for(k in 1:length(r)){
      vel_loc[k]=vel[r[k]]
    }
  
    c=c+1 
    shapiro.test(vel_loc)->aux
    aux$p.value->p_sw[c]
    sf.test(vel_loc)->aux
    aux$p.value->p_sf[c]
    ad.test(vel_loc)->aux
    aux$p.value->p_ad[c]
    cvm.test(vel_loc)->aux
    aux$p.value->p_cvm[c]
    lillie.test(vel_loc)->aux
    aux$p.value->p_lillie[c]
    pearson.test(vel_loc)->aux
    aux$p.value->p_pearson[c]
  
  }

  for(j in 1:ntodo){
     c1=c1+1
     ngroup[c1]=backh
     ra[c1]=x[j]
     dec[c1]=y[j]
     redshift[c1]=vel[j]/300000
     mag_r[c1]=group$mag[j]
     g_r[c1]=group$color[j]
     deltas_gal[c1]=delta[j]
     deltas2_gal[c1]=delta2[j]
     del_cum[c1]=dat_cum$delta[i]
     pval_cum[c1]=dat_cum$pval[i]
     ngal[c1]=dat_cum$ngal[i]
     sw_gal[c1]=p_sw[j]
     sf_gal[c1]=p_sf[j]
     ad_gal[c1]=p_ad[j]
     cvm_gal[c1]=p_cvm[j]
     lillie_gal[c1]=p_lillie[j]
     pearson_gal[c1]=p_pearson[j]
     sw_cum[c1]=dat_cum$p_sw[i] 
     sf_cum[c1]=dat_cum$p_sf[i]
     ad_cum[c1]=dat_cum$p_ad[i]
     cvm_cum[c1]=dat_cum$p_cvm[i]
     lillie_cum[c1]=dat_cum$p_lillie[i]
     pearson_cum[c1]=dat_cum$p_pearson[i]
     col_cum[c1]=dat_cum$color[i] 
     mag_cum[c1]=dat_cum$mag[i]
     ind_cum[c1]=dat_cum$ind[i]

  }
}
close(pb)
vec<-c('ngroup','ra','dec','redshift','mag_r','g_r','deltas_gal','del_cum','pval_cum','ngal','sw_gal','sf_gal','ad_gal','cvm_gal','lillie_gal','pearson_gal','sw_cum','sf_cum','ad_cum','cvm_cum','lillie_cum','pearson_cum','col_cum','mag_cum','ind_cum','deltas2_gal','id')

galaxias<-data.frame(ngroup,ra,dec,redshift,mag_r,g_r,deltas_gal,del_cum,pval_cum,ngal,sw_gal,sf_gal,ad_gal,cvm_gal,lillie_gal,pearson_gal,sw_cum,sf_cum,ad_cum,cvm_cum,lillie_cum,pearson_cum,col_cum,mag_cum,ind_cum,deltas2_gal,clas.gal)
dat_gal=galaxias
write.table(galaxias,file=name.gal,col.names=vec)
#}}}
}
t2<-proc.time()
t<-(t2[3]-t1[3])
paste('El programa duro',toString(t/60),'minutos',sep=' ')->duracion
print(duracion)
}
#}}}
#merclust.par(dat,cum=TRUE,gal=TRUE,rank=TRUE,sep=FALSE,relaxed=FALSE,nrank=100)
#dat debe ser un data frame con: id (numero),ra[° decimales],dec[° decimales],z (redshift),mag(aparente),color
#{{{
merclust.par<-function(dat,cum=TRUE,gal=TRUE,rank=TRUE,relaxed=TRUE,name.groups='estadisticos_grupos.dat',name.gal='estadisticos_galaxias.dat',rank.name='ranking.dat',relaxed.name='ranking_relaxed.dat',folder='folder',nrank=30,name_trainset_cum='trainset_cum.dat',name_trainset_gal='trainset_gal.dat',ntotal=0, ncluster = 8){
#Carpeta para guardar archivos
if(file.exists(folder)==FALSE){
  system(paste('mkdir',folder,sep=' '))
}
mc<-paste(folder,'/merging_clusters',sep='') 
if(file.exists(mc)==FALSE){
  system(paste('mkdir',mc,sep=' '))
}
name.groups<-paste(folder,'/',name.groups,sep='')
name.gal<-paste(folder,'/',name.gal,sep='')
rank.name<-paste(folder,'/',rank.name,sep='')
relaxed.name<-paste(folder,'/',relaxed.name,sep='')

#------------------------------------------------------------------------------
t1<-proc.time()
dat_gal=dat
if(cum == TRUE){
print('Empezando el calculo de estadisticos de los cumulos')
#Estadisticos cumulos
#{{{
dat=dat_gal
if(ntotal==0){
  ntotal=length(dat$ra)
}

ra_cum<-1:2
dec_cum<-1:2
z_cum<-1:2
delta<-1:2
delta_2<-1:2
ngroup<-1:2
ngal<-1:2
pval<-1:2
pval_2<-1:2
ind<-1:2
p_sw<-1:2
p_sf<-1:2
p_ad<-1:2
p_cvm<-1:2
p_lillie<-1:2
p_pearson<-1:2
color<-1:2
mag<-1:2
sph<-1:2
sph.dens<-1:2
gap.sph<-1:2
gamma<-1:2
c=0
#pb <- txtProgressBar(title = "progress bar", min = 0,max = ntotal, width = 82)
for(i in 1:ntotal){
#setTxtProgressBar(pb, i, label=paste( round(i/ntotal*100, 0),"% done"))
  if(length(dat$ra)>0){
print(i)
    groupid=dat$id[1]
    group<-subset(dat,dat$id==groupid)
    dat<-subset(dat,dat$id != groupid)

    ra=group$ra
    dec=group$dec
    vel=group$z*300000

    if(length(ra)>19){
      c=c+1    
      ngroup[c]=groupid
      ra_cum[c]=mean(ra)
      dec_cum[c]=mean(dec)
      z_cum[c]=mean(vel)/300000
      delta[c]=dressler(ra,dec,vel)
      delta_2[c]=dressler2(ra,dec,vel)
      ngal[c]=length(ra)
      pval[c]=mont_par(ra,dec,vel,ncluster)
      pval_2[c]=pval[c]#mont2(ra,dec,vel) #No tiene sentido calcular y tarda mucho
      ind[c]=dressler.iter(ra,dec,vel) 
      shapiro.test(vel)->aux1
      p_sw[c]=aux1$p.val
      sf.test(vel)->aux1
      p_sf[c]=aux1$p.val
      ad.test(vel)->aux1
      p_ad[c]=aux1$p.val
      cvm.test(vel)->aux1
      p_cvm[c]=aux1$p.val
      lillie.test(vel)->aux1
      p_lillie[c]=aux1$p.val
      pearson.test(vel)->aux1
      p_pearson[c]=aux1$p.val
      gamma[c]<- -99 #relaxation(ra=ra,dec=dec,z=vel/300000,mag=group$mag)
      color[c]=mean(group$color)
      l<-sort(group$mag,decreasing=FALSE)
      mag[c]=l[1]
      gap.sph[c]=l[2]-l[1]
      sph[c]=sph.test(ra,dec)
      sph.dens[c]=dens.test(ra,dec)
    }

  }   
}

#close(pb)
grupos<-data.frame(ngroup,ra_cum,dec_cum,z_cum,delta,ngal,pval,ind,p_sw,p_sf,p_ad,p_cvm,p_lillie,p_pearson,delta_2,pval_2,color,mag,sph,sph.dens,gap.sph,gamma)
if(c==1){
  grupos<-grupos[1,]
}
vec<-c('ngroup','ra','dec','z','delta','ngal','pval','ind','p_sw','p_sf','p_ad','p_cvm','p_lillie','p_pearson','delta_2','pval_2','color','mag','sph','sph.dens','gap.sph','gamma')
dat_cum=grupos
write.table(grupos,file=name.groups,col.names=vec,row.names=FALSE)
#}}}
}
if(gal == TRUE){
print('Empezando el calculo de estadisticos de las galaxias')
#Estadisticos Galaxias
#{{{
dat_cum<-read.table(file=name.groups,header=TRUE)
dat=dat_gal
ntotal=length(dat_cum$ngroup)

c1=0
ra<-1:2
ngroup<-1:2
dec<-1:2
redshift<-1:2
mag_r<-1:2
g_r<-1:2
deltas_gal<-1:2
deltas2_gal<-1:2
del_cum<-1:2
pval_cum<-1:2
ngal<-1:2
sw_gal<-1:2 
sf_gal<-1:2
ad_gal<-1:2
cvm_gal<-1:2
lillie_gal<-1:2
pearson_gal<-1:2
sw_cum<-1:2  
sf_cum<-1:2 
ad_cum<-1:2 
cvm_cum<-1:2
lillie_cum<-1:2
pearson_cum<-1:2
col_cum<-1:2  
mag_cum<-1:2
ind_cum<-1:2

pb <- txtProgressBar(title = "progress bar", min = 0,max = ntotal, width = 82)
for(i in 1:ntotal){
setTxtProgressBar(pb, i, label=paste( round(i/ntotal*100, 0),"% done"))
  backh=dat_cum$ngroup[i]
  group<-subset(dat,dat$id==backh)
  ntodo=length(group$id)
  

  x=group$ra
  y=group$dec
  vel=group$z*300000

  delta<-dressler_graf(x,y,vel)
  delta2<-dressler_graf2(x,y,vel)

  p_sw<-1:ntodo
  p_sf<-1:ntodo
  p_ad<-1:ntodo
  p_cvm<-1:ntodo
  p_lillie<-1:ntodo
  p_pearson<-1:ntodo
  c=0
  for(j in 1:ntodo){
    x0=x[j]
    y0=y[j]
    r<-distancia1(x,x0,y,y0)
    vel_loc<-1:length(r)
    for(k in 1:length(r)){
      vel_loc[k]=vel[r[k]]
    }
  
    c=c+1 
    shapiro.test(vel_loc)->aux
    aux$p.value->p_sw[c]
    sf.test(vel_loc)->aux
    aux$p.value->p_sf[c]
    ad.test(vel_loc)->aux
    aux$p.value->p_ad[c]
    cvm.test(vel_loc)->aux
    aux$p.value->p_cvm[c]
    lillie.test(vel_loc)->aux
    aux$p.value->p_lillie[c]
    pearson.test(vel_loc)->aux
    aux$p.value->p_pearson[c]
  
  }

  for(j in 1:ntodo){
     c1=c1+1
     ngroup[c1]=backh
     ra[c1]=x[j]
     dec[c1]=y[j]
     redshift[c1]=vel[j]/300000
     mag_r[c1]=group$mag[j]
     g_r[c1]=group$color[j]
     deltas_gal[c1]=delta[j]
     deltas2_gal[c1]=delta2[j]
     del_cum[c1]=dat_cum$delta[i]
     pval_cum[c1]=dat_cum$pval[i]
     ngal[c1]=dat_cum$ngal[i]
     sw_gal[c1]=p_sw[j]
     sf_gal[c1]=p_sf[j]
     ad_gal[c1]=p_ad[j]
     cvm_gal[c1]=p_cvm[j]
     lillie_gal[c1]=p_lillie[j]
     pearson_gal[c1]=p_pearson[j]
     sw_cum[c1]=dat_cum$p_sw[i] 
     sf_cum[c1]=dat_cum$p_sf[i]
     ad_cum[c1]=dat_cum$p_ad[i]
     cvm_cum[c1]=dat_cum$p_cvm[i]
     lillie_cum[c1]=dat_cum$p_lillie[i]
     pearson_cum[c1]=dat_cum$p_pearson[i]
     col_cum[c1]=dat_cum$color[i] 
     mag_cum[c1]=dat_cum$mag[i]
     ind_cum[c1]=dat_cum$ind[i]

  }
}
close(pb)
vec<-c('ngroup','ra','dec','redshift','mag_r','g_r','deltas_gal','del_cum','pval_cum','ngal','sw_gal','sf_gal','ad_gal','cvm_gal','lillie_gal','pearson_gal','sw_cum','sf_cum','ad_cum','cvm_cum','lillie_cum','pearson_cum','col_cum','mag_cum','ind_cum','deltas2_gal')

galaxias<-data.frame(ngroup,ra,dec,redshift,mag_r,g_r,deltas_gal,del_cum,pval_cum,ngal,sw_gal,sf_gal,ad_gal,cvm_gal,lillie_gal,pearson_gal,sw_cum,sf_cum,ad_cum,cvm_cum,lillie_cum,pearson_cum,col_cum,mag_cum,ind_cum,deltas2_gal)
dat_gal=galaxias
write.table(galaxias,file=name.gal,col.names=vec)
#}}}
}
if(rank== TRUE){
print('Empezando el rankiado de cumulos en merger')
#RANKING
#{{{
nrank=nrank

lim=0.3
lim.gal=0.4
#lim=0.2 #Version2
#lim.gal=0.3 #Version2

trainset<-read.table(file=name_trainset_cum,header=TRUE)
trainset_gal<-read.table(file=name_trainset_gal,header=TRUE)


dat_gals<-read.table(file=name.gal,header=TRUE)
if(length(dat_gals$rf.pred.gal) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
  dat_gals <- subset(dat_gals, select = -c(rf.pred.gal))
}
dat_cum<-read.table(file=name.groups,header=TRUE)
if(length(dat_cum$rf.pred.mer) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
  dat_cum <- subset(dat_cum, select = -c(rf.pred.mer))
}

cont.results=0

  
#cl <- makeCluster(3)
#registerDoParallel(cl)
cl <- makeSOCKcluster(ncluster)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=nrank, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

results_mat<-foreach(jj=1:nrank,.combine='rbind',.packages=c('randomForest','mclust','cosmoFns'), .export=c('mixt.real','gap','dist'),.options.snow=opts) %dopar% 
{
    rf.out<-suppressWarnings(randomForest(id_mer~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
    rf.pred.mer<-predict(rf.out,newdata=dat_cum)
    
    dat<-data.frame(dat_cum,rf.pred.mer)
    dat<-subset(dat,dat$rf.pred.mer > lim)
    
    
    rf.out.gal<-suppressWarnings(randomForest(id~deltas_gal*del_cum*g_r*sw_gal*col_cum,data=trainset_gal,importance=TRUE))
    rf.pred.gal<-predict(rf.out.gal,newdata=dat_gals)
    
    dat_gal<-data.frame(dat_gals,rf.pred.gal)
    dat_gal<-subset(dat_gal,dat_gal$rf.pred.gal>lim.gal)
    
    
    ncum=length(dat$delta)
    
    if(ncum > 0){ 
      aux<-1:(ncum*23)
      aux[1:(ncum*23)]=-99
      mat<-matrix(aux,ncol=23,nrow=ncum)
      for(i in 1:ncum){
        group<-subset(dat_gal,dat_gal$ngroup==dat$ngroup[i])
       
        ra=group$ra
        dec=group$dec
        vel=group$redshift*300000
        mag=group$mag_r
        color=group$g_r
        delta=group$deltas_gal
        peso=group$rf.pred.gal 
        
        delmin=0
        grupo.id<-paste(folder,'/merging_clusters/',toString(group$ngroup[1]),sep='')
        mixt.real(ra,dec,vel,mag,color,delta,peso,delmin,grupo.id)->v
        mat[i,1]=dat$ngroup[i]
        mat[i,2:22]=v
        mat[i,23]=dat$rf.pred.mer[i]
      }    
      mat<-data.frame(mat)
      vec<-c('ngroup','l1','l2','rvir1','rvir2','dvel1','dvel2','mas1','mas2','par1','par2','tot','sigma1_ra','sigma1_dec','sigma2_ra','sigma2_dec','racen1','racen2','deccen1','deccen2','velcen1','velcen2','rf.pred.mer')
      colnames(mat)<-vec
      mat<-subset(mat,mat$mas1>0)
   } 
   if(exists('mat') == FALSE){
     mat <- matrix(-5555, nrow=1,ncol=1)
   }
   return(mat)
}
close(pb)
stopCluster(cl)

rf.pred.mer<-1:length(dat_cum$ngroup)
rf.pred.mer[]=0
dat_cum <- data.frame(dat_cum, rf.pred.mer)

if(results_mat[1,1] != -5555){
  results<-results_mat$ngroup
  ntodo=length(results)
  name<-1:2
  largo<-1:2
  c=0
  for(i in 1:ntodo){
     if(length(results)>0){
        s1<-subset(results,results == results[1])
        c=c+1
        name[c]=s1[1]
        largo[c]=length(s1)
        results<-subset(results,results != results[1])
     }
  }
  
  largo=largo/nrank
  ngroup<-1:length(largo)
  rf.pred.mer<-1:length(largo)
  l1<-1:length(largo)
  l2<-1:length(largo)
  l1_err<-1:length(largo)
  l2_err<-1:length(largo)
  rvir1<-1:length(largo)
  rvir2<-1:length(largo)
  rvir1_err<-1:length(largo)
  rvir2_err<-1:length(largo)
  dvel1<-1:length(largo)
  dvel2<-1:length(largo)
  dvel1_err<-1:length(largo)
  dvel2_err<-1:length(largo)
  m1<-1:length(largo)
  m2<-1:length(largo)
  par1<-1:length(largo)
  par2<-1:length(largo)
  par1_err<-1:length(largo)
  par2_err<-1:length(largo)
  tot<-1:length(largo)
  tot_err<-1:length(largo)
  sigma1_ra<-1:length(largo) 
  sigma2_ra<-1:length(largo) 
  sigma1_ra_err<-1:length(largo) 
  sigma2_ra_err<-1:length(largo) 
  sigma1_dec<-1:length(largo) 
  sigma2_dec<-1:length(largo) 
  sigma1_dec_err<-1:length(largo) 
  sigma2_dec_err<-1:length(largo) 
  m1_err<-1:length(largo)
  m2_err<-1:length(largo)
  ra1<-1:length(largo)
  ra2<-1:length(largo)
  ra1_err<-1:length(largo)
  ra2_err<-1:length(largo)
  dec1<-1:length(largo)
  dec2<-1:length(largo)
  dec1_err<-1:length(largo)
  dec2_err<-1:length(largo)
  z1<-1:length(largo)
  z2<-1:length(largo)
  z1_err<-1:length(largo)
  z2_err<-1:length(largo)
  
  for(i in 1:length(largo)){
    if(length(results_mat$ngroup)>0){
      mat_aux<-subset(results_mat,results_mat$ngroup==results_mat$ngroup[1])
      ind_aux<-which(dat_cum$ngroup == results_mat$ngroup[1])
      results_mat<-subset(results_mat,results_mat$ngroup!=results_mat$ngroup[1])
    
      rf.pred.mer[i]<-mean(mat_aux$rf.pred.mer)
      dat_cum$rf.pred.mer[ind_aux]<-mean(mat_aux$rf.pred.mer)
      ngroup[i]=mat_aux$ngroup[1]
      l1[i]=mean(mat_aux$l1) 
      l1_err[i]=sd(mat_aux$l1) 
      l2[i]=mean(mat_aux$l2) 
      l2_err[i]=sd(mat_aux$l2) 
      rvir1[i]=mean(mat_aux$rvir1) 
      rvir1_err[i]=sd(mat_aux$rvir1) 
      rvir2[i]=mean(mat_aux$rvir2) 
      rvir2_err[i]=sd(mat_aux$rvir2) 
      dvel1[i]=mean(mat_aux$dvel1) 
      dvel1_err[i]=sd(mat_aux$dvel1) 
      dvel2[i]=mean(mat_aux$dvel2) 
      dvel2_err[i]=sd(mat_aux$dvel2) 
      m1[i]=mean(mat_aux$mas1) 
      m1_err[i]=sd(mat_aux$mas1)
      m2[i]=mean(mat_aux$mas2) 
      m2_err[i]=sd(mat_aux$mas2)
      par1[i]=mean(mat_aux$par1) 
      par1_err[i]=sd(mat_aux$par1)
      par2[i]=mean(mat_aux$par2) 
      par2_err[i]=sd(mat_aux$par2)
      tot[i]=mean(mat_aux$tot)
      sigma1_ra[i]=mean(mat_aux$sigma1_ra) 
      sigma1_ra_err[i]=sd(mat_aux$sigma1_ra)
      sigma2_ra[i]=mean(mat_aux$sigma2_ra) 
      sigma2_ra_err[i]=sd(mat_aux$sigma2_ra)
      sigma1_dec[i]=mean(mat_aux$sigma1_dec) 
      sigma1_dec_err[i]=sd(mat_aux$sigma1_dec)
      sigma2_dec[i]=mean(mat_aux$sigma2_dec) 
      sigma2_dec_err[i]=sd(mat_aux$sigma2_dec)
      tot_err[i]=sd(mat_aux$tot)
      ra1[i]=mean(mat_aux$racen1)*180/pi
      ra1_err[i]=sd(mat_aux$racen1)*180/pi
      ra2[i]=mean(mat_aux$racen2)*180/pi
      ra2_err[i]=sd(mat_aux$racen2)*180/pi
      dec1[i]=mean(mat_aux$deccen1)*180/pi
      dec1_err[i]=sd(mat_aux$deccen1)*180/pi
      dec2[i]=mean(mat_aux$deccen2)*180/pi
      dec2_err[i]=sd(mat_aux$deccen2)*180/pi
      z1[i]=mean(mat_aux$velcen1)/300000
      z1_err[i]=sd(mat_aux$velcen1)/300000
      z2[i]=mean(mat_aux$velcen2)/300000
      z2_err[i]=sd(mat_aux$velcen2)/300000
    }
  }
  
  rank=largo
  rankin<-data.frame(ngroup,rank,rf.pred.mer,m1,m1_err,ra1,ra1_err,dec1,dec1_err,z1,z1_err,m2,m2_err,ra2,ra2_err,dec2,dec2_err,z2,z2_err,l1,l1_err,l2,l2_err,tot,tot_err,rvir1,rvir1_err,rvir2,rvir2_err,dvel1,dvel1_err,dvel2,dvel2_err,par1,par1_err,par2,par2_err,sigma1_ra,sigma1_ra_err,sigma2_ra,sigma2_ra_err,sigma1_dec,sigma1_dec_err,sigma2_dec,sigma2_dec_err)

  if(c==1){
     rankin<-rankin[1,]
  }
} else {
  rankin<-'No hay cumulos en merger'
}

write.table(rankin,file=rank.name,row.names=FALSE)
write.table(dat_cum,file=name.groups,row.names=FALSE)
#}}}
}
if(relaxed==TRUE){
print('Empezando el rankiado de cumulos relajados')
#RANKING-RELAXED
#{{{
trainset<-read.table(file=name_trainset_cum,header=TRUE)
nrank=nrank

dat_cum<-read.table(file=name.groups,header=TRUE)
if(length(dat_cum$rf.pred.rel) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
  dat_cum <- subset(dat_cum, select = -c(rf.pred.rel))
}
lim=0.6

cont.rel=0
rf.pred.rel<-1:length(dat_cum$ngroup)
rf.pred.rel[]=0

#cl <- makeCluster(3)
#registerDoParallel(cl)
cl <- makeSOCKcluster(ncluster)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=nrank, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

results_mat<-foreach(jj=1:nrank,.combine='rbind',.packages=c('randomForest','mclust','cosmoFns'),.options.snow=opts) %dopar% 
{
  #rf.out.rel<-suppressWarnings(randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie*sph*sph.dens*gap.sph,data=trainset,importance=TRUE))
  rf.out.rel<-suppressWarnings(randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
  rf.pred.rel<-predict(rf.out.rel,newdata=dat_cum)

  dat<-data.frame(dat_cum,rf.pred.rel)
  dat<-subset(dat,dat$rf.pred.rel>lim)
 
  return(dat)
}
close(pb)
stopCluster(cl)

results <- results_mat
cont.rel <- length(results$rf.pred.rel)

if(cont.rel>0){
  ntodo=length(dat_cum$ngroup)
  rank<-1:2
  for(i in 1:ntodo){
    aux<-subset(results,results$ngroup == dat_cum$ngroup[i])
    rank[i]=length(aux$ngroup)/nrank
    if(rank[i] > 0){rf.pred.rel[i]<-mean(aux$rf.pred.rel)}
  }
  rel_groups<-data.frame(dat_cum,rank)
  if(length(dat_cum$ngroup)==1){
    rel_groups<-rel_groups[1,]
  }
} else {
  rel_groups<-'No hay cumulos relajados'
}

write.table(rel_groups,file=relaxed.name,row.names=FALSE)
dat_cum<-data.frame(dat_cum,rf.pred.rel)
write.table(dat_cum,file=name.groups,row.names=FALSE)
#}}}
}
t2<-proc.time()
t<-(t2[3]-t1[3])
paste('El programa duro',toString(t/60),'minutos',sep=' ')->duracion
print(duracion)
}
#}}}
#relaxation(ra,dec,z,mag) ANDA MAL!!!!!
#{{{
relaxation<-function(ra,dec,z,mag){

#Kernel Gaussiano
  x<-ra-mean(ra)
  y<-dec-mean(dec)
  l<-luminosity(mag,z)
  dat<-data.frame(x,y,l)

  nbin=100
  rmin=0.5
  dx=(max(x)-min(x))/nbin
  dy=(max(y)-min(y))/nbin

  I<-matrix(0,nrow=nbin,ncol=nbin)
  R<-matrix(0,nrow=nbin,ncol=nbin)
  Tita<-matrix(0,nrow=nbin,ncol=nbin)

  for(i in 1:nbin){
    xmin=min(x)+dx*(i-1)
    xmax=xmin+dx
    for(j in 1:nbin){
      ymin=min(y)+dy*(j-1)
      ymax=ymin+dy
      R[i,j]=sqrt((xmin+dx/2)**2+(ymin+dy/2)**2)

      if((xmin+dx/2) > 0 & (ymin+dx/2) > 0){ # 1° cuadrante
        Tita[i,j]=atan2((ymin+dy/2),(xmin+dx/2)) 
      } else if((xmin+dx/2) < 0 & (ymin+dx/2) > 0){ # 2° cuadrante
        Tita[i,j]=atan2((ymin+dy/2),(xmin+dx/2))
      } else if((xmin+dx/2) < 0 & (ymin+dx/2) < 0){ # 3° cuadrante
        Tita[i,j]=atan2((ymin+dy/2),(xmin+dx/2))+2*pi
      } else if((xmin+dx/2) > 0 & (ymin+dx/2) < 0){ # 4° cuadrante
        Tita[i,j]=atan2((ymin+dy/2),(xmin+dx/2))+2*pi
      }

      r<-(dat$x-(xmin+dx/2))**2+(dat$y-(ymin+dy/2))**2
      sub<-data.frame(r,dat$l)
      aux<-subset(sub,sqrt(sub$r) < rmin)

      if(length(aux$r) > 0){
        g<-kernel.gauss(r=aux$r,sigma=0.1)
        I[i,j]=sum(aux$dat.l*g)
      }
    }
  }

#Asymmetry factor

  S2=sum(I*I)
  delta=0

  for(i in 1:nbin){
    ii=nbin-i+1
    for(j in 1:nbin){
      jj=nbin-j+1
      delta=delta+((I[i,j]-I[ii,jj])**2)/2
    }
  }

  alfa=delta/S2

#Ridge factor
ntita=50

  for(k in 1:ntita){
    tita.min=2*pi*(k-1)/ntita
    tita.max=2*pi*k/ntita
    r0[k]=0
    r<-1:2
    l<-1:2
    c=0
    for(i in 1:nbin){
      for(j in 1:nbin){
        if(Tita[i,j] < tita.max & Tita[i,j] > tita.min){
          c=c+1
          r[c]=R[i,j]
          l[c]=I[i,j]
        }
      }
    }
    sub<-data.frame(r,l)
    #rmin=r[which.max(sub$l)]
    #sub<-subset(sub,sub$r > rmin)
    if(length(sub$r) > 0){
      nls(l~I0/(1+(r/r0)**2),data=sub,start=list(I0=mean(sub$l),r0=mean(sub$r)),nls.control(maxiter=100,warnOnly=TRUE))->fit
      r0[k]=coef(fit)[2]
    }
  }

  beta=mean(r0)/min(r0)

#Normalized deviation
nrad=50
drad=max(R)/ntita
r<-1:nrad
l<-1:nrad
  for(k in 1:nrad){
    rad.min=(k-1)/ntita
    rad.max=rad.min+drad
    r0[k]=0
    laux<-1:2
    c=0
    for(i in 1:nbin){
      for(j in 1:nbin){
        if(R[i,j] < rad.max & R[i,j] > rad.min){
          c=c+1
          laux[c]=I[i,j]
        }
      }
    }
    r[k]=rad.min+drad/2
    l[k]=mean(laux)
  }
  sub<-data.frame(r,l)
  #rmin=r[which.max(sub$l)]
  #sub<-subset(sub,sub$r > rmin)
  if(length(sub$r) > 0){
    nls(l~I0/(1+(r/r0)**2),data=sub,start=list(I0=mean(sub$l),r0=mean(sub$r)),nls.control(maxiter=100,warnOnly=TRUE))->fit
    r0=coef(fit)[2]
    I0=coef(fit)[1]
  }

  delta=0
  for(i in 1:nbin){
    for(j in 1:nbin){
       delta=delta+(I[i,j]-I0/(1+(R[i,j]/r0)**2))**2
    }
  }
  delta=delta/S2

  gamma=beta-1.9*alfa-3.58*delta-0.10
  return(gamma)
}
#}}}
#luminosity(mag,z)
#{{{
luminosity<-function(mag,z){
  M_sun=4.42 #banda r
  d=D.L(z)*(10**6) #Pc
  M=mag-5*log10(d/10)
  lum<-10**((M_sun-M)/2.5) #L_sun
  return(lum)
}
#}}}
#kernel.gauss(r,sigma)
#{{{
kernel.gauss<-function(r,sigma){
  g<-(1/(2*pi*sigma))*exp(-r/(2*sigma*sigma))
  return(g)
}
#}}}
#graf3d(ngrupo,movie=FALSE,xray=...,subs=...) Grafico 3d 
#{{{
#library('rgl')
graf3d<-function(name,movie=FALSE,xray,subs){
 

  mat<-read.table(name,header=T)
  z<-mat$vel/300000
  mat<-data.frame(mat,z)

  open3d()
  mat1<-subset(mat,mat$id == '1')
  mat1<-subset(mat1,select=c('ra','dec','z','id','mag'))
  plot3d(x=mat1$ra,y=mat1$dec,z=(mat1$z),radius=0.1*exp(mat1$mag)/exp(max(mat$mag)),type='s',col='red',box=FALSE,xlab='ra',ylab='dec',zlab='z')
  mat1<-subset(mat1,select=c('ra','dec','z'))
  cov1=cov(mat1)
  mean1=c(mean(mat1$ra),mean(mat1$dec),mean(mat1$z))
  plot3d(ellipse3d(cov1,centre=mean1),col='red',alpha=0.5,add=TRUE)
  
  
  mat2<-subset(mat,mat$id == '2')
  mat2<-subset(mat2,select=c('ra','dec','z','id','mag'))
  plot3d(x=mat2$ra,y=mat2$dec,z=mat2$z,radius=0.1*exp(mat2$mag)/exp(max(mat$mag)),type='s',add=TRUE,col='green')
  mat2<-subset(mat2,select=c('ra','dec','z'))
  cov2=cov(mat2)
  mean2=c(mean(mat2$ra),mean(mat2$dec),mean(mat2$z))
  plot3d(ellipse3d(cov2,centre=mean2),col='green',alpha=0.5,add=TRUE)
  
  play3d(spin3d(axis = c(1, 0, 0)), duration = 3)
  play3d(spin3d(axis = c(0, -1, 0)), duration = 1)
  par3d(scale=aspect3d(x=1,y=1,z=1))
#OPTIONAL ARGUMENTS
  if(missing(xray) == FALSE){
    zmin=(600+min(vel_tot))/1000 
    zmax=(600+max(vel_tot))/1000
    for(i in 1:length(xray$ra)){
       rgl.lines(x=c(xray$ra[i],xray$ra[i]),y=c(xray$dec[i],xray$dec[i]),z=c(zmin,zmax),lwd=5,col='black')
    }
  }

 if(missing(subs) == FALSE){
    zmin=(600+min(vel_tot))/1000 
    zmax=(600+max(vel_tot))/1000
    for(i in 1:length(subs$ra)){
       rgl.lines(x=c(subs$ra[i],subs$ra[i]),y=c(subs$dec[i],subs$dec[i]),z=c(zmin,zmax),lwd=5,col='magenta')
    }
  }

  if(movie == TRUE){
    movie_name<-readline(prompt='Write the file name for the cluster movie')
    movie3d(spin3d(axis = c(0, -1, 0)), duration = 12,dir = getwd())
 
    move<-paste(movie_name,'.gif',sep='')
    move<-paste('mv movie.gif',move,sep=' ')
    system(move)
  }
}
#}}}
#mult.merger(group.id,method='distancia',...)
#{{{
#dat_gal<-read.table('estadisticos_galaxias.dat')
#mer<-read.table('ranking.dat')
mult.merger<-function(group.id,method='distancia',subfolder=0,dat_gal,mer,...){

   dat_gal<-subset(dat_gal,dat_gal$ngroup==group.id)
   dat_cum<-subset(mer,mer$ngroup==group.id)

   H0=73
   c=300000
#Metodo de separacion de las galaxias de las 2 subestructuras
if(method=='rvir'){
#Por radio virial
#{{{
   dist<-((dat_gal$ra*pi/180-dat_cum$ra1[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec1[1])**2)
   dist1=sqrt(dist)
   dist<-((dat_gal$ra*pi/180-dat_cum$ra2[1])**2)*cos(dat_cum$dec2[1])*cos(dat_cum$dec2[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec2[1])**2)
   dist2=sqrt(dist)
   aux<-data.frame(dat_gal,dist1,dist2)

   rvir1=dat_cum$rvir1[1]*H0/(c*dat_cum$z1[1])
   rvir2=dat_cum$rvir2[1]*H0/(c*dat_cum$z2[1])

   aux1<-subset(aux,aux$dist1 < rvir1)
   aux2<-subset(aux,aux$dist2 < rvir2)

   ra<-aux1$ra
   dec<-aux1$dec
   z<-aux1$redshift
   mag<-aux1$mag_r
   color<-aux1$g_r
   id<-1:length(aux1$ra)
   id[]=1
   dat1<-data.frame(ra,dec,z,mag,color,id)

   ra<-aux2$ra
   dec<-aux2$dec
   z<-aux2$redshift
   mag<-aux2$mag_r
   color<-aux2$g_r
   id<-1:length(aux2$ra)
   id[]=2

   dat2<-data.frame(ra,dec,z,mag,color,id)

#}}}
} else if(method=='distancia'){
#Por bisectriz o algo asi
#{{{
   dist<-((dat_gal$ra*pi/180-dat_cum$ra1[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec1[1])**2)
   dist1=sqrt(dist)
   dist<-((dat_gal$ra*pi/180-dat_cum$ra2[1])**2)*cos(dat_cum$dec2[1])*cos(dat_cum$dec2[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec2[1])**2)
   dist2=sqrt(dist)
   aux<-data.frame(dat_gal,dist1,dist2)

   dist_12<-((dat_cum$ra1[1]-dat_cum$ra2[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist_12=dist_12+((dat_cum$dec2[1]-dat_cum$dec1[1])**2)
   dist_12=sqrt(dist_12)/2

   aux1<-subset(aux,aux$dist1 < dist_12)
   aux2<-subset(aux,aux$dist2 < dist_12)

   ra<-aux1$ra
   dec<-aux1$dec
   z<-aux1$redshift
   mag<-aux1$mag_r
   color<-aux1$g_r
   id<-1:length(aux1$ra)
   id[]=1
   dat1<-data.frame(ra,dec,z,mag,color,id)

   ra<-aux2$ra
   dec<-aux2$dec
   z<-aux2$redshift
   mag<-aux2$mag_r
   color<-aux2$g_r
   id<-1:length(aux2$ra)
   id[]=2

   dat2<-data.frame(ra,dec,z,mag,color,id)

#}}}
} else if(method=='sustraccion'){
#Por sustraccion o algo asi
#{{{
   dist<-((dat_gal$ra*pi/180-dat_cum$ra1[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec1[1])**2)
   dist1=sqrt(dist)
   dist<-((dat_gal$ra*pi/180-dat_cum$ra2[1])**2)*cos(dat_cum$dec2[1])*cos(dat_cum$dec2[1])
   dist=dist+((dat_gal$dec*pi/180-dat_cum$dec2[1])**2)
   dist2=sqrt(dist)
   aux<-data.frame(dat_gal,dist1,dist2)

   dist_12<-((dat_cum$ra1[1]-dat_cum$ra2[1])**2)*cos(dat_cum$dec1[1])*cos(dat_cum$dec1[1])
   dist_12=dist_12+((dat_cum$dec2[1]-dat_cum$dec1[1])**2)
   dist_12=sqrt(dist_12)/2

   ind<-which(aux$dist1 < dist_12)
   aux1<-aux[-ind,]
   ind<-which(aux$dist2 < dist_12)
   aux2<-aux[-ind,]

   ra<-aux1$ra
   dec<-aux1$dec
   z<-aux1$redshift
   mag<-aux1$mag_r
   color<-aux1$g_r
   id<-1:length(aux1$ra)
   id[]=1
   dat1<-data.frame(ra,dec,z,mag,color,id)

   ra<-aux2$ra
   dec<-aux2$dec
   z<-aux2$redshift
   mag<-aux2$mag_r
   color<-aux2$g_r
   id<-1:length(aux2$ra)
   id[]=2

   dat2<-data.frame(ra,dec,z,mag,color,id)

#}}}
}
if(length(dat1$ra) > 30 & length(dat2$ra) > 30){
   dat<-rbind(dat1,dat2)
} else if(length(dat1$ra) > 30 & length(dat2$ra) < 30){
   dat<-dat1
} else if(length(dat1$ra) < 30 & length(dat2$ra) > 30){
   dat<-dat2
} 

flag=-1
if(length(dat$ra)>30){
#Analisis de las subestructuras por separado
   folder=paste(toString(group.id),'_subestructuras',sep='')
   if(subfolder != 0){
      folder<-paste(subfolder,'/',folder,sep='')
   }
   merclust(dat=dat,folder=folder,nrank=20)

   name<-paste(folder,'/ranking.dat',sep='')
   fl<-read.table(file=name,header=T)
   fl<-subset(fl,fl$rank>0.6)

   flag=0
   if(length(fl$rank)>0){
     if(length(fl$rank) == 2){
        flag=3
     } else {
        if(fl$ngroup[1]==1) {flag=1}
        if(fl$ngroup[1]==2) {flag=2}
     }
   }
}
   return(flag)
}
#}}}
#nose
#{{{
nose<-function(group.id,method='distancia'){
  dat_gal<-read.table('lemze_bien/estadisticos_galaxias.dat',header=T)
  mer<-read.table('ranking_v1.dat',header=T)

#------------------------------Graficos------------------------------------------
  dat_gal<-subset(dat_gal,dat_gal$ngroup==group.id)
  mer<-subset(mer,mer$ngroup==group.id)
  plot(dat_gal$ra,dat_gal$dec)
  points(mer$ra1[1]*180/pi,mer$dec1[1]*180/pi,pch=15,col='red')
  points(mer$ra2[1]*180/pi,mer$dec2[1]*180/pi,pch=15,col='blue')
#--------------------------------------------------------------------------------

  mult.merger(group.id=group.id,dat_gal=dat_gal,mer=mer,method=method)->flag
  folder=paste(toString(group.id),'_subestructuras',sep='')

  cont=1
  cont.fl3=0
  print(paste('flag= ',toString(flag),' in the ',toString(cont),' iteration',sep=''))
  while(flag > 0){
    if(flag==1){
       name<-paste(folder,'/estadisticos_galaxias.dat',sep='')
       dat_gal<-read.table(file=name,header=T)
       name<-paste(folder,'/ranking.dat',sep='')
       mer<-read.table(file=name,header=T)
       pl.sub(mer,pch=(cont+15))
       mult.merger(group.id=1,subfolder=folder,dat_gal=dat_gal,mer=mer,method=method)->flag
       folder=paste(folder,'/1_subestructuras',sep='')
    }
    if(flag==2){
       name<-paste(folder,'/estadisticos_galaxias.dat',sep='')
       dat_gal<-read.table(file=name,header=T)
       name<-paste(folder,'/ranking.dat',sep='')
       mer<-read.table(file=name,header=T)
       pl.sub(mer,pch=(cont+15))
       mult.merger(group.id=2,subfolder=folder,dat_gal=dat_gal,mer=mer,method=method)->flag
       folder=paste(folder,'/2_subestructuras',sep='')
    }
    if(flag==3){
       cont.fl3=cont.fl3+1
       name<-paste(folder,'/estadisticos_galaxias.dat',sep='')
       dat_gal<-read.table(file=name,header=T)
       name<-paste(folder,'/ranking.dat',sep='')
       mer<-read.table(file=name,header=T)
       pl.sub(mer,pch=(cont+15))
       mult.merger(group.id=1,subfolder=folder,dat_gal=dat_gal,mer=mer,method=method)->flag
       folder=paste(folder,'/1_subestructuras',sep='')
    }
    cont=cont+1
    print(paste('flag= ',toString(flag),' in the ',toString(cont),' iteration',sep=''))
  }

  print(paste('Hubo ',toString(cont.fl3),' flag = 3',sep=''))
  return(cont)
}
#}}}
#pl.sub<-function(mer,pch)
#{{{
pl.sub<-function(mer,pch){
  points(mer$ra1[1]*180/pi,mer$dec1[1]*180/pi,pch=pch,col='red')
  points(mer$ra2[1]*180/pi,mer$dec2[1]*180/pi,pch=pch,col='red')
  points(mer$ra1[2]*180/pi,mer$dec1[2]*180/pi,pch=pch,col='blue')
  points(mer$ra2[2]*180/pi,mer$dec2[2]*180/pi,pch=pch,col='blue')
}
#}}}
