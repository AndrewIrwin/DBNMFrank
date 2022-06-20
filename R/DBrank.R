# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#################get the best rank k NMF results in sz.ini different initial values
devfn=function(data,k,sz.ini,hw=F,distn){
  data.nmf=data[,which(colSums(data)!=0)]
  devnmf=1e17
  sd.1=rep(1,sz.ini)
  if(distn=="Normal"){
    for (i in 1:sz.ini){
      nmf=NMF::nmf(t(data.nmf),k,'lee')
      if(NMF::deviance(nmf,method="euclidean")<devnmf){
        if(hw==T){
          h=NMF::basis(nmf)
          w=NMF::coef(nmf)
          sd.1[i]=sqrt(sum((t(w)%*%t(h)-data.nmf)^2)/nrow(data.nmf)/ncol(data.nmf))
        }
        else{
          h=NULL
          w=NULL
        }
        devnmf=NMF::deviance(nmf,method="euclidean")
      }
    }
    sd=mean(sd.1)
  }
  if(distn=="Poisson"){
    for (i in 1:sz.ini){
      nmf=NMF::nmf(t(data.nmf),k,'brunet')
      if(NMF::deviance(nmf,method="KL")<devnmf){
        if(hw==T){
          h=NMF::basis(nmf)
          w=NMF::coef(nmf)
        }
        else{
          h=NULL
          w=NULL
        }
        devnmf=NMF::deviance(nmf,method="KL")
      }
    }
    sd=NULL
  }
  return(list(devnmf=devnmf,h1=h,w1=w,sd1=sd))
}
#############generate bootstrap sample
bst.sample=function(np, nr, un, sd = NULL) {
  if (!is.null(sd)) {
    data.bst = t(matrix(stats::rnorm(np*nr, t(un), sd), np, nr))
    data.bst[which(data.bst<0)]=0
  } else {
    data.bst = t(matrix(stats::rpois(np*nr, t(un)), np, nr))
  }
  return(data.bst = data.bst)
}

###########the hyphothesis test of rank k vs. rank larger than k.
DBtest=function(data,kt,sz,distn,inisz){
  data.rm=data[,colSums(data)!=0]
  np=ncol(data.rm)
  nr=nrow(data.rm)
  dev1=devfn(data.rm,kt-1,inisz,hw=T,distn)
  dev2=devfn(data.rm,kt,inisz,distn=distn)
  dev=dev1$devnmf-dev2$devnmf
  h1=dev1$h1
  w1=dev1$w1
  sd1=dev1$sd1

  devnoibd=rep(0,sz)
  for(j in 1:sz){
    ############mean
    un=t(w1)%*%t(h1)
data.bst=bst.sample(np,nr,un,sd1)
    nmf0=devfn(data.bst,kt-1,1,distn=distn)
    nmf1=devfn(data.bst,kt,1,distn=distn)
    devnoibd[j]=nmf0$devnmf-nmf1$devnmf
  }

  ##################get weight matrix
  ############mean
  data.bst.err=bst.sample(np,nr,un,sd1)

  dev0all=dev1all=rep(0,sz)############try 500

  for(j in 1:inisz){#############try 500
    nmf0=devfn(data.bst.err,kt-1,1,distn=distn)
    nmf1=devfn(data.bst.err,kt,1,distn=distn)
    dev0all[j]=nmf0$devnmf
    dev1all[j]=nmf1$devnmf
  }
  e0=dev0all-min(dev0all)
  e1=dev1all-min(dev1all)
  #  pure error sample
  error=e0[expand.grid(1:sz,1:sz)[,1]]-e1[expand.grid(1:sz,1:sz)[,2]]
  # get convolved bootstrap null distribution
  rst=pmledecon::pmledecon(devnoibd,error,bsz=25)
  pvalue2=sum((dev<rst$sup)*rst$f)*(rst$sup[2]-rst$sup[1])
  return(pvalue2=pvalue2)
}

#' @export
DBrank=function(data,k=1,alpha=0.1,distn="Poisson",sz=50,inisz=50){
  i=k;pvalue2=0
  while(pvalue2<alpha){
    i=i+1
    k2=i
    rst2=DBtest(data,k2,sz,distn,inisz)
    pvalue2=rst2
    print(k2-1)
    print(pvalue2)
  }
  return(list(rank=k2-1,pvalue=pvalue2))
}

