#SPCS ST smoothing method
#Code by Yusong Liu in HEU

#If any part of this code is used in publishable documents, please cite:
#    Yusong Liu, Tongxin Wang, and Ben Duggan et. al, 
#          SPCS: A Spatial and Pattern Combined Smoothing Method for Spatial Transcriptomic Expression (2021)

require("Matrix")
require("rsvd")
require("factoextra")
require("foreach")
require("dplyr")
require("doParallel")

#Pattern neighbors

calcRandPCA<-function(tmat, d=10, seed=42)
{
  # @param tmat A non-negative matrix with samples by features
  # @return A matrix with features by samples projected on PCA space
  set.seed(seed)
  
  rpca_obj <- rpca(tmat, k=d, center=T, scale=F, retx=T, p=10, q=7)
  rpca_obj$x
}

calcPDist <- function(tmat)
{
  get_dist(tmat,method = "pearson")
}

calcPContribution<-function(distMat)
{
  distMat<-as.matrix(distMat)
  rmat<-exp(-4.5*(distMat)^2)
  rmat[distMat>1]<-0
  diag(rmat)<-0
  
  return(Matrix(rmat, sparse=T))
}

findPNeighbors<-function(mat.pcontribution, positions, tau.p=16, nThreads=8)
{
  if(tau.p < 0)
    stop("Parameter tau.p must be 0 or positive.")
  
  bars<-positions$barcodes
  ind.bars<-1:length(bars)
  names(ind.bars)<-bars
  
  cl=makeCluster(nThreads)
  registerDoParallel(cl)
  
  p.neighbors<-foreach(i=1:length(bars),.packages = c("dplyr","Matrix")) %dopar%
    {
      a<-positions[mat.pcontribution[,i]>0,] 
      if(nrow(a)>0){
        b<-a %>% 
          mutate(p.contribution = mat.pcontribution[ind.bars[barcodes],i])
        if(0 == tau.p || nrow(a)<tau.p)
          b %>% mutate(ps.contribution = p.contribution/(sum(p.contribution)+0.00000001))
        else
          (b %>% arrange(desc(p.contribution)))[1:tau.p,] %>% mutate(ps.contribution = p.contribution/(sum(p.contribution)+0.00000001))
      } else
        NULL
    }
  
  stopCluster(cl)
  
  names(p.neighbors)<-bars
  
  return(p.neighbors)
  
}

#Spatial neighbors

GetPointNeighbors<-function(px,py,order=2,is.keepself=F)
{

  x.range<-(px-order):(px+order)
  x.range<-x.range[x.range>=0]
  y.range<-(py-order):(py+order)
  y.range<-y.range[y.range>=0]
  
  x.range<-rep.int(x.range,times=length(y.range))
  x.range<-sort(x.range)
  
  points<-data.frame(x=x.range, y=y.range)
  
  points<-points %>% mutate(distance=abs(x-px)+abs(y-py)) %>%
    filter(distance <= order)
  
  if(!is.keepself)
    points<-points %>% filter(distance>0)
  
  return(points)
  
}

findSNeighbors<-function(positions, tau.s=2, nThreads=8)
{

  bars<-positions$barcodes
  ind.bars<-1:length(bars)
  names(ind.bars)<-bars
  
  cl=makeCluster(nThreads)
  registerDoParallel(cl)
  
  s.neighbors<-foreach(i=1:length(bars),.packages =c("dplyr","Matrix"),.export = c("GetPointNeighbors")) %dopar%
    {
      a<-inner_join(x=positions, y=GetPointNeighbors(px=positions$st.x[i],py=positions$st.y[i],order = tau.s),by=c("st.x"="x","st.y"="y")) 
      if(nrow(a)>0){
        a %>%
        mutate(s.contribution=(1/distance)) %>%
        mutate(ss.contribution=s.contribution/(sum(s.contribution)+0.00000001))
      } else
        NULL
    } 
  
  stopCluster(cl)
  
  names(s.neighbors)<-positions$barcodes
  
  return(s.neighbors)
  
}


#get combined smoothed expression matrix

getCSmoothedExp<-function(exp, s.neighbors, p.neighbors, alpha, beta, nThreads=8)
{
  bars<-names(p.neighbors)
  ind.bars<-1:length(p.neighbors)
  names(ind.bars)<-bars
  
  cl=makeCluster(nThreads)
  registerDoParallel(cl)
  
  s.exp<-foreach(i=1:length(bars),.combine = "cbind",.packages = c("dplyr","Matrix")) %dopar%
    {
      if(!is.null(s.neighbors[[i]]))
        exp[,ind.bars[s.neighbors[[i]]$barcodes]] %*% Matrix(s.neighbors[[i]]$ss.contribution,sparse=T)
      else
        0
    }
  
  p.exp<-foreach(i=1:length(bars),.combine = "cbind",.packages = c("dplyr","Matrix")) %dopar%
    {
      if(!is.null(p.neighbors[[i]]))
        exp[,ind.bars[p.neighbors[[i]]$barcodes]] %*% Matrix(p.neighbors[[i]]$ps.contribution,sparse=T)
      else
        0
    }
  
  stopCluster(cl)
  
  exp.smoothed<-exp * (1-alpha) +
    (s.exp * beta + p.exp * (1-beta))*alpha

  return(exp.smoothed)
}

#Run SPCS
SPCS<-function(data.tbl, coor.tbl, tau.p=16, tau.s=2, alpha=0.6, beta=0.4, nThreads=8)
{
  data.mat<-Matrix(as.matrix(data.tbl))
  
  rpca.res = calcRandPCA(t(data.mat)) 
  pdist.mat = calcPDist(rpca.res)  
  rmat = calcPContribution(pdist.mat)  
  
  positions = coor.tbl
  colnames(positions) = c("st.x","st.y")
  positions$barcodes = rownames(positions)
  
  p.neighbors <- findPNeighbors(rmat, positions, tau.p, nThreads)
  s.neighbors <- findSNeighbors(positions, tau.s, nThreads)
  
  data.tbl.combinedsmooth = getCSmoothedExp(data.mat, s.neighbors, p.neighbors, alpha, beta, nThreads)
  data.tbl.combinedsmooth = as.data.frame(as.matrix(data.tbl.combinedsmooth))
  
  if (sum(is.na(data.tbl.combinedsmooth)) > 0) {
    mean.exp = sum(data.tbl.combinedsmooth[!is.na(data.tbl.combinedsmooth)]) / sum(!is.na(data.tbl.combinedsmooth))
    data.tbl.combinedsmooth[is.na(data.tbl.combinedsmooth)] = mean.exp
  }
  
  return(data.tbl.combinedsmooth)
}
