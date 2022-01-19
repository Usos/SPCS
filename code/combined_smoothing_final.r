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

findSNeighbors<-function(positions, is.hexa=F , tau.s=2, nThreads=8)
{
  if(is.hexa)
    tau.s<-tau.s*2

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

fillBlanks<-function(exp, coord, is.hexa=F, tau.s=2, filling.threshold = 0.5, nThreads = 8)
{

  x.bounds<-c(min(coord$st.x),max(coord$st.x))
  y.bounds<-c(min(coord$st.y),max(coord$st.y))
  x.availible<-x.bounds[1]:x.bounds[2]
  y.availible<-y.bounds[1]:y.bounds[2]
  x.availible<-sort(rep(x.availible,length(y.availible)))
  coord.maybe<-data.frame(st.x=x.availible,st.y=y.availible)

  if(is.hexa)
  {
    tau.s<-tau.s*2
    coord.maybe<- coord.maybe %>% 
      mutate(is.exist = ((st.x %% 2) == (st.y %% 2))) %>%
      select(st.x,st.y)
   }
  
  coord.maybe<-suppressMessages(coord %>% right_join(coord.maybe,copy = T) %>% filter(is.na(barcodes)))
  coord.maybe$barcodes<-paste0("SUPP",1:nrow(coord.maybe))
  
  cl=makeCluster(nThreads)
  registerDoParallel(cl)
  
  s.neighbors<-foreach(i=1:nrow(coord.maybe),.packages = c("dplyr"),.export = c("GetPointNeighbors")) %dopar%
    {
      nbs.maybe<-GetPointNeighbors(px = coord.maybe$st.x[i],py = coord.maybe$st.y[i],order = tau.s)
      if(is.hexa)
        nbs.maybe<-nbs.maybe %>% filter(0 == distance %% 2)
      
      nbs<-coord %>% inner_join(nbs.maybe,by=c("st.x"="x","st.y"="y"))
      if(nrow(nbs)>filling.threshold * nrow(nbs.maybe))
      {
        nbs %>%
        mutate(s.contribution=(1/distance)) %>%
        mutate(ss.contribution=s.contribution/(sum(s.contribution)+0.00000001))
      } else
      NULL
    }
  
  names(s.neighbors)<-coord.maybe$barcodes
  s.neighbors<-s.neighbors[sapply(s.neighbors,function(x){!is.null(x)})]
  
  if(0==length(s.neighbors))
  {
    s.exp<-NULL
    s.neighbors<-NULL
  }
  else
  {
    s.exp<-foreach(i=1:length(s.neighbors),.combine = "cbind",.packages = c("dplyr","Matrix")) %dopar%
      {
        exp[,s.neighbors[[i]]$barcodes] %*% Matrix(s.neighbors[[i]]$ss.contribution,sparse=T)
      }
    colnames(s.exp)<-names(s.neighbors)
  }
  stopCluster(cl)
  
  return(list(exp=s.exp,colData=filter(coord.maybe,barcodes %in% names(s.neighbors))))
}

#Run SPCS
SPCS<-function(data.tbl, coor.tbl, is.hexa = F, is.padding =T, tau.p=16, tau.s=2, alpha=0.6, beta=0.4, filling.thres=0.5, nThreads=8)
{
  data.mat<-Matrix(as.matrix(data.tbl))
  
  message(paste0("SPCS Smoothing on ",ncol(data.tbl)," spots with ",nrow(data.tbl)," genes"))
  message("1) Pre-processing...",appendLF = F)
  
  t0<-proc.time()

  rpca.res = calcRandPCA(t(data.mat)) 
  pdist.mat = calcPDist(rpca.res)  
  rmat = calcPContribution(pdist.mat)  
  
  positions = coor.tbl
  colnames(positions) = c("st.x","st.y")
  positions$barcodes = rownames(positions)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).")
  
  message("2) Detecting neighbors...\n"," --Pattern neighbors...",appendLF = F)
  t0<-proc.time()

  p.neighbors <- findPNeighbors(rmat, positions, tau.p, nThreads)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).\n"," --Spatial neighbors...",appendLF = F)
  t0<-proc.time()

  s.neighbors <- findSNeighbors(positions,is.hexa, tau.s, nThreads)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).")
  
  message("3) Calculating expression...",appendLF = F)
  t0<-proc.time()

  data.tbl.combinedsmooth = getCSmoothedExp(data.mat, s.neighbors, p.neighbors, alpha, beta, nThreads)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).")
  
  if(is.padding)
  {
    message("4) Filling blank spots...",appendLF = F)
    t0<-proc.time()

    filled.blanks<-fillBlanks(data.tbl.combinedsmooth, positions, is.hexa, tau.s, filling.thres, nThreads)

    t1<-proc.time()
    message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).\n",
        nrow(filled.blanks$colData),"blank spots has been filled.")
    data.tbl.combinedsmooth = as.data.frame(as.matrix(cbind(data.tbl.combinedsmooth,filled.blanks$exp)))
    attr(data.tbl.combinedsmooth,"SUPP_ColData")<-filled.blanks$colData
  }
  else
  {
    data.tbl.combinedsmooth = as.data.frame(as.matrix(data.tbl.combinedsmooth))
  }
  
  
  if (sum(is.na(data.tbl.combinedsmooth)) > 0) {
    mean.exp = sum(data.tbl.combinedsmooth[!is.na(data.tbl.combinedsmooth)]) / sum(!is.na(data.tbl.combinedsmooth))
    data.tbl.combinedsmooth[is.na(data.tbl.combinedsmooth)] = mean.exp
  }
  
  return(data.tbl.combinedsmooth)
}