##Modified version of the Eva KF Chan Funtion for Fst calculation (http://www.evachan.org/rscripts.html)
#input files:
# 1) geno:  numeric matrix representing a subset of two populations (i.e. one pairwise comparison)
# 2) subpop:  a character string containing the population membership for each sample.
wc_hapfst <- function(geno, subpop) {
  spop <- unique(as.character(subpop)) #isolate the population identifiers
  r <- length(spop) #calculate the number of populations
  n0 <- n1 <- matrix(NA, ncol=r, nrow=nrow(geno)) #create two matrices that are X SNPs long by Y populations wide
  
  #for each population
  for(i in 1:r) {
    inds <- which(subpop == spop[i]) #identify the member samples
    n0[,i] <- rowSums(geno[,inds]==0,na.rm=T) #calculate the number of member samples with the REF allele
    n1[,i] <- rowSums(geno[,inds]==1,na.rm=T) #calculate the number of member samples with the ALT allele
  }
  
  #This ensuing block uses the count matrices to directly calculate haploid Fst
  #The formula is written in a form to that follows the parameters provided in
  #Weir's 1996 Genetic Data Analysis II: Methods for Discrete Population Genetic Data
  
  ni <- n0 + n1
  
  n_bar <- rowSums(ni,na.rm=T) / r
  
  nc <- ((r * n_bar) - (rowSums((ni * ni),na.rm=T) / (r * n_bar))) / (r - 1)
  
  pi_tilda <- n0/ni
  
  p_bar <-rowSums(((ni * pi_tilda) / (r * n_bar)), na.rm=T )
  
  s_square <- rowSums( ((ni * ((pi_tilda - p_bar)^2)) / ((r - 1) * n_bar)), na.rm=T )
  
  T1<-s_square - ((1 / (n_bar - 1)) * ((p_bar * (1 - p_bar)) - (((r - 1) / r) * s_square)))
  
  T2<-(((nc - 1)/(n_bar - 1)) * (p_bar * (1 - p_bar))) + ((1 + (((r - 1) * (n_bar - nc)) / (n_bar - 1))) * (s_square / r))
  
  theta_hat<-T1/T2
  
  return(theta_hat)
}