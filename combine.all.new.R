##load in data in format gname/-log10(p.motif)/-log10(p.dhs)
motif.and.dnase.stat3.test1 <- read.delim("~/data/stat3_mel1_intersect_two_psorted", header=F)
motif.and.dnase.stat3.test2 <- read.delim("~/data/intersect_stat3_fs_1_two_psorted", header=F)
motif.and.dnase.stat3.test3 <- read.delim("~/data/stat3_mel3_intersect_two_psorted", header=F)

##calculate pvals from log transformed pvals
load.pvals <- function(file){
p1 <- 10**(-1*(file[,2]))##motif p vals
names(p1) <- toupper(file[,1])
n.names <- names(p1)
p2 <- 10**(-1*(file[,3]))##dhs peak_calling p vals
names(p2) <- toupper(file[,1])
return(list(p1,p2))
}

# load.p2 <- function(file){
#   p2 <- 10**(-1*(file[,3]))##dhs peak_calling p vals
#   names(p2) <- toupper(file[,1])
#   #n.names <- names(p2)
#   return(p2)
# }
# p.dhs1 <- load.p2(motif.and.dnase.stat3.test1)

p.motif1 <- unlist(load.pvals(motif.and.dnase.stat3.test1)[1])
p.dhs1 <- unlist(load.pvals(motif.and.dnase.stat3.test1)[2])

p.motif2 <- unlist(load.pvals(motif.and.dnase.stat3.test2)[1])
p.dhs2 <- unlist(load.pvals(motif.and.dnase.stat3.test2)[2])

p.motif3 <- unlist(load.pvals(motif.and.dnase.stat3.test3)[1])
p.dhs3 <- unlist(load.pvals(motif.and.dnase.stat3.test3)[2])

source('~/code/combine.test.R')

####################
##avgpvals function to combine p vals from multiple motifs/dhs for one gene
####################
avgpvals <- function(pvals){
j <- 1
k <- 1
nr <- length(unique(c(names(pvals))))
v <- matrix(NA, nrow = nr, ncol = length(unique(c(names(pvals)))))
#v <- matrix(NA)
unique.names <- unique(c(names(pvals)))
rownames(v) <- unique.names
v[k,j] <- pvals[1]
#names(v[k,j]) <- names(m_d_network_stat3[1])
pname <- names(pvals[1])
for (i in 2:length(pvals)) {
  if (k < nr) {
    if (names(pvals[i]) == pname) {
      #names(v[k]) <- names(m_d_network_stat3[i])
      #v[k] <- (v[k])*(m_d_network_stat3[i])
      #v[k,] <- c(v[k,],rep(NA,1))
      v[k,j+1] <- pvals[i]
      #names(v[k,j+1]) <- names(m_d_network_stat3[i])
      j <- j+1
    } else {
      j <- 1
      k <- k+1
      pname <- names(pvals[i])
      v[k,j] <- pvals[i]
      #names(v[k,j]) <- names(m_d_network_stat3[i])
    }
  }
}  
p.combined <- rep(NA, nr)
for (i in 1:nr){
  nc <- which(!is.na(v[i,]))
  pc <- rep(NA, length(nc))
  for (j in 1:length(nc)){     
    pc[j] <- v[i,j]
  } 
  #p.combined[i] <- combine.test(p = pc, method = "z.transform")##combine p vals using z.transform
  #p.combined[i] <- (-1)*log10(min(pc))##take the best p val
  #p.combined[i] <- (-1)*log10(exp(mean(log(pc))))##geometric mean and take the -log10 of pvals
  #p.combined[i] <- exp(mean(log(pc)))
  p.combined[i] <- prod(pc)
  }
names(p.combined) <- names(v[,1])
return(p.combined)
}

p.motif <- avgpvals(p.motif2)
p.dhs <- avgpvals(p.dhs2)

#######################
##combine motif and dhs p vals with fisher's methods
#######################

np <- list()
np.log <- list()
p.dhs.log <- list()
for (i in 1:length(p.motif)){
  np[i] <- combine.test(p = c(p.motif[i],p.dhs[i]), method = "z.transform")##combine motif p-vals with dhs peak calling p-vals
  np.log[i] <- (-1)*log10(as.numeric(np[i]))
  p.dhs.log[i] <- (-1)*log10(as.numeric(p.dhs[i]))
}

m_d_network_stat3_combined <- as.numeric(np.log) #as.numeric(np)
names(m_d_network_stat3_combined) <- names(p.dhs)
#pname <- "NONAME"

#######################
##use the combined p vals to calculate AUPR
#######################
#m_d_network_stat3_combined <- p.combined
#names(m_d_network_stat3_combined) <- names(p.combined)

gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"
knockout_rna_zscores.all <- read.delim(gs.file)
k_r_network_stat3 <- knockout_rna_zscores.all[,"STAT3"]
names(k_r_network_stat3) <- toupper(rownames(knockout_rna_zscores.all))
k_r_network_stat3 <- k_r_network_stat3[which(k_r_network_stat3>0)]

#unique.names <- unique(c( names(m_d_network_stat3), names(k_r_network_stat3)))
unique.names <- unique(c(names(m_d_network_stat3_combined), names(k_r_network_stat3)))

overlap <- (length(m_d_network_stat3_combined) + length(k_r_network_stat3)) - length(unique.names)

data <- matrix(NA, nr=length(unique.names), nc=2)
rownames(data) <- unique.names
colnames(data) <- c("motif_dnase", "knockout_rna")
for (i in 1:nrow(data)) {
  data[i,1] <- m_d_network_stat3_combined[unique.names[i]]
  data[i,2] <- k_r_network_stat3[unique.names[i]]
}
#data <- data[which(!is.na(data))] #<- 0
#data[which(is.na(data))] <- 0 #1000

source('~/code/AUPR/aupr.R')
gs <- data[,2]
gs[which(is.na(gs))] <- 0 
#gs <- gs[!is.na]
pred <- data[,1]
pred <- pred[which(!is.na(pred))]
#pred <- pred[!is.na]
#colnames(data)
#pred <- pred>0
# gs>0
gs <- gs>0
# gs
calcAupr(pred,gs)

stat3_aupr <- calcAupr(pred,gs)
plot(stat3_aupr$rec, stat3_aupr$prec, type = 'lines', main = 'Stat3_faire_ztransform')
save.image('~/data/combine.RData')