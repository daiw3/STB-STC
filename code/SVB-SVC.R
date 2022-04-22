rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
blk <- as.numeric(args[1])
level <- as.character(args[2])
method <- as.character(args[3])

# ---------------------------------------  Specify the data directly below ----------------------------------------
data.dir = "STB-STC/data/"
load(paste0(data.dir, level, "_", method, "_blk", blk, ".RData")

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

source("STB-STC/code/0-tarv-transform.R") # this transforms the block of OTUs into super-taxon
source("STB-STC/code/0-di_cal_blk.R")) # this calculates depth importance score for each OTU
source("STB-STC/code/0-super_var_form.R")) # this selects cut-off to form super-taxon for STB
source("STB-STC/code/0-super_var_form_mean.R")) # this selects cut-off to form super-taxon for STC


discovery_id1 = read.table(paste0(data.dir, "discovery_id1.txt"), header = T)[,1]
discovery_id2 = read.table(paste0(data.dir, "discovery_id2.txt"), header = T)[,1]
discovery_id = c(discovery_id1, discovery_id2)

# extract OTU data -----------------------------------------------
geno <- geno.divided[, as.character(discovery_id), drop = FALSE]
member <- as.character(discovery_id)
pheno <- pheno[pheno[, "member"] %in% discovery_id, ]
pheno$age <- scale(pheno$age)

# contact DI-TARV based analysis ---------------------------------------------
# calculate depth importance score with random forest only using discovery_id1
geno1 = geno[, as.character(discovery_id1), drop = FALSE]
if((sum(geno1!=0) == 0) & (dim(geno1)[1] == 1)) next # stop("discovery set  only has 1 OTU and all values are 0.")

pheno1 = pheno[pheno$member %in% discovery_id1, ]
start.time = Sys.time()
geno1 = as.data.frame(t(geno1))
if(method == "STC"){
  otu_remove = which(apply(geno1,2,sum) == 0)
  if(length(otu_remove)!=0){
    geno1 = geno1[, -otu_remove]
  }   
} 

num.otus <- dim(geno1)[2]
all.otus <- colnames(geno1)

#---------------------------------------------------------------------------------
############# Ranking ##############
#---------------------------------------------------------------------------------
Z.local <- di_cal(pheno1, geno1, num.trees = 500, num.var.boot = dim(geno1)[2]/2) # marginal effect calculation


#---------------------------------------------------------------------------------
############# Best Cut-off Selection ##############
#---------------------------------------------------------------------------------


if(method == "STC"){
  
  # threshold selection
  geno_ordered <- geno1[,order(Z.local[,1],decreasing = T), drop = FALSE]
  super.var.form <- super_var_form_mean(pheno1, geno_ordered, pheno1$member, file_prefix)
  
}else{
  
  annotation <- rep(file_prefix, dim(geno1)[2])
  geno1 = as.data.frame(t(geno1))
  geno.transed.blk <- TARV_transform(geno1, annotation, Z.local[, 1], direction="dominant")
  geno.transed <- geno.transed.blk
  colnames(geno.transed) <- paste0(file_prefix, c("_agg"), c("+"))
  
  super.var.form <- super_var_form(pheno1, geno.transed, pheno1$member)
}

#---------------------------------------------------------------------------------
############# Marginal Association Test ##############
#---------------------------------------------------------------------------------


## extract information of contributing OTUs of each super-taxon
super.variant.cand <- list()
for(r in seq(dim(super.var.form)[1])){
  if(!is.na(super.var.form[r, 4])){
    z.ord <- order(Z.local[, 1], decreasing = T)
    # select the top SNPs according to the cut-off value
    otus.selected <- all.otus[z.ord]
    otus.selected <- otus.selected[seq(super.var.form[r, 4])]
    otus.selected = cbind(otus.selected, super.var.form$super.var[r])
    otus.selected = as.data.frame(otus.selected)
    colnames(otus.selected) <- c("otu.name", "super.var")
    row.names(otus.selected) <- NULL
    
    super.variant.cand[[r]] <- otus.selected
  }
}


## Super-taxon formation: calculate p-values based on cut-off on the second discovery sample
cutoff_otus <- do.call(rbind, super.variant.cand[which(!is.na(super.var.form[, 4]))])
super.var.otus[[file_prefix]] <- cutoff_otus

out.tmp = cutoff_otus$otu.name
geno2 = geno[out.tmp,as.character(discovery_id2), drop = FALSE]
pheno2 <- pheno[pheno$member %in% discovery_id2, ]

if(method  == "STC"){
  X <- NULL
  num.block <- as.character(unique(cutoff_otus$super.var))
  X.tmp <- NULL
  for(super.blk in num.block){
    trans.mode <- substring(super.blk, nchar(super.blk), nchar(super.blk))
    if(trans.mode == "+"){
      s <- apply(geno2, 2, sum)
      X.gene <- rep(NA, dim(geno2)[2])
      X.gene[s!=0] <- apply(geno2[, s!= 0, drop = FALSE], 2, function(x) sum(x)/sum(x != 0))
      X.gene[is.na(X.gene)] <- 0
    }else{
      s <- apply(geno2, 2, sum)
      X.gene <- rep(NA, dim(geno2)[2])
      X.gene[s!=0] <- apply(geno2[, s!= 0, drop = FALSE], 2, function(x) sum(x)/sum(x != 0))
      X.gene[is.na(X.gene)] <- 0
    }
    X.gene <- data.frame(member=colnames(geno2), X.gene, stringsAsFactors = F)
    colnames(X.gene)[2] <- super.blk
    if(is.null(X.tmp)){
      X.tmp <- X.gene
    }else{
      X.tmp <- merge(X.tmp, X.gene, by="member")
    }
  }
  
  if(is.null(X)){
    X <- X.tmp
  }else{
    X<- merge(X, X.tmp, by="member")
  }
}else{
  X <- NULL
  num.block <- as.character(unique(cutoff_otus$super.var))
  X.tmp <- NULL
  for(super.blk in num.block){
    trans.mode <- substring(super.blk, nchar(super.blk), nchar(super.blk))
    X.gene <- as.numeric(apply(geno2 > 1, 2, any))
    X.gene[is.na(X.gene)] <- 0
    X.gene <- data.frame(member=colnames(geno2), X.gene, stringsAsFactors = F)
    colnames(X.gene)[2] <- super.blk
    if(is.null(X.tmp)){
      X.tmp <- X.gene
    }else{
      X.tmp <- merge(X.tmp, X.gene, by="member")
    }
  }
  
  if(is.null(X)){
    X <- X.tmp
  }else{
    X<- merge(X, X.tmp, by="member")
  }
}


if(any(!X$member %in% pheno2[, "member"])){
  X$member = as.integer(sapply(X$member, function(x) strsplit(x, "X")[[1]][2]))
}

data <- merge(pheno2, X, by="member")
genotmp <- data[, -seq(dim(pheno2)[2]), drop = FALSE]
colnames(genotmp) <- colnames(X)[-1]

## Marginal Logistc Regression for Association Test
marginal.check <- function(geno.idx) {
  form <- as.formula(paste("affected ~", paste0(var.list, collapse="+"), "+ genotmp[, geno.idx]"))
  model <- glm(form, data = data, family = binomial(link = "logit"))
  return(model)
}

var.list <- c("age", "gender") ## covariates
check.res <- lapply(seq(dim(genotmp)[2]), marginal.check)
check.res <- lapply(check.res, summary)
check.res.otus <- sapply(check.res, function(x){
  if(x$aliased[length(x$aliased)]) return(matrix(1, ncol = 1, nrow = 4))
  else return(tail(x$coefficients, n=1))
}) 
colnames(check.res.otus) <- colnames(genotmp)
rownames(check.res.otus) <- colnames(check.res[[1]]$coefficients)
check.res.otus = as.data.frame(t(check.res.otus))

super.var.dis <- check.res.otus

# save all results in a list
res = list(super.var.form = super.var.form, # results on discovery 1
           super.var.dis = super.var.dis,  # results on discovery 2
           otus = super.variant.cand, # selected OTUs
           X = X, # Super-taxon 
           Z = Z.local, # Ranking for OTUs
           geno.transed = geno.transed) # Potential cut-off values
save(res, file = paste0(method, "_results_", blk, ".rda"))

