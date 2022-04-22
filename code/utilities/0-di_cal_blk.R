di_cal <- function(pheno, geno, num.trees = 1000, num.var.boot = NULL){
  # this function calculates depth importance score of SNPs in a given SNP set
  # it constructs a random forest with a given number of trees
  # each tree is built with a randomly selected collection of SNPs w/o replacement from the SNP set
  # if the SNP used to construct tree has missing value, it is imputed based on MAF of the SNP
  
  # The inputs of the function are:
  # pheno: the first column gives eid of samples, and "affected" column indicates case-control status,
  #        and the rest columns are control variables
  # geno: SNPs data in the format of snpMatrix, each column corresponding to a single SNP,
  #        and the column name is the SNP's rsID, row name is the eid of subjects
  # num.trees: number of trees constructed in the random forest, by default it is 1000
  # num.var.boot: number of variables considered in each tree construction,
  #               by default, it is taken as the max of 30 and 1/60 of the SNP set
  
  if(dim(pheno)[1]!=dim(geno)[1]) stop("Error: nrow(geno) != nrow(pheno)\n")
  
  num.otus <- dim(geno)[2]
  all.var <- colnames(geno)
  num.trees <- num.trees
  if(is.null(num.var.boot)){
    num.var.boot <- max(30, floor(num.otus / 70)) #max(floor(num.snps / 70), 30) # 
  } else {
    num.var.boot <- num.var.boot
  }  
  
  wrapper <- function(x){
    library(diTARV)
    
    pheno.lb <- pheno[, c("member", "affected")]
    
    var.idx <- sample(seq(num.otus), num.var.boot, replace = F) # snps to include in a bootstrap sample
    geno.boot <- geno[, var.idx, drop = FALSE]
    otu.idx <- colnames(geno.boot)
    
    
    otusum.boot <- apply(geno.boot, 2, mean)
    
    if(class(pheno.lb[, "member"]) == "integer"){
      geno.boot <- data.frame(as.integer(rownames(geno.boot)), geno.boot, stringsAsFactors=FALSE)
    }else{
      geno.boot <- data.frame(rownames(geno.boot), geno.boot, stringsAsFactors=FALSE)
    }
    
    colnames(geno.boot) <- c("member", otu.idx)
    rownames(geno.boot) <- NULL
    
    data.boot <- merge(pheno.lb, geno.boot, by="member")
    data.boot <- data.boot[, -1]
    
    
    tree.t <- tryCatch(
      {tree_build(data.boot, method = "entropy", is_prune = F) # build tree w/o pruning
      },
      error=function(cond){
        NA
      },
      warning=function(cond){
        NA
      }
    )
    if(!is.na(tree.t)){
      chi.t <- tree.t$chi
      depth.t <- depth_find(tree.t)
      spvl.t <- tree.t$spvl
      di.score.node <- 2^(-depth.t) * chi.t # score for each node
      di.score.node <- di.score.node[!is.na(spvl.t)]
      depth.t <- depth.t[!is.na(spvl.t)]
      
      split.var.t <- split_var_extract(tree.t)
      res.t <- data.frame(split.var.t, di.score.node, stringsAsFactors=FALSE)
      res.t <- res.t[depth.t<=10, ]
      
      split.var.uniq <- unique(res.t[,1])
      di.score.tree <- sapply(split.var.uniq, function(x){ idx<-which(res.t[,1]==x); sum(res.t[idx, 2]) })
      res.tree <- data.frame(split.var.uniq, di.score.tree, stringsAsFactors=FALSE)
      #var.control <- which(split.var.uniq %in% colnames(pheno.lb)[-1])
      #if(length(var.control)>0){ res.tree <- res.tree[-var.control, ] }   
      
      return(res.tree)
    }
    
    return(NULL)
  }
  
  library(snowfall)
  sfInit(parallel = T, cpus = 4)
  # export variables used in the tree construction so that different threads can work on them
  sfExport("num.otus", "num.var.boot", "geno", "pheno")
  index <- seq(num.trees)
  
  start.time <- Sys.time()
  res.test <- sfLapply(index, wrapper)
  print(Sys.time() - start.time)
  
  sfStop()
  
  res.test <- res.test[!sapply(res.test, is.null)]
  
  di.all <- rep(0, num.otus)
  num.used <- rep(0, num.otus)
  
  for(t in seq_along(res.test)){
    idx <- sapply(res.test[[t]][, 1], function(x){which(all.var==x)})
    if(!is.list(idx)){
      di.all[idx] <- di.all[idx] + res.test[[t]][, 2]
      num.used[idx] <- num.used[idx] + 1 # count number of each OTUs used in the forest
    }
  }
  
  di.all.ave <- di.all / num.used
  di.all.ave[is.nan(di.all.ave)] <- 0 # 0/0 is NaN, replace them by 0
  return(cbind(di.all, di.all.ave))
}
