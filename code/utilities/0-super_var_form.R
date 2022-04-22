super_var_form <- function(pheno, geno.transed, member){
  num.var <- dim(geno.transed)[2]
  var.list <- colnames(geno.transed)
  
  p.min.all <- rep(NA, num.var)
  cut.best.all <- rep(NA, num.var)
  mag.effect.all <- matrix(NA, num.var, 4)
  for(var.idx in seq(num.var)){
    geno <- geno.transed[, var.idx]
    cut.off.cand <- unique(geno)
    if(length(cut.off.cand) == 1){ next }
    cut.off.max <- max(cut.off.cand)
    cut.off.cand <- cut.off.cand[-which(cut.off.cand==cut.off.max)]
    
    p.min <- 1
    for(cut.off in cut.off.cand){
      geno.cut <- data.frame(member = member, as.numeric(geno <= cut.off))
      colnames(geno.cut)[1] <- "member"
      data.tmp <- merge(pheno, geno.cut, by = "member")
      data.tmp <- data.tmp[, -1]
      fit <- glm(affected ~ ., data = data.tmp, family = "binomial")
      p.cut <- tail(summary(fit)$coefficient, n = 1)[4]
      if(p.cut <= p.min){
        p.min <- p.cut
        cut.best <- cut.off
        mag.effect <- tail(summary(fit)$coefficient, n = 1)
      }
      # mag.tmp <- rbind(mag.tmp, tail(summary(fit)$coefficient, n = 1))
    }
    
    p.min.all[var.idx] <- p.min
    cut.best.all[var.idx] <- cut.best
    mag.effect.all[var.idx, ] <- mag.effect
  }
  
  colnames(mag.effect.all) <- colnames(summary(fit)$coefficient)
  res.all <- data.frame(var.list=var.list, mag.effect.all, stringsAsFactors = F)
  res.all$z.value <- cut.best.all
  colnames(res.all) <- c("super.var", "est", "se", "cut.off", "P")
  
  return(res.all)
}