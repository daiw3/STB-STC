super_var_form_mean <- function(pheno, geno_ordered, member, file_prefix){
  num.var <- dim(geno_ordered)[2]
  p.min <- 1
  
  for(cut.off in seq(num.var)){
    s = apply(geno_ordered, 1, function(x) sum(x[1:cut.off]))
    geno.trans = rep(0, dim(geno_ordered)[1])
    geno.trans[s!=0] = apply(geno_ordered[s!=0, , drop = FALSE], 1, function(x) sum(x[1:cut.off])/sum(x[1:cut.off] != 0))
    geno.cut <- data.frame(member = member, geno.trans)
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
  
  colnames(mag.effect) <- colnames(summary(fit)$coefficient)
  super.var <- paste0(file_prefix, c("_mean"), c("+"))
  res.all <- data.frame(super.var = super.var, mag.effect, stringsAsFactors = F)
  res.all$z.value <- cut.best
  colnames(res.all) <- c("super.var", "est", "se","cut.off", "P")
  
  return(res.all)
}
