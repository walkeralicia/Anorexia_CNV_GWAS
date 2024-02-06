
# Function to calculate metrics
calculate_burden_function <- function(cnv_data, pheno, predictor_col, cnv_type) {
  
  avg <- pheno %>% group_by(AN) %>% summarise_at(vars(predictor_col), list(mean = mean)) %>% as.data.frame()
  avg_case <- signif(avg[avg$AN==2, c("mean")],3)
  avg_cont <- signif(avg[avg$AN==1, c("mean")],3)
  
  pheno$AN <- as.factor(pheno$AN)
  m <- glm(as.formula(paste0("AN ~ Sex + Age + Array + ", predictor_col)), data = pheno, family = "binomial")
  or <- exp(cbind(OR = coef(m), confint(m))) %>% as.data.frame()
  
  OR <- signif(or[c(predictor_col), c("OR")], 3)
  dp <- nchar(OR) - 2
  lower <- round(or[c(predictor_col), c("2.5 %")], dp)
  upper <- round(or[c(predictor_col), c("97.5 %")], dp)
  CI <- paste0("(", lower, "-", upper, ")")
  pval <- signif(summary(m)$coefficients[c(predictor_col), c("Pr(>|z|)")], 3)
  
  vec <- c(cnv_type, 
           predictor_col, 
           OR, 
           lower, 
           upper, 
           CI,
           pval,
           avg_case,
           avg_cont)
  names(vec) <- c("Type", "Test", "OR", "lower", "upper", "CI", "pval", "Avg_case", "Avg_cont")
  return(vec)
}


