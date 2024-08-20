library(survival)
library(dplyr)

comp <- function(model, data, surv_expr_input) {
  coeff <- data %>% filter(Model == model) %>% select(ENSEMBL, Coef)
  
  sur <- surv_expr_input %>% select(sample, time, status)
  expr <- surv_expr_input %>% select(-c(sample, time, status))
  
  inter <- intersect(names(expr), coeff$ENSEMBL)
  expr <- expr[, inter]
  coeff <- coeff[match(inter, coeff$ENSEMBL), ]
  
  expr <- expr %>% mutate_all(~replace(., is.nan(.), 0))
  
  surv_expr <- cbind(sur, expr)
  surv_expr_train <- surv_expr %>% select(time, status, coeff$ENSEMBL)
  
  tmp_expr <- surv_expr_train %>% select(coeff$ENSEMBL)
  score <- rowSums(tmp_expr * t(coeff$Coef))
  
  temp_surv_df <- cbind(surv_expr_train[1:2], score = score)
  temp_surv_df <- data.frame(sample = rownames(temp_surv_df), temp_surv_df)
  
  fit <- coxph(Surv(time, status) ~ score, data = temp_surv_df)
  sum.surv <- summary(fit)
  c_index <- sum.surv$concordance[1]
  
  res <- data.frame(model = model, c_index = c_index)
  return(res)
}

all <- pbapply::pblapply(unique(data$Model), FUN = function(i) {try(comp(i, data, surv_expr_input), TRUE)})
results <- do.call(rbind, all)

