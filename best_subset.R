library(tidyverse)
library(sets)

explanatory_vars <- c("PESO_kg","ROTA_CERV_IZQ","ROTA_CERV_DER","FLEX_FRONT_CERV_FLEXION",
                      "FLEX_FRONT_CERV_EXTENSION","FLEX_LAT_CERV_DERECHA",
                      "FLEX_LAT_CERV_IZQUIERDA","SCHOBERT_MOD", "FLEX_LUMBAR_LAT_DERECHA",
                      "FLEX_LUMBAR_LAT_IZQUIERDA", "DISTANCIA_INTERMALEOLAR")
2^(length(explanatory_vars))



power_set <- unlist(lapply(1:length(explanatory_vars),
                           # Get all combinations
                           combinat::combn,
                           x = explanatory_vars,
                           simplify = FALSE),
                    recursive = FALSE)
power_set     

paste(power_set[[100]], collapse = " + ")
formula1 <- paste('CONTROLADO', paste(power_set[[100]], collapse = " + "), sep = ' ~ ')
as.formula(formula1)


datos %>%
  group_by(CONTROLADO) %>%
  count()

confusionMatrix(as_factor(rep('SI', len = nrow(datos))), datos$CONTROLADO)

all_log_reg_models <- list()
all_log_reg_models_sensitivity <- NA
all_log_reg_models_specificity <- NA

steps <- c(1:length(power_set))

for(i in steps){
  
  len = length(power_set[[i]])
  
  if(len == 1){
    formula1 <- as.formula(paste('CONTROLADO', power_set[[i]], sep = ' ~ '))
    print(formula1)
    
    log_reg_model <- glm(formula1, data = datos, family = "binomial")
    log_reg_model_pred <- predict(log_reg_model, datos)
    log_reg_model_pred <- as_factor(if_else(log_reg_model_pred < 0.5, 'SI', 'NO'))
    
    Conf_Mat <- confusionMatrix(log_reg_model_pred, datos$CONTROLADO)
    all_log_reg_models_sensitivity[i] = Conf_Mat[['byClass']][['Sensitivity']]
    all_log_reg_models_specificity[i] = Conf_Mat[['byClass']][['Specificity']]
    message = paste0('Specificity: ', all_log_reg_models_specificity[i])
    print(message)
    
  }
  else if(len > 2){
    
    formula1 <- as.formula(paste('CONTROLADO', paste(power_set[[i]], collapse = ' + '), sep = ' ~ '))
    print(formula1)
    
    log_reg_model <- glm(formula1, data = datos, family = "binomial")
    log_reg_model_pred <- predict(log_reg_model, datos)
    log_reg_model_pred <- as_factor(if_else(log_reg_model_pred < 0.5, 'SI', 'NO'))
    
    Conf_Mat <- confusionMatrix(log_reg_model_pred, datos$CONTROLADO)
    all_log_reg_models_sensitivity[i] = Conf_Mat[['byClass']][['Sensitivity']]
    all_log_reg_models_specificity[i] = Conf_Mat[['byClass']][['Specificity']]
    message = paste0('Specificity: ', all_log_reg_models_specificity[i])
    print(message)
    
  }
}

hist(all_log_reg_models_specificity)

best_subset_raw_data <- data.frame(
  'Sensitivity' = all_log_reg_models_sensitivity,
  'Specificity' = all_log_reg_models_specificity
)

best_specificity <- which.max(all_log_reg_models_specificity)

best_subset_raw_data[best_specificity, ]

power_set[[best_specificity]]
