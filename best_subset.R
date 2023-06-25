library(tidyverse)
library(sets)

explanatory_vars <- names(datos)[1:11]

2^(length(explanatory_vars))

set.seed(123)
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

all_log_reg_models_sensitivity <- NA
all_log_reg_models_specificity <- NA
len <- NA

steps <- c(1:length(power_set))

for(i in steps){
  
  len[i] = length(power_set[[i]])
  
  if(len[i] == 1){
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
  else if(len[i] >= 2){
    
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
  'Specificity' = all_log_reg_models_specificity,
  'Num_Variables' = as_factor(len)
)

best_specificity <- which.max(all_log_reg_models_specificity)

best_subset_raw_data[best_specificity, ]

power_set[[best_specificity]]


best_subset_raw_data %>%
  filter(is.na(Specificity) == F) %>%
  group_by(Num_Variables) %>%
  summarize(Max_Spec = max(Specificity),
            Mean_Spec = mean(Specificity),
            SD_Spec = sd(Specificity),
            Max_Sens = max(Sensitivity)) %>%
  unique()


best_subset_raw_data %>%
  ggplot(aes(x = Num_Variables, y = Specificity, group = Num_Variables)) +
  geom_boxplot(aes(fill = Num_Variables)) +
  theme(legend.position = 'bottom')
