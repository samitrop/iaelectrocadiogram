require(ggplot2)
require(precrec)
require(epiR)
require(ROCR)
require(survey)
require(ggpubr)
require(pROC)

################################################################################
### Data reading and pre-processing
################################################################################

data_SMS <- read.csv2(file = 'final_data.csv', fileEncoding = 'CP1252')

data_SMS$Q1 <- ifelse(data_SMS$Q1 == 'Yes', 1, 0)
data_SMS$Q2 <- ifelse(data_SMS$Q2 == 'Yes', 1, 0)
data_SMS$Q3 <- ifelse(data_SMS$Q3 == 'Yes', 1, 0)

### Final models
Prob_AI_ECG_EPI <- exp(-5.4204 + 1.4210 * data_SMS$Q1 + 1.1824 * data_SMS$Q2 + 1.2557 * data_SMS$Q3 + 1.5370 * data_SMS$AI_ECG)/(exp(-5.4204 + 1.4210 * data_SMS$Q1 + 1.1824 * data_SMS$Q2 + 1.2557 * data_SMS$Q3 + 1.5370 * data_SMS$AI_ECG) + 1)
Prob_EPI_only <- exp(-4.7025 + 1.4836 * data_SMS$Q1 + 1.3627 * data_SMS$Q2 + 1.2050 * data_SMS$Q3)/(exp(-4.7025 + 1.4836 * data_SMS$Q1 + 1.3627 * data_SMS$Q2 + 1.2050 * data_SMS$Q3) + 1)
Prob_AI_ECG_only <- exp(-3.614 + 2.678 * data_SMS$AI_ECG)/(exp(-3.614 + 2.678 * data_SMS$AI_ECG) + 1)

### Prediction
pred_AI_ECG_EPI <- factor(ifelse(Prob_AI_ECG_EPI < 0.135, 'Negative', 'Positive'))
pred_EPI_only <- factor(ifelse(Prob_EPI_only < 0.2389333, 'Negative', 'Positive'))
pred_AI_ECG_only <- factor(ifelse(Prob_AI_ECG_only < 0.08445402, 'Negative', 'Positive'))

f1 <- function(x) format(round(x,3), nsmall = 3)


################################################################################
### Prediction metrics calculation
################################################################################

metrics_function <- function(
        n_positive_alerts,
        n_negative_alerts,
        n_positive_alert_positive_soro,
        n_positive_alert_negative_soro,
        n_negative_alert_negative_soro,
        n_negative_alert_positive_soro)
{
    result_matrix <- matrix(NA, nrow = 1000, ncol = 9)

    p_positive_alerts <- n_positive_alerts/(n_positive_alerts + n_negative_alerts)
    
    ### Positive predictive value
    P_positive_soro_positive_alert <- n_positive_alert_positive_soro/(n_positive_alert_positive_soro+n_positive_alert_negative_soro) 
    P_negative_soro_positive_alert <- 1 - P_positive_soro_positive_alert 
    
    ### Negative predictive value
    P_negative_soro_negative_alert <- n_negative_alert_negative_soro/(n_negative_alert_positive_soro + n_negative_alert_negative_soro) 
    P_positive_soro_negative_alert <- 1 - P_negative_soro_negative_alert 
    
    P_positive <- n_positive_alerts/(n_positive_alerts + n_negative_alerts)
    P_negativo <- 1-P_positive
    
    ### Prevalence
    P_positive_soro <- P_positive_soro_positive_alert * P_positive + P_positive_soro_negative_alert * P_negativo
    
    ### Odds ratio
    odds_ratio <- (P_positive_soro_positive_alert/P_negative_soro_positive_alert)/(P_positive_soro_negative_alert/P_negative_soro_negative_alert)
    
    ### Sensitivity
    P_positive_alert_positive_soro <- (P_positive_soro_positive_alert * P_positive)/(P_positive_soro_positive_alert * P_positive + P_positive_soro_negative_alert * P_negativo)
    
    ### Specificity
    P_negative_alert_negative_soro <- (P_negative_soro_negative_alert * P_negativo)/(P_negative_soro_negative_alert * P_negativo + P_negative_soro_positive_alert * P_positive)
    
    ### Positive likelihood ratio
    Positive_likelihood_ratio <- P_positive_alert_positive_soro/(1 - P_negative_alert_negative_soro)
    
    ### Negative likelihood ratio
    Negative_likelihood_ratio <- (1 - P_positive_alert_positive_soro)/P_negative_alert_negative_soro
    
    ### F1 value
    f1_score <- (2*n_positive_alert_positive_soro)/(2*n_positive_alert_positive_soro + n_positive_alert_negative_soro + n_negative_alert_positive_soro)
    
    Estatística <- f1(c(P_positive_soro,
                        P_positive_alert_positive_soro,
                        P_negative_alert_negative_soro,
                        P_positive_soro_positive_alert,
                        P_negative_soro_negative_alert,
                        Positive_likelihood_ratio,
                        Negative_likelihood_ratio,
                        odds_ratio,
                        f1_score))
    
    ### Bootstrap simulation for confidence intervals
    for(i in 1:1000){
        n_positive_alerts_simul <- rbinom(1, size = n_positive_alerts + n_negative_alerts, prob = p_positive_alerts)
        n_negative_alerts_simul <- n_positive_alerts + n_negative_alerts - n_positive_alerts_simul
        p_positive_alerts_simul  <- n_positive_alerts_simul/(n_positive_alerts + n_negative_alerts)
        
        n_positive_alert_positive_soro_simul <-rbinom(1, size = n_positive_alert_positive_soro + n_positive_alert_negative_soro, 
                                                      prob = n_positive_alert_positive_soro/(n_positive_alert_positive_soro + n_positive_alert_negative_soro)) 
        n_positive_alert_negative_soro_simul <- n_positive_alert_positive_soro + n_positive_alert_negative_soro - n_positive_alert_positive_soro_simul
        
        n_negative_alert_positive_soro_simul <-rbinom(1, size = n_negative_alert_positive_soro + n_negative_alert_negative_soro, 
                                                      prob = n_negative_alert_positive_soro/(n_negative_alert_positive_soro + n_negative_alert_negative_soro)) 
        n_negative_alert_negative_soro_simul <- n_negative_alert_positive_soro + n_negative_alert_negative_soro - n_negative_alert_positive_soro_simul
        
        P_positive_soro_positive_alert_simul <- n_positive_alert_positive_soro_simul/(n_positive_alert_positive_soro_simul+n_positive_alert_negative_soro_simul) ### Positivo
        P_negative_soro_positive_alert_simul <- 1 - P_positive_soro_positive_alert_simul ### Falso positivo
        P_negative_soro_negative_alert_simul <- n_negative_alert_negative_soro_simul/(n_negative_alert_positive_soro_simul + n_negative_alert_negative_soro_simul) ### Negativo
        P_positive_soro_negative_alert_simul <- 1 - P_negative_soro_negative_alert_simul ### Falso negativo
        P_positive_simul <- n_positive_alerts_simul/(n_positive_alerts_simul + n_negative_alerts_simul)
        P_negativo_simul <- 1-P_positive_simul
        P_positive_soro_simul <- P_positive_soro_positive_alert_simul * P_positive_simul + P_positive_soro_negative_alert_simul * P_negativo_simul
        
        OR_simul <- (P_positive_soro_positive_alert_simul/P_negative_soro_positive_alert_simul)/(P_positive_soro_negative_alert_simul/P_negative_soro_negative_alert_simul)
        
        P_positive_alert_positive_soro_simul <- (P_positive_soro_positive_alert_simul * P_positive_simul)/(P_positive_soro_positive_alert_simul * P_positive_simul + P_positive_soro_negative_alert_simul * P_negativo_simul)
        
        P_negative_alert_negative_soro_simul <- (P_negative_soro_negative_alert_simul * P_negativo_simul)/(P_negative_soro_negative_alert_simul * P_negativo_simul + P_negative_soro_positive_alert_simul * P_positive_simul)
        
        Positive_likelihood_ratio_simul <- P_positive_alert_positive_soro_simul/(1 - P_negative_alert_negative_soro_simul)
        Negative_likelihood_ratio_simul <- (1 - P_positive_alert_positive_soro_simul)/P_negative_alert_negative_soro_simul
        
        f1_score_simul <- (2*n_positive_alert_positive_soro_simul)/(2*n_positive_alert_positive_soro_simul + n_positive_alert_negative_soro_simul + n_negative_alert_positive_soro_simul)
        
        result_matrix[i,] <- c(P_positive_soro_simul,
                               P_positive_alert_positive_soro_simul,
                               P_negative_alert_negative_soro_simul,
                               P_positive_soro_positive_alert_simul,
                               P_negative_soro_negative_alert_simul,
                               Positive_likelihood_ratio_simul,
                               Negative_likelihood_ratio_simul,
                               OR_simul,
                               f1_score_simul)
    }
    
    ### Bootstrap percentiles
    q025 <- f1(apply(result_matrix, 2, quantile, prob = 0.025))
    q975 <- f1(apply(result_matrix, 2, quantile, prob = 0.975))
    
    ### Table of results
    data_result <- data.frame(`Estimativa (IC95%)` = paste0(Estatística, ' (', q025, ' ; ', q975, ')'),
                              row.names = c('Prevalence', 'Sensitivity', 'Specificity',
                                            'Positive predictive value',
                                            'Negative predictive value',
                                            'Positive likelihood ratio',
                                            'Negative likelihood ratio',
                                            'Odds ratio',
                                            'F1 score'))   
    
    data_result
}

### General
metrics_general <- metrics_function(n_positive_alerts = 6795,
                                    n_negative_alerts = 52005,
                                    n_positive_alert_positive_soro = 899,
                                    n_positive_alert_negative_soro = 1588,
                                    n_negative_alert_negative_soro = 925,
                                    n_negative_alert_positive_soro = 92)



### Hyperendemic region
metrics_Hyper <- metrics_function(n_positive_alerts = 5712,
                                  n_negative_alerts = 25524,
                                  n_positive_alert_positive_soro = 863,
                                  n_positive_alert_negative_soro = 1220,
                                  n_negative_alert_negative_soro = 616,
                                  n_negative_alert_positive_soro = 83)


### Endemic region
metrics_Endemic <- metrics_function(n_positive_alerts = 1083,
                                    n_negative_alerts = 26481,
                                    n_positive_alert_positive_soro = 62,
                                    n_positive_alert_negative_soro = 368,
                                    n_negative_alert_negative_soro = 283,
                                    n_negative_alert_positive_soro = 9)


### Male
metrics_Male <- metrics_function(n_positive_alerts = 2567,
                                 n_negative_alerts = 21539,
                                 n_positive_alert_positive_soro = 337,
                                 n_positive_alert_negative_soro = 598,
                                 n_negative_alert_negative_soro = 321,
                                 n_negative_alert_positive_soro = 44)


### Female
metrics_Female <- metrics_function(n_positive_alerts = 4228,
                                   n_negative_alerts = 30466,
                                   n_positive_alert_positive_soro = 588,
                                   n_positive_alert_negative_soro = 990,
                                   n_negative_alert_negative_soro = 578,
                                   n_negative_alert_positive_soro = 48)


### Adults
metrics_Adults <- metrics_function(n_positive_alerts = 3224,
                                   n_negative_alerts = 31580,
                                   n_positive_alert_positive_soro = 373,
                                   n_positive_alert_negative_soro = 824,
                                   n_negative_alert_negative_soro = 563,
                                   n_negative_alert_positive_soro = 43)


### Elderly
metrics_Elderly <- metrics_function(n_positive_alerts = 3571,
                                    n_negative_alerts = 20425,
                                    n_positive_alert_positive_soro = 552,
                                    n_positive_alert_negative_soro = 763,
                                    n_negative_alert_negative_soro = 331,
                                    n_negative_alert_positive_soro = 49)


### ECG: normal
metrics_Normal <- metrics_function(n_positive_alerts = 2407,
                                   n_negative_alerts = 27267,
                                   n_positive_alert_positive_soro = 208,
                                   n_positive_alert_negative_soro = 735,
                                   n_negative_alert_negative_soro = 506,
                                   n_negative_alert_positive_soro = 41)


### ECG: minor abnormalities
metrics_Minor <- metrics_function(n_positive_alerts = 2293,
                                  n_negative_alerts = 16801,
                                  n_positive_alert_positive_soro = 279,
                                  n_positive_alert_negative_soro = 557,
                                  n_negative_alert_negative_soro = 278,
                                  n_negative_alert_positive_soro = 30)


### ECG: major abnormalities
metrics_Major <- metrics_function(n_positive_alerts = 1919,
                                  n_negative_alerts = 6364,
                                  n_positive_alert_positive_soro = 420,
                                  n_positive_alert_negative_soro = 261,
                                  n_negative_alert_negative_soro = 89,
                                  n_negative_alert_positive_soro = 21)


### Reported Chagas
metrics_Reported <- metrics_function(n_positive_alerts = 1412,
                                     n_negative_alerts = 411,
                                     n_positive_alert_positive_soro = 507,
                                     n_positive_alert_negative_soro = 58,
                                     n_negative_alert_negative_soro = 6,
                                     n_negative_alert_positive_soro = 26)


### Non-reported Chagas
metrics_non_reported <- metrics_function(n_positive_alerts = 5383,
                                         n_negative_alerts = 51594,
                                         n_positive_alert_positive_soro = 418,
                                         n_positive_alert_negative_soro = 1530,
                                         n_negative_alert_negative_soro = 893,
                                         n_negative_alert_positive_soro = 66)


### General (only EPI)
metrics_General_EPI <- metrics_function(n_positive_alerts = 5591,
                                        n_negative_alerts = 53209,
                                        n_positive_alert_positive_soro = 836,
                                        n_positive_alert_negative_soro = 1344,
                                        n_negative_alert_negative_soro = 1143,
                                        n_negative_alert_positive_soro = 181)


### Hyperendemic region (only EPI)
metrics_Hyper_EPI <- metrics_function(n_positive_alerts = 4800,
                                      n_negative_alerts = 26463,
                                      n_positive_alert_positive_soro = 788,
                                      n_positive_alert_negative_soro = 1062,
                                      n_negative_alert_negative_soro = 774,
                                      n_negative_alert_positive_soro = 158)


### Endemic region (only EPI)
metrics_Endemic_EPI <- metrics_function(n_positive_alerts = NULL,
                                        n_negative_alerts = NULL,
                                        n_positive_alert_positive_soro = NULL,
                                        n_positive_alert_negative_soro = NULL,
                                        n_negative_alert_negative_soro = NULL,
                                        n_negative_alert_positive_soro = NULL)

n_positive_alerts <- 791
n_negative_alerts <- 26773
p_casos <- n_positive_alerts/(n_positive_alerts + n_negative_alerts)

n_positivo_doente <- 48
n_positivo_sadio <- 282
n_negativo_sadio <- 369
n_negativo_doente <- 23


################################################################################
### ROC analysis (area under roc curve)
################################################################################

### Point estimation
roc_function <- function(data_set, n_positive_alerts, n_negative_alerts,
                       n_positive_alerts_EPI_only, n_negative_alerts_EPI_only,
                       n_positive_alerts_AI_ECG_only, n_negative_alerts_AI_ECG_only){
    data_roc <- data.frame(cutoff = c(seq(0.0000001, 0.60, 0.015), 0.70, 0.80, 0.90, 0.999),
                           sens_EPI_only = NA, spec_EPI_only = NA,
                           sens_EPI_AI_ECG = NA, spec_EPI_AI_ECG = NA,
                           sens_AI_ECG_only = NA, spec_AI_ECG_only = NA)
    
    Q1 <- data_set$Q1
    Q2 <- data_set$Q2
    Q3 <- data_set$Q3
    
    Prob_EPI_AI_ECG <- exp(-5.4204 + 1.4210 * Q1 + 1.1824 * Q2 + 1.2557 * Q3 + 1.5370 * AI_ECG)/(exp(-5.4204 + 1.4210 * Q1 + 1.1824 * Q2 + 1.2557 * Q3 + 1.5370 * AI_ECG) + 1)
    Prob_EPI_only <- exp(-4.7025 + 1.4836 * Q1 + 1.3627 * Q2 + 1.2050 * Q3)/(exp(-4.7025 + 1.4836 * Q1 + 1.3627 * Q2 + 1.2050 * Q3) + 1)
    Prob_AI_ECG_only <- exp(-3.614 + 2.678 * AI_ECG)/(exp(-3.614 + 2.678 * AI_ECG) + 1)
    
    data_set$Weight <- NA
    data_set$Weight[which(data_set$Alert == 'Negative')] <- n_negative_alerts/(n_negative_alerts + n_positive_alerts)
    data_set$Weight[which(data_set$Alert == 'Positive')] <- 1 - n_negative_alerts/(n_negative_alerts + n_positive_alerts)
    
    data_set$Weight_EPI_only <- NA
    data_set$Weight_EPI_only[which(data_set$Alert_EPI_only == 'Negative')] <- n_negative_alerts_EPI_only/(n_negative_alerts_EPI_only + n_positive_alerts_EPI_only)
    data_set$Weight_EPI_only[which(data_set$Alert_EPI_only == 'Positive')] <- 1 - n_negative_alerts_EPI_only/(n_negative_alerts_EPI_only + n_positive_alerts_EPI_only)
    
    data_set$Weight_AI_ECG_only <- NA
    data_set$Weight_AI_ECG_only[which(data_set$Alert_AI_ECG_only == 'Negative')] <- n_negative_alerts_AI_ECG_only/(n_negative_alerts_AI_ECG_only + n_positive_alerts_AI_ECG_only)
    data_set$Weight_AI_ECG_only[which(data_set$Alert_AI_ECG_only == 'Positive')] <- 1 - n_negative_alerts_AI_ECG_only/(n_negative_alerts_AI_ECG_only + n_positive_alerts_AI_ECG_only)
    
    for(i in 1:nrow(data_roc)){
        class_EPI_AI_ECG <- factor(c('Positive', 'Negative', ifelse(Prob_EPI_AI_ECG < data_roc$cutoff[i], 'Negative', 'Positive')))
        class_EPI_AI_ECG <- relevel(class_EPI_AI_ECG, ref = 'Negative')
        data_set$class_EPI_AI_ECG <- class_EPI_AI_ECG[-c(1,2)]
        
        design <- svydesign(ids = ~1, data = data_set, weights = ~Weight)
        
        # Tabela de contingência ponderada
        tab <- svytable(~ Serology + class_EPI_AI_ECG, design)
        
        # Calcular sensibilidade e specificidade
        sens <- tab[2,2] / sum(tab[2,])
        sp  <- tab[1,1] / sum(tab[1,])
        
        data_roc$sens_EPI_AI_ECG[i] <- sens
        data_roc$spec_EPI_AI_ECG[i] <- sp
        
        class_EPI_only <- factor(c('Positive', 'Negative', ifelse(Prob_EPI_only < data_roc$cutoff[i], 'Negative', 'Positive')))
        class_EPI_only <- relevel(class_EPI_only, ref = 'Negative')
        data_set$class_EPI_only <- class_EPI_only[-c(1,2)]
        
        design <- svydesign(ids = ~1, data = data_set, weights = ~Weight_EPI_only)
        
        # Tabela de contingência ponderada
        tab <- svytable(~ Serology + class_EPI_only, design)
        
        # Calcular sensibilidade e specificidade
        sens <- tab[2,2] / sum(tab[2,])
        sp  <- tab[1,1] / sum(tab[1,])
        
        data_roc$sens_EPI_only[i] <- sens
        data_roc$spec_EPI_only[i] <- sp
        
        class_AI_ECG_only <- factor(c('Positive', 'Negative', ifelse(Prob_AI_ECG_only < data_roc$cutoff[i], 'Negative', 'Positive')))
        class_AI_ECG_only <- relevel(class_AI_ECG_only, ref = 'Negative')
        data_set$class_AI_ECG_only <- class_AI_ECG_only[-c(1,2)]
        
        design <- svydesign(ids = ~1, data = data_set, weights = ~Weight_AI_ECG_only)
        
        # Tabela de contingência ponderada
        tab <- svytable(~ Serology + class_AI_ECG_only, design)
        
        # Calcular sensibilidade e specificidade
        sens <- tab[2,2] / sum(tab[2,])
        sp  <- tab[1,1] / sum(tab[1,])
        
        data_roc$sens_AI_ECG_only[i] <- sens
        data_roc$spec_AI_ECG_only[i] <- sp
        
    }
    area_EPI_AI_ECG <- pracma::trapz(data_roc$spec_EPI_AI_ECG, data_roc$sens_EPI_AI_ECG)
    area_EPI_only <- pracma::trapz(data_roc$spec_EPI_only, data_roc$sens_EPI_only)
    area_AI_ECG_only <- pracma::trapz(data_roc$spec_AI_ECG_only, data_roc$sens_AI_ECG_only)
    c(area_EPI_AI_ECG, area_EPI_only, area_AI_ECG_only)
}

### Bootstrap confidence intervals

boot_AUC <- function(data_set, n_positive_alerts, n_negative_alerts,
                     n_positive_alerts_EPI_only, n_negative_alerts_EPI_only,
                     n_positive_alerts_AI_ECG_only, n_negative_alerts_AI_ECG_only){
    boot_vec <- matrix(NA, nrow = 500, ncol = 3)
    for(i in 1:500){
        data_set2 <- data_set[sample(1:nrow(data_set), nrow(data_set), replace = TRUE),]
        boot_vec[i,] <- roc_function(data_set2, n_positive_alerts, n_negative_alerts,
                                   n_positive_alerts_EPI_only, n_negative_alerts_EPI_only,
                                   n_positive_alerts_AI_ECG_only, n_negative_alerts_AI_ECG_only)
        print(i)
    }
    boot_vec
}

### Calculation 
f1 <- roc_function(data_SMS, n_positive_alerts = 5712, n_negative_alerts = 25524,
                 n_positive_alerts_EPI_only = 5591, n_negative_alerts_EPI_only = 53209,
                 n_positive_alerts_AI_ECG_only = 10576, n_negative_alerts_AI_ECG_only = 48578)

b1 <- boot_AUC(data_SMS, n_positive_alerts = 5712, n_negative_alerts = 25524,
               n_positive_alerts_EPI_only = 5591, n_negative_alerts_EPI_only = 53209,
               n_positive_alerts_AI_ECG_only = 10576, n_negative_alerts_AI_ECG_only = 48578)
t2 <- Sys.time()

quantile(b1[,1], c(0.025, 0.975))
quantile(b1[,2], c(0.025, 0.975))
quantile(b1[,3], c(0.025, 0.975))


################################################################################
### ROC and PR curves
################################################################################

roc_function_para_plot <- function(data_set, n_positive_alerts, n_negative_alerts,
                                 n_positive_alerts_only_EPI, n_negative_alerts_only_IA,
                                 n_positive_alerts_only_ECG_AI, n_negative_alerts_only_ECG_AI){
    
    data_roc <- data.frame(cutoff = seq(0.0000001, 0.999, 0.01),
                           sens = NA, spec = NA, PPV = NA,
                           sens_EPI_only = NA, spec_EPI_only = NA, PPV_EPI_only = NA,
                           sens_AI_ECG_only = NA, spec_AI_ECG_only = NA, PPV_AI_ECG_only = NA)
    
    Prob_AI_ECG_EPI <- with(data_set, exp(-5.4204 + 1.4210 * Q1 + 1.1824 * Q2 + 1.2557 * Q3 + 1.5370 * AI_ECG)/(exp(-5.4204 + 1.4210 * Q1 + 1.1824 * Q2 + 1.2557 * Q3 + 1.5370 * AI_ECG) + 1))
    Prob_EPI_only <- with(data_set, exp(-4.7025 + 1.4836 * Q1 + 1.3627 * Q2 + 1.2050 * Q3)/(exp(-4.7025 + 1.4836 * Q1 + 1.3627 * Q2 + 1.2050 * Q3) + 1))
    Prob_AI_ECG_only <- with(data_set, exp(-3.614 + 2.678 * AI_ECG)/(exp(-3.614 + 2.678 * AI_ECG) + 1))
    
    data_set$Weigth <- NA
    data_set$Weigth[which(data_set$Alert == 'Negative')] <- n_negative_alerts/(n_negative_alerts + n_positive_alerts)
    data_set$Weigth[which(data_set$Alert == 'Positive')] <- 1 - n_negative_alerts/(n_negative_alerts + n_positive_alerts)
    
    for(i in 1:nrow(data_roc)){
        classification <- factor(c('Positive', 'Negative', ifelse(Prob_AI_ECG_EPI < data_roc$cutoff[i], 'Negative', 'Positive')))
        classification <- relevel(classification, ref = 'Negative')
        data_set$classification <- classification[-c(1,2)]
        
        design <- svydesign(ids = ~1, data = data_set, weights = ~Weigth)
        
        # Tabela de contingência ponderada
        tab <- svytable(~ Serology + classification, design)
        
        # Calcular sensibilidade e specificidade
        sens <- tab[2,2] / sum(tab[2,])
        spec  <- tab[1,1] / sum(tab[1,])
        PPV <- tab[2,2] / sum(tab[,2])
        
        data_roc$sens[i] <- sens
        data_roc$spec[i] <- spec
        data_roc$PPV[i] <- PPV
        
        classification_EPI_only <- factor(c('Positive', 'Negative', ifelse(Prob_EPI_only < data_roc$cutoff[i], 'Negative', 'Positive')))
        classification_EPI_only <- relevel(classification_EPI_only, ref = 'Negative')
        data_set$classification_EPI_only <- classification_EPI_only[-c(1,2)]
        
        classification_AI_ECG_only <- factor(c('Positive', 'Negative', ifelse(Prob_AI_ECG_only < data_roc$cutoff[i], 'Negative', 'Positive')))
        classification_AI_ECG_only <- relevel(classification_AI_ECG_only, ref = 'Negative')
        data_set$classification_AI_ECG_only <- classification_AI_ECG_only[-c(1,2)]
        
        
        design <- svydesign(ids = ~1, data = data_set, weights = ~Weigth)
        
        # Tabela de contingência ponderada
        tab <- svytable(~ Serology + classification_EPI_only, design)
        
        # Calcular sensibilidade e specificidade
        sens_EPI_only <- tab[2,2] / sum(tab[2,])
        spec_EPI_only  <- tab[1,1] / sum(tab[1,])
        PPV_EPI_only <- tab[2,2] / sum(tab[,2])
        
        
        data_roc$sens_EPI_only[i] <- sens_EPI_only
        data_roc$spec_EPI_only[i] <- spec_EPI_only
        data_roc$PPV_EPI_only[i] <- PPV_EPI_only
        
        # Tabela de contingência ponderada
        tab <- svytable(~ Serology + classification_AI_ECG_only, design)
        
        # Calcular sensibilidade e specificidade
        sens_AI_ECG_only <- tab[2,2] / sum(tab[2,])
        spec_AI_ECG_only  <- tab[1,1] / sum(tab[1,])
        PPV_AI_ECG_only <- tab[2,2] / sum(tab[,2])
        
        data_roc$sens_AI_ECG_only[i] <- sens_AI_ECG_only
        data_roc$spec_AI_ECG_only[i] <- spec_AI_ECG_only
        data_roc$PPV_AI_ECG_only[i] <- PPV_AI_ECG_only
        
    }
    data_roc
}


f1_general <- roc_function_para_plot(data_SMS, n_positive_alerts = 5712, n_negative_alerts = 25524,
                                   n_positive_alerts_only_EPI = 5591, n_negative_alerts_only_IA = 53209,
                                   n_positive_alerts_only_ECG_AI = 10576, n_negative_alerts_only_ECG_AI = 48578)

g1 <- ggplot(data = f1_general)+
    geom_line(aes(x = 1-spec, y = sens, col = 'ECG-AI + Risk factors'), size = 1)+
    geom_line(aes(x = 1-spec_EPI_only, y = sens_EPI_only, col = 'Risk factors'), size = 1)+
    geom_line(aes(x = 1-spec_AI_ECG_only, y = sens_AI_ECG_only, col = 'ECG-AI'), size = 1)+
    scale_color_manual(
        values = c("ECG-AI + Risk factors" = "#F8766D",
                   "Risk factors" = "#00BA38",
                   "ECG-AI" = "#619CFF"),
        breaks = c("ECG-AI + Risk factors",
                   "Risk factors",
                   "ECG-AI")
    )+
    theme_classic(base_size = 14)+
    scale_x_continuous(breaks = seq(0,1,0.1))+
    scale_y_continuous(breaks = seq(0,1,0.1))+
    theme(
        legend.position = c(0.75, 0.2),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", colour = "grey70")
    )+
    theme(legend.title=element_blank())+
    geom_abline(intercept = 0, slope = 1, size = 1, lty = 2)+
    xlab('False Positive Rate (FPR)')+
    ylab('True Positive Rate (TPR)')


g2 <- ggplot(data = f1_general)+
    geom_line(aes(x = sens, y = PPV, col = 'ECG-AI + Risk factors'), size = 1)+
    geom_line(aes(x = sens_EPI_only, y = PPV_EPI_only, col = 'Risk factors'), size = 1)+
    geom_line(aes(x = sens_AI_ECG_only, y = PPV_AI_ECG_only, col = 'ECG-AI'), size = 1)+
    scale_color_manual(
        values = c("ECG-AI + Risk factors" = "#F8766D",
                   "Risk factors" = "#00BA38",
                   "ECG-AI" = "#619CFF"),
        breaks = c("ECG-AI + Risk factors",
                   "Risk factors",
                   "ECG-AI")
    )+
    theme_classic(base_size = 14)+
    scale_x_continuous(breaks = seq(0,1,0.1))+
    scale_y_continuous(breaks = seq(0,1,0.1))+
    theme(
        legend.position = 'none'
    )+
    xlab('Recall')+
    ylab('Precision')

ggarrange(g1, g2, ncol = 2, nrow = 1)