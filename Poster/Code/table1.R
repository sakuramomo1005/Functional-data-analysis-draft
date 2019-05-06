# table 0505
# draw the table
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week13')
load('data_trans.RData')
data = data_trans
placebo = data[data$group == 1, ]
drug = data[data$group == 2, ]
n_drug = dim(drug)[1]; n_placebo = dim(placebo)[1]

drug_response = round(table(drug[drug$responder == 1, ]$cluster)/ n_drug* 100)
drug_noresponse = round(table(drug[drug$responder == 0, ]$cluster)/ n_drug* 100)
placebo_response = round(table(placebo[placebo$responder == 1, ]$cluster)/ n_placebo* 100)
placebo_noresponse = round(table(placebo[placebo$responder == 0, ]$cluster)/ n_placebo* 100)

drug_response = (table(drug[drug$responder == 1, ]$cluster)/ n_drug* 100)
drug_noresponse = (table(drug[drug$responder == 0, ]$cluster)/ n_drug* 100)
placebo_response =  (table(placebo[placebo$responder == 1, ]$cluster)/ n_placebo* 100)
placebo_noresponse =  (table(placebo[placebo$responder == 0, ]$cluster)/ n_placebo* 100)

drug_table = merge(data.frame(drug_response), 
                   data.frame(drug_noresponse), by = 'Var1', all = TRUE)
drug_table[is.na(drug_table)] = 0
drug_table$Total = apply(drug_table[,2:3],1,sum)

placebo_table = merge(data.frame(placebo_response), 
                      data.frame(placebo_noresponse), by = 'Var1', all = TRUE)
placebo_table[is.na(placebo_table)] = 0
placebo_table$Total = apply(placebo_table[,2:3],1,sum)
sum(drug_table$Total)
sum(placebo_table$Total)

table1 = cbind(round(drug_table[,2:4]),
round(placebo_table[,2:4]))
table1[2,5] = 8; table1[3,5] = 32
table1 = cbind(rownames(table1), table1)
table1 = data.frame(table1)
colnames(table1) = c('Cluster','% Responders', '% Non- Responders','Total', 
                     '% Responders', '% Non- Responders','Total')
table1$Cluster = as.numeric(table1$Cluster)
table1 = rbind(table1, c('Overall', apply(table1[,2:7],2,sum)))

save(table1, file = 'table1.RData')
# Fluoxetine, n = 196
# Placebo, n = 162