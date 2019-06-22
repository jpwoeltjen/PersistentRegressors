## code to prepare `DATASET` dataset goes here
table_deltas = c(-0.999, -0.975, -0.95, -0.925,
                 -0.9, -0.875, -0.85, -0.825,
                 -0.8, -0.775, -0.75, -0.725,
                 -0.7, -0.675, -0.65, -0.625,
                 -0.6, -0.575, -0.55, -0.525,
                 -0.5, -0.475, -0.45, -0.425,
                 -0.4, -0.375, -0.35, -0.325,
                 -0.3, -0.275, -0.25, -0.225,
                 -0.2, -0.175, -0.15, -0.125,
                 -0.1, -0.075, -0.05, -0.025)





for (index in 1:length(table_deltas)){
  d = table_deltas[index]
  name = paste('df_gls_delta',d, sep='')
  dir_name = paste('/Users/jan/Desktop/PersistentRegressors/ci/',name, '.csv', sep='')
  df = read.csv(dir_name, header =TRUE, sep=",", dec="." , row.names = 1, check.names=TRUE)
  if (index==1){
    bonferroniQ = data.frame(row.names = rownames(df))

  }
  list = NULL
  for (i in 1:dim(df)[1]){
    list[[i]] = c(df[i,1],df[i,2])
  }
  bonferroniQ[1:61,index] = list(list)


}
  colnames(bonferroniQ)=table_deltas

usethis::use_data(bonferroniQ,overwrite = TRUE)
