
name = 'pretest_ci'
dir_name = paste('/Users/jan/Desktop/PersistentRegressors/ci/',name, '.csv', sep='')
df = read.csv(dir_name, header =FALSE, sep=",", dec=".", check.names=TRUE)

deltas = c(df[,1], df[,4])

pretest = data.frame( c(as.vector(df[,2]), as.vector(df[,5])), c(as.vector(df[,3]), as.vector(df[,6])), row.names = deltas)

colnames(pretest)= c('cl','cu')

usethis::use_data(pretest, overwrite = TRUE)
