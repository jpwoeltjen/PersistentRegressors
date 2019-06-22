
name = 'df_gls95'
dir_name = paste('/Users/jan/Desktop/PersistentRegressors/ci/',name, '.csv', sep='')
df_gls95 = read.csv(dir_name, header =TRUE, sep=",", dec="." , row.names = 1, check.names=TRUE)

usethis::use_data(df_gls95, overwrite = TRUE)

name = 'df_gls90'
dir_name = paste('/Users/jan/Desktop/PersistentRegressors/ci/',name, '.csv', sep='')
df_gls90 = read.csv(dir_name, header =TRUE, sep=",", dec="." , row.names = 1, check.names=TRUE)

usethis::use_data(df_gls90, overwrite = TRUE)

name = 'df_gls80'
dir_name = paste('/Users/jan/Desktop/PersistentRegressors/ci/',name, '.csv', sep='')
df_gls80 = read.csv(dir_name, header =TRUE, sep=",", dec="." , row.names = 1, check.names=TRUE)

usethis::use_data(df_gls80, overwrite = TRUE)

