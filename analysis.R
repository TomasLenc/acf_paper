rm(list=ls())

# get general parameters
source('config.R')

# make sure the directory with reports is empty
reports_path <- file.path(experiment_path, 'reports')
if (!dir.exists(reports_path)) dir.create(reports_path)

# Only noise
results_fname <- sprintf('date-%s_time-%s_report_only_noise', 
                         Sys.Date(), format(Sys.time(),'%H-%M'))

rmarkdown::render(sprintf('report_only_noise.Rmd'), 
                  output_file=sprintf('%s/%s.html', reports_path, results_fname))

# Syncrange
results_fname <- sprintf('date-%s_time-%s_report_syncrange', 
                         Sys.Date(), format(Sys.time(),'%H-%M'))

rmarkdown::render(sprintf('report_syncrange.Rmd'), 
                  output_file=sprintf('%s/%s.html', reports_path, results_fname))

