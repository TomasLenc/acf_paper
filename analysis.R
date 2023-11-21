rm(list=ls())

# get general parameters
source('config.R')

data_dirs <- c(
    # 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.2'
    # 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4'
    'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'
    # 'maxlag-11lags_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'
    # 'maxlag-4.8_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'
    # 'maxlag-halfTrial_meterRel-0.4_meterUnrel-0.6_1.0_1.4'
)

i_data_dir <- 1

for (i_data_dir in c(1:length(data_dirs))) {
    
    data_path = file.path(experiment_path, 'data', data_dirs[i_data_dir])
    fig_path = file.path(experiment_path, 'figures', data_dirs[i_data_dir])
    
    if (! dir.exists(fig_path)) {
        dir.create(fig_path)
    }

    results_fname <- sprintf('report_main')
    rmarkdown::render(sprintf('report_main.Rmd'),
                      output_file=sprintf('%s/%s.html', data_path, results_fname))
    

    
}

# run selection-independent reports
results_fname <- sprintf('report_sel_comparisons')
rmarkdown::render(sprintf('report_sel_comparisons.Rmd'),
                  output_file=sprintf('%s/data/%s.html', experiment_path, results_fname))
