rm(list=ls())

# get general parameters
source('lib/R/config.R')

data_dirs <- c(
    # 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-false_keepBand-false',
    'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-false'
    # 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-true'
)

i_data_dir <- 1

for (i_data_dir in c(1:length(data_dirs))) {
    
    data_path = file.path(experiment_path, 'results', data_dirs[i_data_dir])
    fig_path = file.path(experiment_path, 'figures', data_dirs[i_data_dir])
    
    if (! dir.exists(fig_path)) {
        dir.create(fig_path)
    }

    results_fname <- sprintf('report_main')
    rmarkdown::render(sprintf('lib/R/report_main.Rmd'),
                      output_file=sprintf('%s/%s.html', data_path, results_fname))
    

    
}
