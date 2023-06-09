
library(gghalves)




col_name <- 'z_meter_acf_subtr'
col_name_orig <- 'z_meter_acf_orig'   


df$x <- as.numeric(df$snr)
df$x_orig_start <- df$x-0.3
df$x_orig_end <- df$x+0.3

df$y <- df[, col_name] - df[, col_name_orig]

ggplot(df, aes(x, y)) + 
    # geom_segment(aes(x=x_orig_start, xend=x_orig_end,
    #                  y=.data[[col_name_orig]], yend=.data[[col_name_orig]]),
    #              size=3, alpha=0.5, color='grey70') +
    geom_hline(yintercept=0, linewidth=3, color='red', alpha=0.4) + 
    # geom_point(alpha=0.3, position=position_jitter(height=0, width=0.2)) +
    geom_half_violin(aes(group=snr), side='l', scale='width',
                     fill='grey60', color=NA, alpha=0.6, trim=T) + 
    geom_half_point(aes(group=snr), side='r', alpha=0.1, color='black') + 
    ylab(paste(col_name, '(dist from true)')) + 
    facet_wrap(~ max_lag, nrow=2) + 
    theme_cowplot()


df_bias_var <- df %>% 
    group_by(max_lag, snr) %>% 
    summarise(bias=mean(y), 
              bias_ci=sd(y) / sqrt(n()) * qnorm(1 - 0.05), 
              variance=var(y))

pos <- position_dodge(0.1)
p1 <- ggplot(df_bias_var, aes(snr, bias, group=max_lag, color=max_lag)) + 
    geom_point(position=pos) + 
    geom_errorbar(aes(ymin=bias-bias_ci, ymax=bias+bias_ci), position=pos, linewidth=1, width=0) + 
    geom_line(position=pos) + 
    geom_hline(yintercept=0, linewidth=3, color='red', alpha=0.4) + 
    theme_cowplot() + 
    theme(
        legend.position = 'none'
    )

p2 <- ggplot(df_bias_var, aes(snr, variance, group=max_lag, color=max_lag)) + 
    geom_point() + 
    geom_line() + 
    geom_hline(yintercept=0, linewidth=3, color='red', alpha=0.4) + 
    theme_cowplot()

library(ggpubr)
ggarrange(p1, p2, common.legend=TRUE, legend="right")













fname <- '/datadisk/projects_backed_up/autocorrelation/data/maxlag-2.4_meterUnrel-0.2/irType-erp2_nrep-500_snrVsNlagsNfrex.csv'

df <- read.csv(fname)

# df <- df %>% select(-X)

# df$max_freq <- ifelse(df$max_lag == 2.4, 5, 30)
# df <- df %>% relocate(max_freq, .after=max_lag)

write.csv(df, fname, row.names=F)











fname <- '/datadisk/projects_backed_up/autocorrelation/data/maxfreq-30_excl5/irType-erp2_nrep-50_emphVsNoise.csv'
df <- read.csv(fname)

df <- filter(df, z_snr > 0)

ggplot(df, aes(log10(z_snr), z_meter_fft_raw)) + 
    # geom_point(alpha=0.1, color='#701d85') + 
    geom_point(aes(color=as.factor(pat)), alpha=0.1) + 
    geom_hline(yintercept = 0, color='#5e5e5e') + 
    # coord_cartesian(xlim=c(0,100)) + 
    scale_x_reverse() +
    theme_cowplot()




















