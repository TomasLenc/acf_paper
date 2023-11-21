
my_theme <- function(){
    theme(
        axis.line.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=fontsize), 
        axis.title = element_text(size=fontsize), 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(),
        legend.text = element_text(size=fontsize), 
        legend.title = element_text(size=fontsize), 
        strip.background = element_blank(), 
        strip.text = element_text(size=fontsize, face='bold')
    )
}


col_scale <- scale_color_brewer(palette='Paired')

plot_data <- function(df, col_name, col_name_orig, col_max){
    
    df$x <- as.numeric(df$snr)
    df$x_orig_start <- df$x-0.3
    df$x_orig_end <- df$x+0.3
    
    # df$y <- df[, col_name] - df[, col_name_orig]
    
    p <- ggplot(df, aes(x, .data[[col_name]])) + 
        # geom_segment(aes(x=x_orig_start, xend=x_orig_end,
        #                  y=.data[[col_name_orig]], yend=.data[[col_name_orig]]),
        #              size=3, alpha=0.5, color='grey70') +
        geom_hline(yintercept=0, linewidth=1, color='black', alpha=1, linetype='dotted') + 
        geom_hline(aes(yintercept=.data[[col_name_orig]]), linewidth=3, color='red', alpha=0.4) + 
        # geom_point(alpha=0.3, position=position_jitter(height=0, width=0.2)) +
        geom_half_violin(aes(group=snr), side='l', scale='width',
                         fill='grey60', color=NA, alpha=1, trim=T) + 
        geom_half_point(aes(group=snr), side='r', alpha=0.1, color='black') + 
        ylab(paste(col_name, '(dist from true)')) + 
        facet_wrap(~ .data[[col_max]], nrow=2) + 
        theme_cowplot()
    
    return(p)
    
}


plot_bias_var <- function(df, col_name, col_name_orig, col_max){
    
    df$x <- as.numeric(df$snr)
    
    df$y <- df[, col_name] - df[, col_name_orig]
    
    df_bias_var <- df %>% 
        group_by(.data[[col_max]], snr) %>% 
        summarise(bias=mean(y), 
                  bias_ci=sd(y) / sqrt(n()) * qnorm(1 - 0.05), 
                  variance=var(y))
    
    pos <- position_dodge(0.1)
    p1 <- ggplot(df_bias_var, aes(snr, bias, group=.data[[col_max]], color=.data[[col_max]])) + 
        geom_hline(yintercept=0, linewidth=1, color='red', alpha=1, linetype='dotted') + 
        geom_point(position=pos) + 
        geom_errorbar(aes(ymin=bias-bias_ci, ymax=bias+bias_ci),
                      position=pos, linewidth=1, width=0) + 
        geom_line(position=pos) + 
        col_scale + 
        theme_cowplot() + 
        theme(
            legend.position = 'none'
        )
    
    p2 <- ggplot(df_bias_var, aes(snr, variance, group=.data[[col_max]], color=.data[[col_max]])) + 
        geom_hline(yintercept=0, linewidth=1, color='red', alpha=1, linetype='dotted') + 
        geom_point() + 
        geom_line() + 
        col_scale + 
        theme_cowplot()
    
    p <- ggarrange(p1, p2, nrow=2, common.legend=TRUE, legend="right")
    
    return(p)
}



plot_marginal_distribution <- function(df, var, color, alpha=0.5, bins=30) {
    ggplot(df, aes(x=.data[[var]])) +
        geom_histogram(bins=bins, fill=color, alpha=alpha, position="identity") +
        # geom_density(alpha=alpha, size = 0.1, fill=color, color=NA) +
        guides(fill='none') +
        theme_void() +
        theme(plot.margin = margin())
}


plot_correlation_scatter <- function(df, var_x, var_y, corr_result=NULL, color=col_eeg, ...){
    # ... are kwargs for the marginal plot
    
    label_pos_x <- max(df[,var_x])
    label_pos_y <- max(df[,var_y])
    if (is.null(corr_result)) {
        label_text <- ''
    } else {
        corr_method <- corr_result$method
        if (corr_method == "Spearman's rank correlation rho") corr_symbol <- 'rho'
        if (corr_method == "Kendall's rank correlation tau") corr_symbol <- 'tau'
        if (corr_method == "Pearson's product-moment correlation") corr_symbol <- 'r'
        label_text <- sprintf('%s=%.2f (p=%s)', 
                              corr_symbol,
                              corr_result$estimate, 
                              format.pval(corr_result$p.value, digits=2, eps=1e-3))
    }
    # prepare scatterplot of the two variables
    scatterplot <- ggplot(df, aes(.data[[var_x]], .data[[var_y]])) + 
        geom_point(shape=19, size=3, alpha=0.3, col=color) + 
        # scale_color_manual(values=col_generator) + 
        annotate('text', x=label_pos_x, y=label_pos_y, label=label_text, hjust='right') + 
        theme_cowplot() 
    # prepare marginal histograms
    x_hist <- plot_marginal_distribution(df, var_x, color, alpha=0.8, ...)
    y_hist <- plot_marginal_distribution(df, var_y, color, alpha=0.8, ...) + coord_flip()
    # align histograms with scatterplot
    aligned_x_hist <- align_plots(x_hist, scatterplot, align = "v")[[1]]
    aligned_y_hist <- align_plots(y_hist, scatterplot, align = "h")[[1]]
    # arrange plots
    plt <- plot_grid(
        aligned_x_hist, 
        NULL, 
        scatterplot,
        aligned_y_hist,
        ncol=2,
        nrow=2,
        rel_heights=c(0.2, 1),
        rel_widths=c(1, 0.2)
    )
    return(plt)
}



plot_subtr_effect <- function(df, pats, method) {
    # function that plots meter zscore as a function of SNR, to compare the estimate without/with 
    # noise subtraction
    df <- filter(df, pat %in% {{pats}})
    df_summary <- summarySE(data=df,
                            measurevar='z_meter', 
                            groupvars=c('selected_for', 'method', 'snr', 'pat', 'subtr'))
    pd_summary <- position_dodge(0.5)
    cols <- c('no'='#2b7557', 'yes'='#711885')
    ggplot(filter(df, selected_for=={{method}} & method=={{method}}), aes(snr, z_meter, color=subtr)) +
        geom_hline(yintercept = 0, color='#000000') + 
        geom_hline(aes(yintercept=z_meter_orig), color='#5e5e5e', linetype='dotted', linewidth=1) + 
        geom_point(alpha=0.05, position=position_jitter(height=0, width=0.1)) + 
        geom_point(data=filter(df_summary, selected_for=={{method}} & method=={{method}}), 
                   size=2, position=pd_summary) + 
        geom_errorbar(data=filter(df_summary, selected_for=={{method}} & method=={{method}}), 
                      aes(ymin=z_meter-ci, ymax=z_meter+ci), position=pd_summary, width=0, linewidth=1) + 
        scale_color_manual(name='subtracted', values=cols) + 
        scale_x_discrete(labels=round(as.numeric(levels(df$snr)), 2)) + 
        scale_y_continuous(limits=c(-1, 1)) + 
        theme_cowplot()
}



plot_dist_to_truth <- function(df, 
                               init_params_fft=list('min_'=0, 'max_'=1, 'x50_'=1, 'slope_'=1.5),
                               init_params_acf=list('min_'=0, 'max_'=1, 'x50_'=0.4, 'slope_'=1.2),
                               df_eeg_z_snr=NULL){
    # convert to factors
    df$pat <- factor(df$pat)
    
    # log transform SNR values 
    df$snr_log <- log10(df$snr)
    
    # fix negative values before taking log 
    df <- filter(df, z_snr > 0)
    df$z_snr_log <- log10(df$z_snr)
    
    # # fisher-transform the correlation coefficients
    # df$r_acf_raw <- r_to_z(df$r_acf_raw)
    # df$r_acf_subtr <- r_to_z(df$r_acf_subtr)
    # df$r_fft_raw <- r_to_z(df$r_fft_raw)
    # df$r_fft_subtr <- r_to_z(df$r_fft_subtr)
    
    var_name <- 'z_snr_log'
    
    # fit sigmoids
    m_acf <- fit_sigm(df[,var_name], df$r_acf_subtr, init_params=init_params_acf)
    m_fft <- fit_sigm(df[,var_name], df$r_fft_subtr, init_params=init_params_fft)
    
    # print confidence intervals on paramter estimates 
    confint(m_acf)
    confint(m_fft)
    
    # plot 
    p <- ggplot() +
        geom_hline(yintercept = 0, color='#000000') + 
        geom_point(data=df, aes(.data[[var_name]], r_acf_subtr), color='#bd3a3a', alpha=0.2) + 
        geom_point(data=df, aes(.data[[var_name]], r_fft_subtr), color='#2987b3', alpha=0.2) + 
        geom_function(fun=sigm, args=coef(m_acf), aes(color='acf'), linewidth=2) +
        geom_function(fun=sigm, args=coef(m_fft), aes(color='fft'), linewidth=2) +
        scale_color_manual(name='method', values=c('acf'='#802626', 'fft'='#164961')) + 
        ylab('r with ground truth') + 
        theme_cowplot() +
        scale_x_reverse() 
    
    if (!is.null(df_eeg_z_snr)) {
        df_eeg_z_snr$z_snr_log <- log10(df_eeg_z_snr$z_snr)
        df_summary <- summarySE(df_eeg_z_snr, measurevar='z_snr_log', na.rm=TRUE)
        # df_summary$z_snr_log <- log10(df_summary$z_snr_log)
        # df_summary$ci <- log10(df_summary$ci)
        p <- p + 
            geom_vline(xintercept=df_summary$z_snr_log, color='black', linewidth=1, linetype='dotted') + 
            annotate('rect', 
                     xmin=df_summary$z_snr_log-df_summary$ci, 
                     xmax=df_summary$z_snr_log+df_summary$ci, 
                     ymin=-Inf, 
                     ymax=Inf, 
                     alpha=0.1, 
                     fill='black')
    }    
    
    return(p)
}





plot_dist_to_truth_density <- function(df, y_var_name, z_snr_var_name='z_snr',
                                       df_eeg_z_snr=NULL, z_snr_var_name_eeg='z_snr'){
    
    # convert to factors
    df$pat <- factor(df$pat)
    
    # log transform SNR values 
    df$snr_log <- log10(df$snr)
    
    # fix negative values before taking log 
    df <- filter(df, {{z_snr_var_name}} > 0)
    df$z_snr_log <- log10(df[, z_snr_var_name])
    
    var_name <- 'z_snr_log'
    
    if (!is.null(df_eeg_z_snr)) {
        # this may produce nans (if z_snr < 0)
        df_eeg_z_snr$z_snr_log <- log10(df_eeg_z_snr[, z_snr_var_name_eeg])
        df_summary <- summarySE(df_eeg_z_snr, measurevar='z_snr_log', na.rm=TRUE)
        bin_step = 4/12
        bin_centres <- sort(unique(c(seq(from=df_summary$z_snr_log - bin_step/2, to=-1, by=-bin_step), 
                                     seq(from=df_summary$z_snr_log + bin_step/2, to=3, by=+bin_step))))
        title_str <- sprintf('mean log z_snr in eeg: %.2f', df_summary$z_snr_log)
    } else {
        bin_centres <- seq(from=-1, to=3, length.out=12)
        title_str <- ''
    }
    
    bin_labels <- format(bin_centres[-length(bin_centres)] + diff(bin_centres) / 2, digits=1)
    df$z_snr_log_binned <- cut(df$z_snr_log, bin_centres)
    df <- df[!is.na(df$z_snr_log_binned), ]
    
    ggplot(df, aes(x=.data[[y_var_name]], y=z_snr_log_binned, fill=stat(x))) + 
        geom_density_ridges_gradient(scale=0.9, size=0.3, rel_min_height=0.001, 
                                     jittered_points=TRUE,
                                     position=position_points_jitter(width=0.0, height=0.05),
                                     point_size=0.4, point_alpha=0.3, point_shape=20, alpha=0.7) + 
        scale_fill_viridis_c(name="r", option="C", limits=c(-1.2, 1.2), breaks=c(-1, 0, 1)) +
        # scale_fill_gradientn(colours=coolwarm(100), limits=c(-1.2, 1.2), breaks=c(-1, 0, 1)) + 
        scale_y_discrete(limits=rev, labels=rev(bin_labels)) + 
        scale_x_continuous(limits=c(-1.2, 1.2)) + 
        theme_ridges() + 
        coord_flip() + 
        geom_vline(xintercept=0, color='black', linewidth=1) + 
        ggtitle(title_str)
    
}


plot_feat_ind_lowhigh <- function(df, 
                          col_name_eeg='z_meter_fft_subtr', col_name_coch='z_meter_fft_sound'){
    
    df$x <- as.numeric(df$tone)
    df$x_coch_start <- df$x-0.3
    df$x_coch_end <- df$x+0.3
    
    df_summary <- summarySE(df, groupvars=c('rhythm', 'tone'), measurevar=col_name_eeg)
    
    df_summary$x <- as.numeric(df_summary$tone)
    
    pd <- position_dodge(0.0)
    
    p <- ggplot(df, aes(x, .data[[col_name_eeg]], color=tone)) + 
        geom_segment(data=df[df$subject==df$subject[1],],
                     inherit.aes=FALSE,
                     aes(x=x_coch_start, xend=x_coch_end,
                         y=.data[[col_name_coch]], yend=.data[[col_name_coch]],
                         color=tone),
                     linewidth=3, alpha=0.5) +
        scale_color_manual(name='tone', values=cols) +
        scale_x_continuous(breaks=c(1:length(levels(df$tone))), 
                           labels=levels(df$tone)) +  
        scale_y_continuous(breaks=c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5)) +
        new_scale_color() + new_scale_fill() +
        geom_line(aes(group=paste(subject, rhythm)), color='grey80', position=pd) + 
        geom_point(aes(fill=tone, color=tone, group=paste(subject, rhythm)), size=2, shape=21, position=pd) +
        scale_color_manual(name='tone', values=darken(cols_individual, 0.1)) + 
        scale_fill_manual(name='tone', values=cols_individual) + 
        geom_hline(yintercept=0) +
        new_scale_color() + new_scale_fill() + 
        geom_point(data=df_summary, aes(color=tone), size=3) + 
        geom_errorbar(data=df_summary, 
                      aes(color=tone, ymin=.data[[col_name_eeg]]-ci,
                          ymax=.data[[col_name_eeg]]+ci), 
                      linewidth=1, width=0.2) + 
        scale_color_manual(name='tone', values=cols) + 
        facet_wrap(~rhythm, nrow=1) + 
        theme_cowplot() + 
        my_theme() + 
        theme(axis.ticks.y = element_line(), 
              legend.position = 'none',
              strip.text = element_blank(),
              strip.background = element_blank())
    
    return(p)
}

