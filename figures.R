
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


plot_feature_per_factor <- function(df, col_name_x='rhythm', col_name_wrap=NULL, col_name_id_var='ses',
                                    SE_type='within',  col_name_eeg='z_meterRel',
                                    col_name_coch='z_meterRel_hilbert', point_alpha=0.2){
    if (is.null(col_name_x)) {
        df$tmp <- ''
        col_name_x <- 'tmp'
    }
    
    if (!is.factor(df[,col_name_x])) df[,col_name_x] <- factor(df[,col_name_x])
    x_labs <- levels(df[,col_name_x])
                     
    df$x <- as.numeric(df[,col_name_x])
    df$x_coch_start <- df$x-0.3
    df$x_coch_end <- df$x+0.3
    
    if (length(x_labs) == 1 & is.null(col_name_wrap)){
        df_summary <- summarySE(df, measurevar=col_name_eeg, na.rm=TRUE)
        df_summary[,col_name_x] <- factor('')
    } else {
        if (SE_type == 'within') {
            df_summary <- summarySEwithin(df,
                                          withinvars=c(col_name_x, col_name_wrap), 
                                          idvar=col_name_id_var, 
                                          measurevar=col_name_eeg, 
                                          na.rm=TRUE)
        } else if (SE_type == 'between') {
            df_summary <- summarySE(df,
                                  groupvars=c(col_name_x, col_name_wrap), 
                                  measurevar=col_name_eeg, 
                                  na.rm=TRUE)
        } else {
            stop('SE type not recognized...')
        }
    }

    df_summary$x <- as.numeric(df_summary[,col_name_x])
    
    p <- ggplot(df, aes(x, .data[[col_name_eeg]]))
    if (!is.null(col_name_coch)){
        p <- p + 
            geom_segment(data=df[df$ses==df$ses[1], ],
                         inherit.aes=FALSE,
                         aes(x=x_coch_start, xend=x_coch_end, 
                             y=.data[[col_name_coch]], yend=.data[[col_name_coch]]),
                         linewidth=3, color='black', alpha=0.5) 
    }
    p <- p + 
        geom_point(color=alpha(col_eeg, point_alpha), fill=alpha(col_eeg, point_alpha),
                   size=3, stroke=0, shape=19, position=position_jitter(width=0.1, height=0)) +
        scale_x_continuous(breaks=c(1:length(x_labs)), labels=x_labs) +  
        geom_hline(yintercept=0) +
        new_scale_color() + new_scale_fill() + 
        geom_point(data=df_summary, color='black', size=3) + 
        geom_errorbar(data=df_summary, 
                      aes(ymin=.data[[col_name_eeg]]-ci, ymax=.data[[col_name_eeg]]+ci), 
                      size=1, width=0, color='black') + 
        theme_cowplot() + 
        my_theme() + 
        theme(axis.text.x = element_text(size=fontsize, angle=30, hjust=1, vjust=1))
    if (!is.null(col_name_wrap)){
        p <- p + 
            facet_wrap(~.data[[col_name_wrap]]) + 
            # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
            annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
    }
    return(p)
}




plot_pairwise_matrix <- function(p_matrix){
    # This function plots matrix with pvalues corresponding to pairwise comparisons between conditions. 
    df <- as.data.frame.table(p_matrix)
    names(df)[3] <- 'p'
    df$p_text <- format.pval(df$p, digits=3, eps=1e-1)
    df$p_text[df$p_text=='NA'] <- NA
    fontsize <- ifelse(length(unique(df$Var2)) > 9, 3, 6)
    df %>%
        ggplot(aes(Var1, Var2, fill=p)) + 
        geom_tile(color='white', lwd=0.1, linetype=1) +
        geom_text(aes(label=p_text), size=fontsize) + 
        scale_fill_gradient2(mid="red", limit=c(0,0.05),name="p") + 
        theme_cowplot() + 
        theme(
            axis.title = element_blank(), 
            axis.text.x = element_text(angle=60, hjust=1)
        )  + 
        coord_fixed()
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
        geom_point(shape=19, size=1, alpha=0.3, col=color) + 
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


plot_feat_over_time <- function(df, col_name_eeg, col_name_time='date', x_jitter=0, y_jitter=0, 
                                point_alpha=0.1, ncol=1, nrow=NULL){
    p <- ggplot(df, aes(.data[[col_name_time]], .data[[col_name_eeg]])) + 
        geom_point(alpha=point_alpha, size=2, color=col_eeg, 
                   position=position_jitter(width=x_jitter, height=y_jitter)) + 
        facet_wrap(~rhythm, ncol=ncol, nrow=nrow) + 
        theme_cowplot() + 
        my_theme() + 
        theme(
            axis.text.x = element_text(size=fontsize, angle=30, hjust=1, vjust=1),
            axis.line.x = element_line(),
            axis.line.y = element_line()
        )
    if (is.Date(df[, col_name_time])) {
        p <- p + 
            scale_x_date(date_labels = "%b %Y") 
    }
    return(p)
}


plot_feat_over_trials <- function(df, yvar, m_summary=NULL){
    
    p <- plot_feat_over_time(df, yvar, col_name_time='trial',
                             x_jitter=0.2, point_alpha=0.01, ncol=2) + 
        scale_x_continuous(breaks=c(1:7)) + 
        theme(
            axis.text.x = element_text(angle=0)
        )
    if (!is.null(m_summary)){
        df_line_fit <- data.frame(
            rhythm = c('unsyncopated', 'syncopated'), 
            intercept = c(m_summary$coefficients['(Intercept)', 'Estimate'], 
                          m_summary$coefficients['(Intercept)', 'Estimate'] + 
                              m_summary$coefficients['rhythmsyncopated', 'Estimate']), 
            slope = c(m_summary$coefficients['trial', 'Estimate'], 
                      m_summary$coefficients['trial', 'Estimate'] + 
                          m_summary$coefficients['trial:rhythmsyncopated', 'Estimate'])
        )
        df_line_fit <- fix_factors(df_line_fit)
        
        p <- p + geom_abline(data=df_line_fit, aes(intercept=intercept, slope=slope), 
                             color='grey50', linetype='solid', linewidth=2) 
    }
    return(p)
}


plot_feat_over_sliding_win <- function(df, yvar, m_summary=NULL){
    
    win_times <- unique(df$win_time)
    
    p <- plot_feat_over_time(df, yvar, col_name_time='win_time',
                             x_jitter=0.2, point_alpha=0.01, ncol=2) + 
        scale_x_continuous(breaks=c(win_times[1], win_times[length(win_times)])) + 
        theme(
            axis.text.x = element_text(angle=0)
        )
    if (!is.null(m_summary)){
        df_line_fit <- data.frame(
            rhythm = c('unsyncopated', 'syncopated'), 
            intercept = c(m_summary$coefficients['(Intercept)', 'Estimate'], 
                          m_summary$coefficients['(Intercept)', 'Estimate'] + 
                              m_summary$coefficients['rhythmsyncopated', 'Estimate']), 
            slope = c(m_summary$coefficients['win_time', 'Estimate'], 
                      m_summary$coefficients['win_time', 'Estimate'] + 
                          m_summary$coefficients['win_time:rhythmsyncopated', 'Estimate'])
        )
        df_line_fit <- fix_factors(df_line_fit)
        
        p <- p + geom_abline(data=df_line_fit, aes(intercept=intercept, slope=slope), 
                             color='grey50', linetype='solid', linewidth=2) 
    }
    return(p)
}
