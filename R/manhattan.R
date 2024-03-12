# silence R CMD checks for data.table columns
BETA = BP = CHR = P = SE = SNP = chr_len = tot = x = i.tot = highlight = secondary_highlight = NULL

#' @title Manhattan plot
#' @description
#' Create a Manhattan plot with ggplot2 geom_point.
#'
#' @param gwas a data.table with a minimum of columns SNP, CHR, BP, and P
#' @param highlight_snps (optional) a character vector of SNPs to highlight
#' @param highlight_win (optional) a numeric, the number of kb either side of the highlight_snps to also highlight (i.e create peaks)
#' @param annotate_snps (optional) a character vector of SNPs to annotate
#' @param colours (optional) a character vector colour codes to be replicated along the chromosomes
#' @param highlight_colour (optional) a character colour code; the colour to highlight points in
#' @param highlight_shape (optional) a numeric shape code; the shape of the highlight points (see ggplot2 shape codes)
#' @param highlight_alpha (optional) a numeric value between 0 and 1; the alpha of the highlighted points colour
#' @param sig_line_1 (optional) a numeric value (-log10(P)) for where to draw a horizontal line
#' @param sig_line_2 (optional) a numeric value (-log10(P)) for where to draw a second horizontal line
#' @param y_limits (optional) a numeric length 2 vector c(min-Y, max-Y)
#' @param title (optional) a string title
#' @param subtitle (optional) a string subtitle
#' @param hit_table (optional) a logical, whether to display a table of top hits (lowest P values)
#' @param max_table_hits (optional) an integer, how many top hits to show in the table
#' @param downsample (optional) a numeric between 0 and 1, tghe proportion by which to downsample by, e.g. 0.5 will remove 50% of points above the downsample_pval threshold (can help increase plotting speed with minimal impact on plot appearance)
#' @param base_text_size an integer, `base_size` for the ggplot2 theme
#' @param downsample_pval (optional) a numeric between 0 and 1, the p-values affected by downsampling, default >0.1
#' @return a ggplot
#' @export
#' @importFrom data.table setkey
#' @importFrom ggpubr ggarrange
#' @importFrom ggrepel geom_label_repel
#' @importFrom gridExtra tableGrob
#' @import ggplot2
#'
manhattan <- function(gwas,
                      highlight_snps = NULL,
                      highlight_win = 100,
                      annotate_snps = NULL,
                      colours = c("#d9d9d9", "#bfbfbf"),
                      highlight_colour = "#e15758",
                      highlight_shape = 16,
                      highlight_alpha = 1.0,
                      sig_line_1 = 5e-8,
                      sig_line_2 = NULL,
                      y_limits = c(NULL,NULL),
                      title = NULL,
                      subtitle = NULL,
                      base_text_size = 14,
                      hit_table = FALSE,
                      max_table_hits = 10,
                      downsample = 0.1,
                      downsample_pval = 0.1) {

  # checks
  stopifnot("gwas input must be a data.frame like object" = inherits(gwas, "data.frame"))
  stopifnot("Input data must contain at least columns: snp, chr, bp, and p" = all(c("snp","chr","bp","p") %in% colnames(gwas)))
  stopifnot("highlight_snps should be a character vector, a subset of the 'snp' column" = is.null(highlight_snps) | is.character(highlight_snps))
  stopifnot("Annotate should be a character vector, a subset of the 'snp' column" = is.null(annotate_snps) | is.character(annotate_snps))
  stopifnot("highlight_win must be numeric (kb either side of each hit to highlight; default=0)" = is.numeric(highlight_win))
  stopifnot("downsample percentage must be numeric between 0 and 1" = is.numeric(downsample))
  stopifnot("downsample_pval must be numeric between 0 and 1" = is.numeric(downsample_pval))

  # #testing
  # gwas = genepi.utils::generate_random_gwas_data(100000)
  # highlight_snps = c("15:135277414[b37]G,T")
  # annotate_snps = c("15:135277414[b37]G,T")

  # create factor from CHR
  gwas[, chr := factor(chr, levels=c(as.character(1:25), setdiff(as.character(unique(chr)), as.character(1:25))))]

  # down-sample the insignificant dots for plotting
  if(downsample > 0) {
    set.seed(2023)
    insig_rows <- which(gwas[["p"]] > downsample_pval)
    insig_pval <- gwas[insig_rows, ][["p"]]
    prob       <- (insig_pval - min(insig_pval, na.rm=T)) / (max(insig_pval, na.rm=T) - min(insig_pval, na.rm=T)) # high P more likely to be removed
    n_remove   <- round(length(insig_pval)*downsample)
    remove_idxs<- sample(insig_rows, size=n_remove, prob=prob, replace=TRUE) # replace=FALSE takes ages, so just call unique after and accept exact % wont be quite accurate
    gwas <- gwas[-unique(remove_idxs), ]
  }

  # prepare x axis
  data.table::setkey(gwas, "chr", "bp")
  gwas[gwas[, list(chr_len = as.numeric(max(bp))), by = "chr"]
           [, tot := data.table::shift(cumsum(chr_len), fill=0)],
       x := bp + i.tot,
       on = "chr"]
  x_ticks <- gwas[, list("x_ticks"= (max(x) + min(x)) / 2), by = "chr"]

  # prepare y axis
  if(is.null(y_limits)) {
    # data max
    max_p <- max(-log10(gwas$p),na.rm=T)
    # data & line1 max
    if(!is.null(sig_line_1)) {
      max_p <- max(c(max_p, -log10(sig_line_1)),na.rm=T)
    }
    # data & line2 max
    if(!is.null(sig_line_2)) {
      max_p <- max(c(max_p, -log10(sig_line_2)),na.rm=T)
    }
    y_limits <- c(0, ceiling(max_p))
  }
  log10P <- expression(paste("-log"[10], plain(P)))

  # base plot
  plot <- ggplot2::ggplot(gwas, ggplot2::aes(x=x, y=-log10(p), color=chr)) +
    ggplot2::geom_point(alpha=0.8, size=0.2) +
    ggplot2::scale_x_continuous(label=x_ticks$chr, breaks=x_ticks$x_ticks, expand=c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(limits=y_limits, expand=c(0, 0)) +
    ggplot2::scale_color_manual(values=rep(colours, length.out=length(levels(gwas$chr)))) +
    ggplot2::theme_classic(base_size = base_text_size) +
    ggplot2::theme(
      legend.position ="none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    )

  # add highlighted variants within requested kb BP window
  if(!is.null(highlight_snps) & length(highlight_snps)>0) {
    gwas[snp %in% highlight_snps, highlight := TRUE]
    if(highlight_win > 0) {
      highlight_data <- gwas[snp %in% highlight_snps, ]
      for(i in 1:nrow(highlight_data)) {
        bp_  <- highlight_data[[i, "bp"]]
        chr_ <- highlight_data[[i, "chr"]]
        gwas[bp > (bp_ - highlight_win*1000) & bp < (bp_ + highlight_win*1000) & chr == chr_, secondary_highlight := TRUE]
      }
    }
    plot <- plot +
      ggplot2::geom_point(data = gwas[secondary_highlight==TRUE, ], color=highlight_colour, fill=highlight_colour, alpha=highlight_alpha) +
      ggplot2::geom_point(data = gwas[highlight==TRUE, ], color=highlight_colour, fill=highlight_colour, shape=highlight_shape, alpha=highlight_alpha, size=3)
  }

  # add horizontal significance lines
  if(!is.null(sig_line_1)) {
    plot <- plot + ggplot2::geom_hline(yintercept = -log10(sig_line_1), linetype = "dashed", color = "darkgrey")
  }
  if(!is.null(sig_line_2)) {
    plot <- plot + ggplot2::geom_hline(yintercept = -log10(sig_line_2), linetype = "dashed", color = "darkgrey")
  }

  # add SNP label annotations
  if(!is.null(annotate_snps) & length(annotate_snps)>0) {
    gwas[SNP %in% annotate_snps, annotate := TRUE]
    label_x_nudge <- max(gwas[["x"]], na.rm=TRUE) / 22
    plot <- plot +
      ggrepel::geom_label_repel(data=gwas[annotate==TRUE, ], ggplot2::aes(label=SNP), colour="black", label.size = 0.1, nudge_x=label_x_nudge )
  }

  # add titles and labels
  plot <- plot + ggplot2::labs(title=title, subtitle=subtitle, x="Chromosome", y=log10P)

  # add the hit_table if requested
  if(hit_table) {
    table <- hit_table(gwas, max_table_hits)
    plot <-  ggpubr::ggarrange(plot, table, nrow=1, widths=c(0.75, 0.25))
  }

  # return the plot
  return(plot)
}


#' @title Create top hits table
#' @param gwas a data.table with at least columns SNP & P; optionally BETA and SE
#' @param n number of top hits to display
#' @return a gridExtra::tableGrob
#' @importFrom utils head
#' @noRd
#'
hit_table <- function(gwas, n) {
  # base columns and BETA and SE if provided
  cols <- c(c("snp","p"), names(gwas)[names(gwas) %in% c("beta","se")])
  table_data <- utils::head(gwas[order(p), cols, with=FALSE], n=n)
  table_data <- table_data[, lapply(.SD, signif, digits=2), by=snp]
  table_data[, snp := strtrim(snp, 23)]
  table <- gridExtra::tableGrob(table_data, rows=NULL)
  return(table)
}


#' @title Miami plot
#' @description
#' Create a Miami plot. Please look carefully at the parameters as these largely map to the
#' `manhattan()` parameters, the main difference being that you need to supply a 2 element
#' list of the parameter, one for the upper and one for the lower plot aspect of the Miami
#' plot. Some parameters are not duplicated however - see the example defaults below.
#' @param gwases a list of 2 data.tables
#' @inheritParams manhattan
#' @return a ggplot
#' @export
#' @importFrom ggpubr annotate_figure
#'
miami <- function(gwases,
                  highlight_snps   = list("top"=NULL, "bottom"=NULL),
                  highlight_win    = list("top"=100, "bottom"=100),
                  annotate_snps    = list("top"=NULL,"bottom"=NULL),
                  colours          = list("top"=c("#d9d9d9","#bfbfbf"),"bottom"=c("#bfbfbf","#d9d9d9")),
                  highlight_colour = list("top"="#e15758","bottom"="#4f79a7"),
                  highlight_shape  = list("top"=16,"bottom"=16),
                  sig_line_1       = list("top"=5e-8,"bottom"=5e-8),
                  sig_line_2       = list("top"=NULL,"bottom"=NULL),
                  y_limits         = list("top"=c(NULL,NULL),"bottom"=c(NULL,NULL)),
                  title            = NULL,
                  subtitle         = list("top"=NULL,"bottom"=NULL),
                  base_text_size   = 14,
                  hit_table        = FALSE,
                  max_table_hits   = 10,
                  downsample       = 0.1,
                  downsample_pval  = 0.1) {

  # checks
  stopifnot("Ensure inputs are lists of length two (see documentation for exceptions)" = all(sapply(list(gwases,highlight_snps,highlight_win,annotate_snps,colours,highlight_colour,highlight_shape,sig_line_1,sig_line_2,y_limits),
                                                                                             function(x) is.list(x) & length(x)==2) ))

  # create upper plot
  plot_upper <- manhattan(gwas             = gwases[[1]],
                          highlight_snps   = highlight_snps[[1]],
                          highlight_win    = highlight_win[[1]],
                          annotate_snps    = annotate_snps[[1]],
                          colours          = colours[[1]],
                          highlight_colour = highlight_colour[[1]],
                          highlight_shape  = highlight_shape[[1]],
                          sig_line_1       = sig_line_1[[1]],
                          sig_line_2       = sig_line_2[[1]],
                          y_limits         = y_limits[[1]],
                          hit_table        = FALSE,
                          max_table_hits   = max_table_hits,
                          downsample       = downsample,
                          downsample_pval  = downsample_pval)

  # create lower plot
  plot_lower <- manhattan(gwas             = gwases[[2]],
                          highlight_snps   = highlight_snps[[2]],
                          highlight_win    = highlight_win[[2]],
                          annotate_snps    = annotate_snps[[2]],
                          colours          = colours[[2]],
                          highlight_colour = highlight_colour[[2]],
                          highlight_shape  = highlight_shape[[2]],
                          sig_line_1       = sig_line_1[[2]],
                          sig_line_2       = sig_line_2[[2]],
                          y_limits         = y_limits[[2]],
                          hit_table        = FALSE,
                          max_table_hits   = max_table_hits,
                          downsample       = downsample,
                          downsample_pval  = downsample_pval)

  # remove x axis labels
  plot_upper <- plot_upper +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  # flipping overwrites manhattan() things so need to recalculate
  if(is.null(y_limits[[2]])) {
    # data max
    max_p <- max(-log10(gwases[[2]]$p),na.rm=T)
    # data & line1 max
    if(!is.null(sig_line_1[[2]])) {
      max_p <- max(c(max_p, -log10(sig_line_1[[2]])),na.rm=T)
    }
    # data & line2 max
    if(!is.null(sig_line_2[[2]])) {
      max_p <- max(c(max_p, -log10(sig_line_2[[2]])),na.rm=T)
    }
    y_limits[[2]] <- c(ceiling(max_p), 0)
  } else {
    # user usually provides in the normal way c(0, upper)
    y_limits[[2]] <- rev(y_limits[[2]])
  }

  # flip the bottom plot
  plot_lower <- plot_lower +
    ggplot2::scale_y_reverse(limits=y_limits[[2]], expand=c(0, 0)) +
    ggplot2::scale_x_continuous(expand=c(0.01, 0.01), position = "top") +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.line.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   plot.background=ggplot2::element_rect(fill='transparent', color=NA),
                   panel.background=ggplot2::element_rect(fill='transparent'))

  # combine plots, with or without tables
  if(hit_table) {
    upper_table <- hit_table(gwases[[1]], max_table_hits)
    lower_table <- hit_table(gwases[[2]], max_table_hits)
    plot <- ggpubr::ggarrange(plot_upper,NULL,upper_table,NULL,NULL,NULL, plot_lower,NULL,lower_table,
                              nrow=3, ncol=3, heights=c(1,0,1), widths=c(0.75,0,0.25),
                              labels=c(subtitle[[1]],"","",subtitle[[2]],"",""),
                              vjust = 0.5,
                              font.label=list(size=12, color="black", face="plain"))
  } else {
    plot <- ggpubr::ggarrange(plot_upper, NULL, plot_lower, ncol=1, heights=c(1,0,1),
                              labels = c(subtitle[[1]],"",subtitle[[2]]),
                              vjust = 0.5,
                              font.label=list(size=12, color="black", face="plain"))
  }

  # add overall title
  if(!is.null(title)) {
    plot <- ggpubr::annotate_figure(plot, top=ggpubr::text_grob(title, face="bold", size=base_text_size+2))
  }

  # return the complete plot
  return(plot)
}


