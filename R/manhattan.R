gwas_in = data.table::fread("/Users/xx20081/Documents/local_data/hermes_progression/biostat_disc/raw/BIOSTAT_Discovery.allcause.gz")
set.seed(2023)
gwas_in = gwas_in[sample(1:nrow(gwas_in), size=100000), c("CHR","POS","OTHER_ALLELE","EFFECT_ALLELE","EAF","BETA","SE","P","EUR")]
gwas_in[, BETA := gwas_in[sample(nrow(gwas_in)), BETA]]
gwas_in[, P    := gwas_in[sample(nrow(gwas_in)), P]]
gwas_in[, P    := gwas_in[sample(nrow(gwas_in)), P]]
gwas_in[, SNP  := gwas_in[, paste0(CHR,":",POS,"[b37]",OTHER_ALLELE,",",EFFECT_ALLELE)]]
data.table::setnames(gwas_in, c("CHR","POS","OTHER_ALLELE","EFFECT_ALLELE","EAF","BETA","SE","P","EUR"), c("CHR","BP","OA","EA","EAF","BETA","SE","P","EUR_EAF"))
data.table::setkey(gwas_in, CHR, BP)
gwas_path = "/Users/xx20081/git/genepi.utils/inst/extdata/example3_gwas_sumstats.tsv"
data.table::fwrite(gwas_in, gwas_path, sep="\t")
gwas_in = data.table::fread(gwas_path)

highlight_snps = gwas_in[SNP=="4:32205845[b37]C,T", ][["SNP"]]
annotate_snps = gwas_in[SNP=="4:32205845[b37]C,T", ][["SNP"]]


p <- manhattan(gwas,
               highlight_snps = highlight_snps,
               highlight_win = 250,
               annotate_snps = annotate_snps,
               hit_table = TRUE,
               title = "HELLO",
               subtitle = "subdsadsda")
p


m <- miami(gwases           = list(gwas, data.table::copy(gwas)),
           highlight_snps   = list("top"=highlight_snps, "bottom"=highlight_snps),
           highlight_win    = list("top"=100, "bottom"=100),
           annotate_snps    = list("top"=annotate_snps,"bottom"=annotate_snps),
           title            = "Miami",
           hit_table        = TRUE)
m






manhattan <- function(gwas,
                      highlight_snps = NULL,
                      highlight_win = 100,
                      annotate_snps = NULL,
                      colours = c("#d9d9d9", "#bfbfbf"),
                      highlight_colour = "#e15758",
                      highlight_shape = 16,
                      sig_line_1 = 5e-8,
                      sig_line_2 = NULL,
                      y_limits = c(0, 7.5),
                      title = NULL,
                      subtitle = NULL,
                      hit_table = FALSE,
                      max_table_hits = 10,
                      downsample = 0.1,
                      downsample_pval = 0.1) {

  # checks
  stopifnot("gwas input must be a data.frame like object" = inherits(gwas, "data.frame"))
  stopifnot("Input data must contain at least columns: SNP, CHR, BP, and P" = all(c("SNP","CHR","BP","P") %in% colnames(gwas)))
  stopifnot("highlight_snps should be a character vector, a subset of the 'SNP' column" = is.null(highlight_snps) | is.character(highlight_snps))
  stopifnot("Annotate should be a character vector, a subset of the 'SNP' column" = is.null(annotate_snps) | is.character(annotate_snps))
  stopifnot("highlight_win must be numeric (kb either side of each hit to highlight; default=0)" = is.numeric(highlight_win))
  stopifnot("downsample percentage must be numeric between 0 and 1" = is.numeric(downsample))
  stopifnot("downsample_pval must be numeric between 0 and 1" = is.numeric(downsample_pval))

  # create factor from CHR
  gwas[, CHR := as.factor(CHR)]
  gwas[, CHR_int := as.integer(CHR)]

  # down-sample the insignificant dots for plotting
  if(downsample > 0) {
    set.seed(2023)
    insig_rows <- which(gwas[["P"]] > downsample_pval)
    insig_pval <- gwas[insig_rows, ][["P"]]
    gwas <- gwas[-sample(x    = insig_rows,
                         size = round(length(insig_pval)*downsample),
                         prob = (insig_pval - min(insig_pval, na.rm=T)) / (max(insig_pval, na.rm=T) - min(insig_pval, na.rm=T))), ]
  }

  # prepare x axis
  data.table::setkey(gwas, CHR_int, BP)
  gwas[gwas[, .(chr_len = as.numeric(max(BP))), by = CHR_int]
           [, tot := data.table::shift(cumsum(chr_len), fill=0)],
       x := BP + i.tot,
       on = "CHR_int"]
  x_ticks <- gwas[, .("x_ticks"= (max(x) + min(x)) / 2), by = CHR_int]

  # prepare y axis
  if(is.null(y_limits)) {
    y_limits <- c(min(-log10(gwas$P),na.rm=T), max(-log10(gwas$P),na.rm=T))
  }
  log10P <- expression(paste("-log"[10], plain(P)))

  # base plot
  plot <- ggplot2::ggplot(gwas, ggplot2::aes(x=x, y=-log10(P), color=CHR)) +
    ggplot2::geom_point(alpha=0.8, size=0.2) +
    ggplot2::scale_x_continuous(label=x_ticks$CHR_int, breaks=x_ticks$x_ticks, expand=c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(limits=y_limits, expand=c(0, 0)) +
    ggplot2::scale_color_manual(values=rep(colours, length.out=length(levels(gwas$CHR)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position ="none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    )

  # add highlighted variants within requested kb BP window
  if(!is.null(highlight_snps)) {
    gwas[SNP %in% highlight_snps, highlight := TRUE]
    if(highlight_win > 0) {
      highlight_data <- gwas[SNP %in% highlight_snps, ]
      for(i in 1:nrow(highlight_data)) {
        bp <- highlight_data[[i, "BP"]]
        chr<- highlight_data[[i, "CHR"]]
        gwas[BP > (bp - highlight_win*1000) & BP < (bp + highlight_win*1000) & CHR == chr, highlight := TRUE]
      }
    }
    plot <- plot +
      ggplot2::geom_point(data = gwas[highlight==TRUE, ], color=highlight_colour, shape=highlight_shape, alpha=0.8, size=1.5)
  }

  # add SNP label annotations
  if(!is.null(annotate_snps)) {
    gwas[SNP %in% annotate_snps, annotate := TRUE]
    label_x_nudge <- max(gwas[["x"]], na.rm=TRUE) / 22
    plot <- plot +
      ggrepel::geom_label_repel(data=gwas[annotate==TRUE, ], ggplot2::aes(label=SNP), colour="black", label.size = 0.1, nudge_x=label_x_nudge )
  }

  # add horizontal significance lines
  if(!is.null(sig_line_1)) {
    plot <- plot + ggplot2::geom_hline(yintercept = -log10(sig_line_1), linetype = "dashed", color = "darkgrey")
  }
  if(!is.null(sig_line_2)) {
    plot <- plot + ggplot2::geom_hline(yintercept = -log10(sig_line_2), linetype = "dashed", color = "darkgrey")
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

hit_table <- function(gwas, n) {
  # base columns and BETA and SE if provided
  cols <- c(c("SNP","P"), names(gwas)[names(gwas) %in% c("BETA","SE")])
  table_data <- head(gwas[order(P), cols, with=FALSE], n=n)
  table_data <- table_data[, lapply(.SD, signif, digits=2), by=SNP]
  table_data[, SNP := strtrim(SNP, 23)]
  table <- gridExtra::tableGrob(table_data, rows=NULL)
  return(table)
}


miami <- function(gwases,
                  highlight_snps   = list("top"=NULL, "bottom"=NULL),
                  highlight_win    = list("top"=100, "bottom"=100),
                  annotate_snps    = list("top"=NULL,"bottom"=NULL),
                  colours          = list("top"=c("#d9d9d9","#bfbfbf"),"bottom"=c("#bfbfbf","#d9d9d9")),
                  highlight_colour = list("top"="#e15758","bottom"="#4f79a7"),
                  highlight_shape  = list("top"=16,"bottom"=16),
                  sig_line_1       = list("top"=5e-8,"bottom"=5e-8),
                  sig_line_2       = list("top"=NULL,"bottom"=NULL),
                  y_limits         = list("top"=c(0, 7.5),"bottom"=c(0, 7.5)),
                  title            = NULL,
                  # subtitle         = list("top"=NULL,"bottom"=NULL), # doesn't look good at the minute
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

  # flip the bottom plot
  plot_lower <- plot_lower +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_continuous(expand=c(0.01, 0.01), position = "top") +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.line.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   plot.background=ggplot2::element_rect(fill='transparent', color=NA),
                   panel.background=ggplot2::element_rect(fill='transparent'))

  # combine plots
  if(hit_table) {
    upper_table <- hit_table(gwases[[1]], max_table_hits)
    lower_table <- hit_table(gwases[[2]], max_table_hits)
    plot <- ggpubr::ggarrange(plot_upper,NULL,upper_table,NULL,NULL,NULL, plot_lower,NULL,lower_table,
                              nrow=3, ncol=3, heights=c(1,-0.05,1), widths=c(0.75,0,0.25))
                              # doesnt line up.... labels=c(subtitle[[1]],NULL,NULL,NULL,NULL,NULL,subtitle[[2]],NULL,NULL))
  } else {
    plot <- ggpubr::ggarrange(plot_upper, NULL, plot_lower, ncol=1, heights=c(1, -0.05, 1))
                              # labels = subtitle)
  }

  # add overall title
  if(!is.null(title)) {
    plot <- ggpubr::annotate_figure(plot, top=ggpubr::text_grob(title, face = "bold"))
  }

  return(plot)
}


