# cat("Plotting QQ:", "[", this_study, this_outcome, "]\n")

# # gather the data for the QQ plot
# qq_dat <- dt |>
#   dplyr::filter(!is.na(P),
#                 !is.nan(P),
#                 !is.null(P),
#                 is.finite(P),
#                 P < 1.0,
#                 P > 0.0) |>
#   dplyr::mutate(chisq = qchisq(1-P, 1)) |>
#   dplyr::group_by(QC_STATUS) |>
#   dplyr::mutate(observed = -log10(sort(P, decreasing=FALSE)),
#                 expected = -log10(stats::ppoints(dplyr::n())),
#                 clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:dplyr::n(), shape2 = dplyr::n():1)),
#                 cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:dplyr::n(), shape2 = dplyr::n():1))) |>
#   dplyr::ungroup()
#
#
# # generate QQ-plot
#
# # calculate the lambda_gc and generate labels df
# lambda_gc_pre  = round( median(qq_dat$chisq[qq_dat$QC_STATUS == "PRE-QC"]) / qchisq(0.5, 1), digits=2)
# lambda_gc_post = round( median(qq_dat$chisq[qq_dat$QC_STATUS == "POST-QC"]) / qchisq(0.5, 1), digits=2)
# ann_text <- data.frame(expected = c(5,5),
#                        observed = c(2,2),
#                        label = paste0("\u03BB = ", c(lambda_gc_pre, lambda_gc_post)))
#
# # axis labels
# log10Pe <- expression(paste("Expected -log"[10], plain(P)))
# log10Po <- expression(paste("Observed -log"[10], plain(P)))
#
# # plot
# plot <- qq_dat |>
#   ggplot(aes(x = expected, y = observed)) +
#   geom_point(size = 0.5, color=wes_palette("GrandBudapest2")[4]) +
#   geom_ribbon(mapping = aes(x = expected, ymin = clower, ymax = cupper),
#               alpha = 0.1) +
#   geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
#   geom_text(data = ann_text, aes(x=expected, y=observed, label=label), size=6) +
#   labs(title = paste0("QQ plot - [", this_study, "][", this_outcome, "]"),
#        x = log10Pe,
#        y = log10Po,
#        color = NULL) +
#   theme_light(base_size = 22) +
#   theme(aspect.ratio = 1) +
#   facet_grid(cols = vars(QC_STATUS))
#
# # save plot
# grDevices::png(filename=paste0(base_fig_path, "qq_plot.png"), bg="white", height=600, width=1200, units="px")
# print(plot)
# grDevices::dev.off()
