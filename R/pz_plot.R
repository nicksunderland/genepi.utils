# pz_plot <- function(gwas,)
#
# # generate PZ-plot
# plot <- dt |>
#   dplyr::filter(!is.na(P),
#                 !is.nan(P),
#                 !is.null(P),
#                 is.finite(P),
#                 P < 1.0,
#                 P > 0.0,
#                 !is.na(BETA),
#                 !is.nan(BETA),
#                 !is.null(BETA),
#                 is.finite(BETA),
#                 !is.na(SE),
#                 !is.nan(SE),
#                 !is.null(SE),
#                 is.finite(SE)) |>
#   dplyr::group_by(QC_STATUS) |>
#   dplyr::mutate(observed = -log10(P),
#                 expected = -log10(2*pnorm(abs(BETA / SE), lower.tail=FALSE))) |>
#   dplyr::ungroup() |>
#   ggplot(aes(x = expected, y = observed)) +
#   geom_point(size = 0.5, color=wesanderson::wes_palette("GrandBudapest2")[4]) +
#   geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
#   labs(title = paste0("PZ plot - [", this_study, "][", this_outcome, "]"),
#        x = "P.ztest (-log10)",
#        y = "P (-log10)",
#        color = NULL) +
#   theme_light(base_size = 22) +
#   theme(aspect.ratio = 1) +
#   facet_grid(cols = vars(QC_STATUS))
#
# # save plot
# grDevices::png(filename=paste0(base_fig_path, "pz_plot.png"), bg="white", height=600, width=1200, units="px")
# print(plot)
# grDevices::dev.off()
