#' #' @title QTL class
#' #' @description
#' #' A short description...
#' #' @param dat a data.table or file path to a data set containing at least columns ...
#' #' @param p_val a numeric, p-value threshold to filter on
#' #' @param join_key a string, column name to join on to other datasets
#' #'
#' #' @return an S3 QTL class data structure
#' #' @export
#' #'
#' QTL <- function(dat, p_val, join_key=NULL) {
#'
#'   dat <- import_table(dat)
#'
#'   structure(
#'     .Data = list(
#'       data = dat,
#'       p_val = p_val,
#'       join_key = join_key),
#'     class = "QTL"
#'   )
#' }
#'
#'
#' #' @title Create a drug target proxy instrument
#' #' @description
#' #' Create a drug target proxy instrument base on several association statistics
#' #' @param gwas_gene .
#' #' @param gene_chr .
#' #' @param gene_start .
#' #' @param gene_end .
#' #' @param gene_flanks .
#' #' @param build .
#' #' @param clump .
#' #' @param clump_ref .
#' #' @param p1 .
#' #' @param p2 .
#' #' @param r2 .
#' #' @param kb .
#' #' @param join_key .
#' #' @param QTL_list .
#' #' @param concordance .
#' #' @return a data.table
#' #' @export
#' #'
#' drug_target_proxy <- function(gwas_gene,
#'                               gene_chr,
#'                               gene_start,
#'                               gene_end,
#'                               gene_flanks = 5e5,
#'                               build = "GRCh37",
#'                               clump = TRUE,
#'                               clump_ref = which_1000G_reference("GRCh37"),
#'                               p1    = 5e-8,
#'                               p2    = 1,
#'                               r2    = 0.2,
#'                               kb    = 250,
#'                               join_key = "RSID",
#'                               QTL_list    = list(),
#'                               concordance = data.table::data.table(data_name_1=character(), data_name_2=character(), concordant=logical())) {
#'
#'   # silence R CMD checks
#'   index = pthresh = data_name_1 = data_name_2 = concordant = clumping = pass = NULL
#'
#'   # checks
#'   stopifnot("`gene_chr` must be a character, in-line with standardise_gwas() CHR column typing" = is.character(gene_chr))
#'   stopifnot("`gene_start`, `gene_end`, `gene_end` and `gene_flanks` must be numeric whole numbers" = all(sapply(c(gene_start,gene_end,gene_end,gene_flanks), function(x) is.numeric(x) & x%%1==0)))
#'   build <- match.arg(build, c("GRCh37", "GRCh38"))
#'
#'   # filter for variants within the gene regions; if clumping extend the region by the clumping window,
#'   # as we can avoids clumping the whole gwas then filtering, instead only clump the region
#'   # -begin-|--clump-win--|--gene-flanks--|-------actual-gene-------|--gene-flanks--|--clump-win--|-end- #
#'   if(clump) {
#'     flank <- gene_flanks + kb
#'   # -begin-|--gene-flanks--|-------actual-gene-------|--gene-flanks--|-end- #
#'   } else {
#'     flank <- gene_flanks
#'   }
#'
#'   # get only variants in the gene region
#'   gene_region <- gwas_gene[CHR==gene_chr &
#'                            BP > (gene_start-flank) &
#'                            BP < (gene_end+flank), ]
#'   if(nrow(gene_region)==0) {
#'     msg <- paste0("After applying gene region filtering [",gene_chr,":",gene_start-flank,"-",gene_end+flank,"] no variants remained")
#'     stop(msg)
#'   }
#'
#'   # run clumping and code result
#'   # clump_outcomes <- c("Not performed"=1, "Index variant"=2, "Clumped variant"=3, "Not in reference"=4, "internal_code_error"=5)
#'   if(clump) {
#'     stopifnot("`plink2` path must be valid - check output of `genepi.utils::which_plink2()` and see the clump vignette for setup details" = file.exists(genepi.utils::which_plink2()))
#'     stopifnot("`plink_ref` path must be valid - check output of `genepi.utils::which_1000G_reference()` and see the clump vignette for setup details" = file.exists(paste0(genepi.utils::which_1000G_reference(build=build), ".pvar")))
#'     gene_region <- genepi.utils::clump(gene_region,
#'                                        p1 = p1,
#'                                        p2 = p2,
#'                                        r2 = r2,
#'                                        kb = kb,
#'                                        plink2 = genepi.utils::which_plink2(),
#'                                        plink_ref = genepi.utils::which_1000G_reference(build))
#'     # code the clumping result
#'     # gene_region[, clumping := data.table::fcase(is.na(index), factor(4, clump_outcomes, names(clump_outcomes)),
#'     #                                             index,        factor(2, clump_outcomes, names(clump_outcomes)),
#'     #                                             !index,       factor(3, clump_outcomes, names(clump_outcomes)),
#'     #                                             default   =   factor(5, clump_outcomes, names(clump_outcomes)))]
#'     gene_region[, clumping := ifelse(index, TRUE, ifelse(!is.na(index), FALSE, NA))]
#'
#'   } else {
#'     # code that we didn't do any clumping
#'     # gene_region[, clumping := factor(1, clump_outcomes, names(clump_outcomes))]
#'     gene_region[, clumping := NA]
#'
#'   }
#'
#'   # run variant thresholding
#'   gene_region[, pthresh := P < p1]
#'
#'   # rename
#'   # data.table::setnames(gene_region, names(gene_region), paste0(names(gene_region), "_gene"))
#'
#'   # use QTL data if provided
#'   if(length(QTL_list)>0) {
#'
#'     # name if not provided
#'     if(is.null(names(QTL_list))) {
#'
#'       QTL_list <- paste0("qtl_", 1:length(QTL_list))
#'
#'     }
#'
#'     # merge the QTL data
#'     data.table::setkeyv(gene_region, join_key)
#'     for(i in seq_along(QTL_list)) {
#'
#'       # the QTL name
#'       qtl_name <- names(QTL_list)[[i]]
#'
#'       # the QTL threshold
#'       qtl_thresh <- QTL_list[[i]]$p_val
#'
#'       # get the QTL data
#'       qtl_data <- QTL_list[[i]]$data
#'
#'       # rename
#'       # data.table::setnames(qtl_data, names(qtl_data), paste0(names(qtl_data), "_", qtl_name))
#'
#'       # set the join key, default is the main gwas_gene join_key
#'       if(!is.null(QTL_list[[i]]$join_key)) {
#'         data.table::setkeyv(qtl_data, QTL_list[[i]]$join_key)
#'       } else {
#'         data.table::setkeyv(qtl_data, paste0(join_key, "_", qtl_name))
#'       }
#'
#'       # merge the datasets and harmonise
#'       gene_region <- harmonise(gene_region, qtl_data, gwas1_trait = "", gwas2_trait = qtl_name)
#'
#'       # code P
#'       gene_region[, paste0("pthresh_", qtl_name) := get(paste0("P_", qtl_name)) < qtl_thresh]
#'
#'     } # end QTL list processing
#'
#'
#'
#'     # assess BETA directions
#'     for(j in 1:nrow(concordance)) {
#'
#'         # output column name
#'         col_name <- ifelse(concordance[j, data_name_1]=="",
#'                            paste0("direction_gwas_", concordance[j, data_name_2]),
#'                            ifelse(concordance[j, data_name_2]=="",
#'                                   paste0("direction_", concordance[j, data_name_1], "_gwas"),
#'                                   paste0("direction_", concordance[j, data_name_1], "_", concordance[j, data_name_2])))
#'
#'         # beta col names
#'         beta_col1 <- ifelse(concordance[j, data_name_1]=="", "BETA", paste0("BETA_", concordance[j, data_name_1]))
#'         beta_col2 <- ifelse(concordance[j, data_name_2]=="", "BETA", paste0("BETA_", concordance[j, data_name_2]))
#'
#'         # code whether the desired direction
#'         gene_region[, (col_name) := ifelse(rep(concordance[j, concordant], nrow(gene_region)),
#'                                            # looking for same direction BETAs
#'                                            sign(get(beta_col1)) == sign(get(beta_col2)),
#'                                            # looking for opposite direction BETAs
#'                                            sign(get(beta_col1)) != sign(get(beta_col2)))]
#'
#'     }
#'
#'   } # end QTL list > 1
#'
#'   # assess overall instrument criteria
#'   criteria_cols <- names(gene_region)[grepl("^(pthresh|direction)", names(gene_region))]
#'   gene_region[, pass := all(.SD),, by=1:nrow(gene_region), .SDcols=criteria_cols]
#'   if(clump) {
#'     gene_region[clumping==FALSE | is.na(clumping), pass := FALSE]
#'   }
#'   data.table::setorder(gene_region, -pass)
#'
#'   # return
#'   return(gene_region)
#' }
#'
#'
#'
#' #' @title Plot the drug proxy instrument
#' #' @param dat a data.table, output from `drug_target_proxy()``
#' #' @param remove a character vector, one or more of c("gwas_pthresh", "clumping")
#' #' @return a ggplot
#' #' @export
#' #'
#' plot_drug_proxy_instrument <- function(dat, remove=c("gwas_pthresh", "clumping")) {
#'
#'   # silence R CMD checks
#'   pthresh = clumping = value = num_passing = variable = RSID = pass = NULL
#'
#'   # checks
#'   remove <- match.arg(remove, c("gwas_pthresh", "clumping"), several.ok=TRUE)
#'
#'   data_cols <- c("P", "BP", "RSID", "SNP", "BETA")
#'   pivot_cols <- c("clumping", "pass",
#'                   names(dat)[grepl("^pthresh",names(dat))],
#'                   names(dat)[grepl("^direction",names(dat))])
#'
#'   if("gwas_pthresh" %in% remove) {
#'     plot_dat <- dat[pthresh==TRUE, ]
#'   }
#'   if("clumping" %in% remove) {
#'     plot_dat <- dat[clumping==TRUE, ]
#'   }
#'
#'   plot_dat <- dat[pthresh==TRUE & clumping==TRUE, ]
#'
#'   plot_dat <- data.table::melt(plot_dat, id.vars=data_cols, measure.vars=pivot_cols)
#'
#'   plot_dat <- stats::na.omit(plot_dat)
#'
#'   ordering <- plot_dat[, list("num_passing" = sum(value)), by="variable" ][order(num_passing), variable]
#'   plot_dat[, variable := factor(variable, levels = c("pass", as.character(ordering[ordering!="pass"])))]
#'
#'   # Selection steps
#'   p_steps <- ggplot2::ggplot(data = plot_dat,
#'                        mapping = ggplot2::aes(x = RSID,
#'                                               y = variable,
#'                                               fill = value)) +
#'     ggplot2::geom_tile(color = "gray") +
#'     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust=1)) +
#'     ggplot2::coord_equal() +
#'     ggplot2::labs(x = "Variant", y = "Selection step", fill = "Pass")
#'
#'
#'   # Manhattan
#'   log10P <- expression(paste("-log"[10], plain(P)))
#'   p_manhattan <- ggplot2::ggplot(data = dat,
#'                                  mapping = ggplot2::aes(x = BP,
#'                                                         y = -log10(P),
#'                                                         color = pass)) +
#'     ggplot2::geom_point() +
#'     ggplot2::geom_abline(slope=0, intercept = -log10(5e-8), linetype="dotted") +
#'     ggplot2::labs(x = "Chromosome position", y = log10P, color = "Pass")
#'
#'   # Combined plot
#'   p <- ggpubr::ggarrange(NULL, p_manhattan, NULL, p_steps, ncol=1,
#'                          heights = c(0.10, 1, 0.05, 1),
#'                          vjust = 2,
#'
#'                          labels = c("Gene region", "", "Variant selection"))
#'
#'   return(p)
#' }
#'
#'
#'
#'
#' #
#' #
#' # # GLP1R significant at <5e8 ?????? P<5e-8 & also clump here?
#' # gwas_glp1r <- gwas_t2dm[P<5e-8 & CHR==glp1r_chr & BP<(glp1r_end_b38+gene_win) & BP>(glp1r_start_b38-gene_win), ][, Phenotype := "t2dm (GLP1R gene)"]
#' # # gwas_glp1r <- clump(gwas_glp1r, p1=5e-8, p2=1, r2=0.2, kb=1)
#' # # gwas_glp1r <- gwas_glp1r[index==TRUE, ]
#' #
#' # # HbA1c variants - build38
#' # # Gill - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10317892/pdf/125_2023_Article_5925.pdf
#' # # GLP1R b37 chr6:39016574–39055519 (they are quoting the Ensembl hg19 in the introduction)
#' # # GLP1R b38 chr6:39048781-39091303
#' # # paper_glp1r_snps <- data.table(rsid = c("rs10305420", "rs34179517", "rs9296291", "rs10305457"),
#' # #                                chr  = rep("6", times=4),
#' # #                                bp_b37 = c(39016636, 39281456, 39056929, 39034095),
#' # #                                bp_b38 = c(39048860, 39313680, 39089153, 39066319))
#' # # paper_glp1r_snps[, within_gene_b37 := paper_glp1r_snps$bp_b37 > (39016574) & paper_glp1r_snps$bp_b37 < (39055519)]
#' # # paper_glp1r_snps[, within_gene_winb37 := paper_glp1r_snps$bp_b37 > (39016574-gene_win) & paper_glp1r_snps$bp_b37 < (39055519+gene_win)]
#' # # paper_glp1r_snps[, within_gene_b38 := paper_glp1r_snps$bp_b38 > glp1r_start_b38 & paper_glp1r_snps$bp_b38 < glp1r_end_b38]
#' # # paper_glp1r_snps[, within_gene_winb38 := paper_glp1r_snps$bp_b38 > (glp1r_start_b38-gene_win) & paper_glp1r_snps$bp_b38 < (glp1r_end_b38+gene_win)]
#' #
#' # # add HbA1c effect - by RSID as build37
#' # #ieugwasr::check_access_token()
#' # #gwas_hba1c <- TwoSampleMR::extract_instruments("ukb-b-19953", p1=5e-08, clump=FALSE, p2=1)
#' # #associations(paste0("6:",glp1r_start_b37,"-",glp1r_end_b37), "ebi-a-GCST90014006", r2=1) |> data.table::as.data.table()
#' # # path_hba1c <- "/Users/xx20081/Downloads/HbA1c_METAL_European.txt"
#' # # path_hba1c <- "/Users/xx20081/Downloads/GCST90025974_buildGRCh37.tsv"
#' # path_hba1c <- "/Users/xx20081/Downloads/MAGIC1000G_HbA1c_EUR.tsv"
#' # gwas_hba1c <- fread(path_hba1c, nThread=12)
#' # gwas_hba1c <- liftover_gwas(gwas_hba1c, to=38, snp_col="variant", chr_col="chromosome", pos_col="base_pair_location", ea_col="effect_allele", oa_col="other_allele")
#' # # HbA1c significant at <0.05 pvalue<0.05 & p_value<0.05
#' # gwas_hba1c <- gwas_hba1c[p_value<0.1 & chromosome==glp1r_chr & base_pair_location<(glp1r_end_b38+gene_win) & base_pair_location>(glp1r_start_b38-gene_win), ]
#' # mapping    <- list(SNP="variant", CHR="chromosome", BP="base_pair_location", EA="effect_allele", OA="other_allele", EAF="effect_allele_frequency", P="p_value", BETA="beta", SE="standard_error", N="sample_size")
#' # gwas_hba1c <- standardise_gwas(gwas_hba1c, input_format=mapping, build="GRCh38", drop=TRUE, populate_rsid="b38_dbsnp156")[, Phenotype := "hba1c"]
#' #
#' # # GTex data v7b38
#' # # gtex_sig_files <- list.files("/Users/xx20081/Downloads/GTEx_Analysis_v8_eQTL", pattern = "v8.signif_variant_gene_pairs.txt.gz", full.names = T)
#' # # gtex_data <- lapply(gtex_sig_files, function(x) {
#' # #     d <- fread(x)
#' # #     d[, c("CHR","BP","OA","EA","build") := tstrsplit(variant_id, "_", fixed=TRUE)]
#' # #     d[, CHR := sub("chr","",CHR, ignore.case = TRUE)]
#' # #     d[, BP := as.numeric(BP)]
#' # #     d <- d[CHR==glp1r_chr & BP < (glp1r_end_b38+gene_win) & BP > (glp1r_start_b38-gene_win), ]
#' # #     if(nrow(d)>0) {
#' # #       mapping <- list(SNP="variant_id",CHR="CHR",BP="BP",OA="OA",EA="EA",EAF="maf",BETA="slope",SE="slope_se",P="pval_nominal")
#' # #       d <- standardise_gwas(d, mapping, drop=TRUE, build="GRCh38", populate_rsid="b38_dbsnp156")
#' # #       return(d)
#' # #     }else{
#' # #       return(NULL)
#' # #     }
#' # #   }) |>
#' # #   `names<-`(sub(".v8.signif_variant_gene_pairs.txt.gz","",basename(gtex_sig_files))) |>
#' # #   rbindlist(idcol="tissue", fill=TRUE)
#' # # fwrite(gtex_data, "/Users/xx20081/Downloads/tmp_gtexv8_data.tsv", sep="\t")
#' # gtex_data <- fread("/Users/xx20081/Downloads/tmp_gtexv8_data.tsv", nThread=12)
#' #
#' # # merge HbA1c and gtex to find likely variants
#' # # directionally concordant associations with gene expression if
#' # # - lower HbA1c (negative beta) and higher GLP1R (positive slope), or
#' # # - higher HbA1c (positve beta) and lower GLP1R (negative slope)
#' # # concordant Gtex and HbA1c results
#' # gwas_hba1c[, CHR := as.character(CHR)]
#' # gtex_data[, CHR := as.character(CHR)]
#' # snps_sig_and_concord_effects <- genepi.utils::harmonise(gwas_hba1c, gtex_data, "hba1c", "gtex")
#' # snps_sig_and_concord_effects <- snps_sig_and_concord_effects[((BETA_gtex>0 & BETA_hba1c<0) | (BETA_gtex<0 & BETA_hba1c>0)), ]
#' #
#' # # merge likely variants with T2DM/GLP1R variants - these are the instrument
#' # data.table::setkey(gwas_glp1r, RSID)
#' # data.table::setkey(snps_sig_and_concord_effects, RSID_hba1c)
#' # merge_data <- gwas_glp1r[snps_sig_and_concord_effects, nomatch = NULL]
#' # merge_data <- unique(merge_data, by="RSID")
#' #

