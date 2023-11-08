#' @title Set dbSNP directory
#' @param path path to the `dbsnp` directory
#' @return NULL, updated config file
#' @export
#' @importFrom yaml read_yaml write_yaml
#'
set_dbsnp_directory <- function(path) {

  stopifnot("`path` must be a valid directory" = dir.exists(path))
  stopifnot("`path` must be a valid directory named `dbsnp`" = basename(path)=="dbsnp")

  config_path <- system.file("config.yml", package="genepi.utils")

  config <- yaml::read_yaml(config_path)

  config[["dbsnp_directory"]] <- path

  yaml::write_yaml(config, config_path)

  message(paste0("dbSNP data directory set to: ", path))
  subdirs <- list.dirs(path, full.names=F)[-1]
  message(paste0(" - available builds:\n", paste0("   - ", subdirs, collapse = "\n")))

  invisible(path)
}

#' @title Get dbSNP directory
#' @return a string file path, the currently set dbSNP directory path
#' @export
#'
which_dbsnp_directory <- function() {

  config_path <- system.file("config.yml", package="genepi.utils")

  config <- yaml::read_yaml(config_path)

  return(config[["dbsnp_directory"]])
}


#' @title Get available dbSNP builds
#' @param build a dbSNP build
#' @return a list of available dbSNP builds - name(dbSNP build): value(directory_path)
#' @export
#'
which_dbsnp_builds <- function(build=NULL) {

  config_path <- system.file("config.yml", package="genepi.utils")

  config <- yaml::read_yaml(config_path)

  build_dirs <- list.dirs(config[["dbsnp_directory"]], full.names=T)[-1]

  names(build_dirs) <- basename(build_dirs)

  if(is.null(build)) {

    return(build_dirs)

  } else {

    if(!build %in% names(build_dirs)) {

      msg <- paste0("Error: `", build, "` is not a valid build name. Valid options: ", paste0(names(build_dirs), collapse = ", "))
      warning(msg)

    } else {

      return(build_dirs[build])

    }
  }
}


#' @title Generate random GWAS data
#' @description
#' Generates rows of synthetic GWAS summary stats data. Useful for developing plotting and other
#' methods. No attempt is made to make this data at all realistic.
#' @param n number of fake variants to generate
#' @param seed seed, for reproducibility
#' @return a data.table with columns SNP, CHR, BP, OA, EA, EAF, BETA, P, EUR_EAF
#' @export
#' @importFrom stats rnorm runif
#'
generate_random_gwas_data <- function(n, seed=2023) {

  EAF = EUR_EAF = SNP = NULL

  set.seed(seed)

  gwas <- data.table::data.table(
    "CHR"= sample(c(as.character(1:22), "X"), size=n, replace=TRUE),
    "BP"= sample(1:200000000, size=n, replace=FALSE),
    "OA" = sample(c("A","C","T","G","AAC","GTTC","TAT"), size=n, replace=TRUE, prob=c(0.95,0.95,0.95,0.95,0.05,0.05,0.05)),
    "EA" = sample(c("A","C","T","G","AAC","GTTC","TAT"), size=n, replace=TRUE, prob=c(0.95,0.95,0.95,0.95,0.05,0.05,0.05)),
    "EAF"= stats::rnorm(n),
    "BETA"= stats::rnorm(n, sd=1.5),
    "P"= stats::runif(n)^2,
    "EUR_EAF"= stats::rnorm(n)
  )

  gwas[, EAF     := (EAF - min(EAF)) / (max(EAF) - min(EAF))]
  gwas[, EUR_EAF := (EUR_EAF - min(EUR_EAF)) / (max(EUR_EAF) - min(EUR_EAF))]
  gwas[, SNP     := paste0(CHR,":",BP,"[b37]",OA,",",EA)]

  return(gwas)
}


