#' set_dbsnp_directory
#'
#' @param path path to the `dbsnp` directory
#'
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
