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

#' which_dbsnp_directory
#'
#' @return a string file path, the currently set dbSNP directory path
#' @export
#'
which_dbsnp_directory <- function() {

  config_path <- system.file("config.yml", package="genepi.utils")

  config <- yaml::read_yaml(config_path)

  return(config[["dbsnp_directory"]])
}


#' which_dbsnp_builds
#'
#' @param build a dbSNP build
#'
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


