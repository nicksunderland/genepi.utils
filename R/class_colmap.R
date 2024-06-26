#' @title Column object
#' @param name the standard column name
#' @param alias a character vector of aliases (other column names) for this column
#' @param type a character, an atomic R type
#' @slot name the standard column name
#' @slot alias a character vector of aliases (other column names) for this column
#' @slot type a character, an atomic R type
#' @return an S7 class genepi.utils::Column object
#' @export
Column <- new_class(
  #==============================
  # Column class name
  #==============================
  name       = "Column",
  package    = "genepi.utils",
  properties = list(
    name  = class_character,
    alias = class_character,
    type  = class_character
  ),
  validator = function(self) {
    std_cols <- c('rsid','chr','bp','ea','oa','eaf','p','beta','se','or','ncase', 'info','or_se','or_lb','or_ub',
                  'beta_lb','beta_ub','z','q_stat', 'q_p', 'strand', 'imputed', 'i2','nstudies','n','effects')
    if (!self@name %in% std_cols) {
      msg <- paste0("standard column name not recognised `", self@name, "` \nOptions: ", paste0(std_cols, collapse = ", "))
      stop(msg)
    }
    stopifnot("standard column type not recognised; use atomic R types" = self@type %in% c('character','integer','numeric','logical'))
  }
)

#' @title ColumnMap object
#' @description
#' A mapping to the standardised column names used in this package. Available names: 'rsid', 'chr', 'bp', 'ea', 'oa', 'eaf', 'p', 'beta',
#' 'se', 'or', 'or_se', 'or_lb', 'or_ub', 'beta_lb', 'beta_ub', 'z', 'q_stat', 'i2', 'nstudies', 'n'
#' @param x either a list of `Column` class objects, a valid string for a pre-defined map: default, metal, ieu_ukb, ieugwasr, ns_map,
#' gwama, giant, or a named character vector or list (standard name = old name)
#' @slot map a list of `Column` class objects
#' @return an S7 class genepi.utils::ColumnMap object
#' @export
ColumnMap <- new_class(
  #==============================
  # Column class name
  #==============================
  name       = "ColumnMap",
  package    = "genepi.utils",
  properties = list(
    map   = class_list
  ),
  constructor = function(x) {

    # passing a genepi.utils::ColumnMap object
    if(inherits(x, 'genepi.utils::ColumnMap')) {

      map <- x@map

    # passing a list of genepi.utils::Column objects
    } else if(all(sapply(x, function(x0) inherits(x0, 'genepi.utils::Column')))) {

      map <- x

    # passing single/list/vector of column names
    } else {

      # possible column names and their meta-data
      base_map = list(
        rsid    = Column(name='rsid',    type='character', alias=c('rsid','RSID','SNP','MarkerName','MARKER','rs_number')),
        chr     = Column(name='chr',     type='character', alias=c('chr','CHR','chromosome')),
        bp      = Column(name='bp',      type='integer',   alias=c('bp','BP','POS','position')),
        ea      = Column(name='ea',      type='character', alias=c('ea','EA','A1','Allele1','ALLELE1','Tested_Allele','reference_allele')),
        oa      = Column(name='oa',      type='character', alias=c('oa','OA','A2','Allele2','ALLELE0','Other_Allele','nea','other_allele')),
        eaf     = Column(name='eaf',     type='numeric',   alias=c('eaf','EAF','Freq1','A1FREQ','Freq_Tested_Allele')),
        beta    = Column(name='beta',    type='numeric',   alias=c('beta','BETA','Effect')),
        se      = Column(name='se',      type='numeric',   alias=c('se','SE','StdErr')),
        beta_lb = Column(name='beta_lb', type='numeric',   alias=c('beta_lb','beta_95L')),
        beta_ub = Column(name='beta_ub', type='numeric',   alias=c('beta_ub','beta_95U')),
        z       = Column(name='z',       type='numeric',   alias=c('z','Z')),
        p       = Column(name='p',       type='numeric',   alias=c('p','P','P-value','p-value','P_BOLT_LMM_INF')),
        n       = Column(name='n',       type='integer',   alias=c('n','N','n_samples')),
        ncase   = Column(name='ncase',   type='integer',   alias=c('ncase','NCASE','N_CASE')),
        nstudies= Column(name='nstudies',type='integer',   alias=c('nstudies','n_studies')),
        or      = Column(name='or',      type='numeric',   alias=c('or','OR')),
        or_se   = Column(name='or_se',   type='numeric',   alias=c('or_se','OR_SE')),
        or_lb   = Column(name='or_lb',   type='numeric',   alias=c('or_lb','OR_LB')),
        or_ub   = Column(name='or_ub',   type='numeric',   alias=c('or_ub','OR_UB')),
        q_stat  = Column(name='q_stat',  type='numeric',   alias=c('q_stat','q_statistic')),
        q_p     = Column(name='q_p',     type='numeric',   alias=c('q_p','q_p-value')),
        i2      = Column(name='i2',      type='numeric',   alias=c('i2')),
        info    = Column(name='info',    type='numeric',   alias=c('info','INFO')),
        imputed = Column(name='imputed', type='logical',   alias=c('imputed','Imputed')),
        strand  = Column(name='strand',  type='character', alias=c('strand','STRAND','Strand','Orientation')),
        effects = Column(name='effects', type='character', alias=c('effects'))
      )

      # pre-specified maps
      defined_maps <- list(
        default  = c('RSID','CHR','BP','EA','OA','EAF','P','BETA','SE','OR','OR_SE','OR_LB','OR_UB'),
        metal    = c('MarkerName','Allele1','Allele2','Freq1','P-value','Effect','StdErr'),
        ieu_ukb  = c('SNP','ALLELE1','ALLELE0','A1FREQ','P_BOLT_LMM_INF','BETA','SE'),
        ieugwasr = c('rsid','chr','position','ea','nea','eaf','p','beta','se','n'),
        ns_map   = c('MARKER','CHR','POS','A1','A2','EAF','P','BETA','SE'),
        gwama    = c('rs_number','chromosome','position','reference_allele','other_allele','eaf','p-value','beta','se','BETA_95L','beta_95U','z','q_statistic','q_p-value','i2','n_studies','n_samples','effects'),
        giant    = c('SNP','CHR','POS','Tested_Allele','Other_Allele','Freq_Tested_Allele','P','BETA','SE','N','INFO')
      )

      # passing a single string
      if(length(x)==1 && x %in% names(defined_maps)) {

        x <- defined_maps[[x]]

      }

      # passing a named list or named character vector
      if(!is.null(names(x))) {

        stopifnot("All names must be standard columns" = all(names(x) %in% names(base_map)))
        map <- base_map[names(x)]
        map <- Map(function(map, name) {
          map@alias <- unique(c(name, map@alias)) # name to the front
          return(map)
        }, map, unname(x))

      # passing a list or character vector, check aliases
      } else {

        map <- lapply(x, function(x0) {
          i <- which(sapply(base_map, function(cm) x0 %in% cm@alias))
          if(length(i) > 0) {
            map <- base_map[[i]]
            map@alias <- unique(c(x0, map@alias)) # name to the front
            return(map)
          }
        })
        map <- map[!sapply(map, is.null)]

      }

    }

    # assign to object
    names(map) <- sapply(map, function(n) n@name)
    new_object(S7_object(), map=map)
  },
  validator = function(self) {
    stopifnot("@map must be a list of class genepi.utils::Column elements" = all(sapply(self@map, function(x) inherits(x, 'genepi.utils::Column'))))
  }
)


raw_names <- new_generic('raw_names', 'map')
method(raw_names, ColumnMap) <- function(map, dat_cols) {

  sapply(lapply(map@map, function(x) x@alias), function(y) y[which(y %in% dat_cols)[1]])

}

method(`[[`, ColumnMap) <- function(map, i) { map@map[[i]] }
method(`$`, ColumnMap) <- function(map, i) { map@map[[i]] }
method(`[`, ColumnMap) <- function(map, i) {
  map@map <- map@map[i]
  map
}

