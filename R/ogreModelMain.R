# get_lab_values ---------------------------------------------------------------

#' Reads all lab values from ODBC-source
#'
#' Appends also site code (e.g., "NEU") and substance name apart from
#' standard-fields the fields "CensorCode", "QualityControlLevelID", and
#' "UnitsAbbreviation" are included.  If one measurements exists twice, only
#' higher QualityControlLevelID is included
#'
#' @param odbc_name Name of the odbc source
#'
#' @export
#'
#' @importFrom kwb.odm odm_DataValues
#' @importFrom kwb.odm odm_Samples
#' @importFrom kwb.odm odm_Sites
#' @importFrom kwb.odm odm_Units
#' @importFrom kwb.odm odm_Variables
#'
get_lab_values <- function(odbc_name)
{
  # Function returning a modified version of a kwb.odm-function. In the returned
  # version, the function kwb.db::selectFromDb() will be called with
  # "stringsAsFactors = FALSE".
  no_factor_function <- function(FUN) {
    function(...) FUN(..., stringsAsFactors = FALSE)
  }

  get_sites <- no_factor_function(kwb.odm::odm_Sites)
  get_variables <- no_factor_function(kwb.odm::odm_Variables)
  get_samples <- no_factor_function(kwb.odm::odm_Samples)
  get_values <- no_factor_function(kwb.odm::odm_DataValues)
  get_units <- no_factor_function(kwb.odm::odm_Units)

  # get info from db
  sites <- get_sites(db = odbc_name, select = 1:3)

  substances <- get_variables(db = odbc_name, select = c(1:3, 5))

  samples <- get_samples(db = odbc_name, select = 1:2)

  values <- get_values(
    db = odbc_name, select = paste0(
      "ValueID, DataValue, LocalDateTime, DateTimeUTC, UTCOffset, ",
      "SampleID, SiteID, VariableID, CensorCode, QualityControlLevelID"
    ),
    orderBy_QualityControlLevelID = 1,
    as.is = TRUE # keep dates and times as character
  )

  # Explicitly convert character to POSIXct
  values$LocalDateTime <- as.POSIXct(values$LocalDateTime, tz = "Europe/Berlin")
  values$DateTimeUTC <- as.POSIXct(values$DateTimeUTC, tz = "UTC")

  units_variables <- get_units(db = odbc_name, select = c(1, 4))

  # merge info in one data.frame

  x_all <- merge(values, sites, by = "SiteID")

  x_all <- merge(x_all, substances, by = "VariableID")

  x_all <- merge(
    x_all, units_variables, by.x = "VariableUnitsID", by.y = "UnitsID"
  )

  x_all <- merge(x_all, samples, by = "SampleID")

  # as merged dataframe x_all is not sorted by QualityControlLevelID anymore, do it again
  x_all <- x_all[order(x_all$QualityControlLevelID), ]

  # keep only higher QualityControlLevelID if measurements are double
  y <- stats::aggregate(
    ValueID ~ SampleID + VariableID, data = x_all, FUN = length
  )

  y2 <- stats::aggregate(
    ValueID ~ SampleID + VariableID, data = x_all, FUN = "[", 1
  )

  ValueIDs_LQ <- y2$ValueID[which(y[, 3] > 1)]

  indices_LQ <- match(ValueIDs_LQ, x_all$ValueID)

  if (length(indices_LQ) > 0) {
    x_all <- x_all[- indices_LQ, ]
  }

  x_all
}

# only_new_dl_metals -----------------------------------------------------------
#' removes metal samples below detection limit (dl),
#'
#' when dl was too high (old analytical method). Works by VariableCode
#'
#' @param x_in name of input data.frame
#'
#' @export
only_new_dl_metals <- function(x_in)
{
  # get detection limits that were lowered during project
  DL_metals_new <- getNewDetectionLimits()

  # delete samples below detection limit, where detection limit is higher than
  # current
  for (i in seq_along(DL_metals_new$VariableCode)) {

    indices <- which(
      x_in$VariableCode == DL_metals_new[i, 1] &
      x_in$CensorCode == "lt" &
      x_in$DataValue > DL_metals_new[i, 2]
    )

    if (length(indices) > 0) {
      x_in <- x_in[- indices, ]
    }
  }

  x_in
}

# only_composite ---------------------------------------------------------------
#' removes rows in data.frame with SampleType ! = "Composite"
#'
#' @param x_in name of input data.frame
#'
#' @export
only_composite <- function(x_in)
{
  x_in[which(x_in$SampleType == "Composite"), ]
}

# remove_group -----------------------------------------------------------------
#' removes measurements of Variables of a specific group
#'
#' @param x_in name of input data.frame
#' @param group Variable group to be removed (string)
#'
#' @export
remove_group <- function(x_in, group)
{
  variables <- kwb.ogre::OGRE_VARIABLES()

  codes <- variables$VariableCode[variables$VariableGroupName == paste(group)]

  x_in[which(x_in$VariableCode %in% codes == "FALSE"), ]
}

# no_Panke ---------------------------------------------------------------------
#' removes rows in data.frame with site code = "PNK"
#'
#' @param x_in name of input data.frame
#'
#' @export
no_Panke <- function(x_in)
{
  x_in[which(x_in$SiteCode != "PNK"), ]
}

# Panke ------------------------------------------------------------------------
#' keeps only rows in data.frame with site code = "PNK"
#'
#' @param x_in name of input data.frame
#'
#' @export
Panke <- function(x_in)
{
  x_in[which(x_in$SiteCode == "PNK"), ]
}

# reduce_codes -----------------------------------------------------------------
#' removes lines with censor codes
#'
#' other than "lt" and "nc" (e.g., "???" are removed)
#'
#' @param x_in name of input data.frame
#'
#' @export
reduce_codes <- function(x_in)
{
  x_in[which(x_in$CensorCode %in% c("lt", "nc")), ]
}

# only_representative_subst ----------------------------------------------------
#' removes Variables, which have not
#'
#' at least one measurement (can also be below detection limit) per catchment
#' type
#'
#' @param x_in name of input data.frame
#'
#' @export
only_representative_subst <- function(x_in)
{
  VariableIDs <- unique(x_in$VariableID)
  siteIDs <- unique(x_in$SiteID)

  for (i in seq_along(VariableIDs)) {

    indices <- which(x_in$VariableID == VariableIDs[i])
    site_match <- siteIDs %in% x_in$SiteID[indices]

    if (sum(site_match) < length(siteIDs)) {

      x_in <- x_in[- indices, ]
    }
  }

  x_in
}

# non_detect -------------------------------------------------------------------
#' lists substances without detection in any sample
#'
#' Apart from entire dataset x_in (first column), lists are given for each site
#' individually (following columns)
#'
#' @param x_in name of input data.frame
#'
#' @export
non_detect <- function(x_in)
{
  site_IDs <- unique(x_in$SiteID)
  y <- order(site_IDs)
  site_IDs <- site_IDs[y]

  mat_nodetect <- matrix(nrow = 200, ncol = length(site_IDs) + 1)

  #find substances which are < dl in all samples
  indices_lt <- which(x_in$CensorCode == "lt")
  subst_lt <- unique(x_in$VariableName[indices_lt])
  subst_nc <- unique(x_in$VariableName[- indices_lt])
  subst_lt <- setdiff(subst_lt,subst_nc)
  mat_nodetect[seq_along(subst_lt), 1] <- subst_lt
  max_length <- length(subst_lt)

  #find substances which are < dl for all samples of one site

  for (i in seq_along(site_IDs)) {

    indices_lt <- which(x_in$CensorCode == "lt" & x_in$SiteID == site_IDs[i])
    indices_nc <- which(x_in$CensorCode == "nc" & x_in$SiteID == site_IDs[i])
    subst_lt <- unique(x_in$VariableName[indices_lt])
    subst_nc <- unique(x_in$VariableName[indices_nc])
    subst_lt <- setdiff(subst_lt, subst_nc)
    mat_nodetect[seq_along(subst_lt), i + 1] <- subst_lt
    max_length <- max(length(subst_lt), max_length)
  }

  mat_nodetect <- mat_nodetect[seq_len(max_length), ]

  non_detected <- as.data.frame(mat_nodetect, stringsAsFactors = FALSE)

  site_names <- x_in$SiteCode[match(site_IDs, x_in$SiteID)]

  stats::setNames(non_detected, c("All_sites", site_names))
}

# detect -----------------------------------------------------------------------
#' removes substances, without detection from data.frame
#'
#' requires list of these substances as single vector
#'
#' @param x_in name of input data.frame
#' @param x_nd vector with substances < detection
#'
#'
#' @export
detect <- function(x_in, x_nd)
{
  x_in[which(is.na(match(x_in$VariableName, x_nd))), ]
}

# adapt_nondetect --------------------------------------------------------------
#' adapts values of single results < detection limit.
#'
#' requires list of these substances as data.frame for each site (one column per
#' site). If substance is always < dl at one site, results are set to zero. If
#' substance is sometimes > dl, resluts are set to a factor*dl
#'
#' @param x_in name of input data.frame
#' @param x_nd vector with substances < detection (one column per site)
#' @param factor multiplier of detection limit if smaller dl (e.g., 0, 0.5 or 1)
#'
#' @export
adapt_nondetect <- function(x_in, x_nd, factor = 0.5)
{
  site_names <- colnames(x_nd)[-1]

  for (i in seq_len(dim(x_nd)[2] - 1)) {

    #set samples <dl for substances which are never detected = 0
    y <- match(x_in$VariableName, x_nd[, i + 1])
    indices <- which(x_in$SiteCode == site_names[i] & y > 0)
    x_in$DataValue[indices] <- 0
  }

  #set samples <dl for substances which are partially detected = factor*dl
  indices <- which(x_in$CensorCode == "lt" & x_in$DataValue != 0)
  x_in$DataValue[indices] <- factor * x_in$DataValue[indices]

  x_in
}

# annual_mean_conc -------------------------------------------------------------
#' estimates annual mean concentrations per site.
#'
#' Apart from a matrix with mean concentrations for each substance and site (= 0
#' if always below detection limit), matrices with N, standard error, standard
#' deviation, RMSE, as well as measured min and max are calculated. Result is
#' given as a list, of these matrices. Different methods can be chosen:
#'  Method 1: arithmetic mean
#'  Method 2: functions and RMSE from file, arithmetic mean for substances
#'  without functions
#'
#' @param x_in name of input data.frame
#' @param method estimation method: 1 = arithmetic mean; 2 = functions and RMSE
#'   from file, arithmetic mean for substances without functions
#' @param data.dir file directory where correlation data and rain series for
#'   method 2 are located
#'
#' @export
annual_mean_conc <- function(x_in, method, data.dir)
{
  #order by VariableID
  y <- order(x_in$VariableID)
  x_in <- x_in[y, ]

  #find existing sites
  site_IDs <- unique(x_in$SiteID)
  y <- order(site_IDs)
  site_IDs <- site_IDs[y]
  indices <- match(site_IDs, x_in$SiteID)
  site_names <- x_in$SiteCode[indices]

  #output structure
  x_out_mean <- data.frame(
    VariableID = unique(x_in$VariableID),
    VariableName = unique(x_in$VariableName),
    stringsAsFactors = FALSE
  )

  indices <- match(x_out_mean$VariableID, x_in$VariableID)
  x_out_mean$UnitsAbbreviation <- x_in$UnitsAbbreviation[indices]

  x_out_stdev <- x_out_mean

  if (method == 1) {

    #apply method 1 (arithmetic mean for all sites)
    x_out <- annual_mean_conc_method1(
      x_in = x_in,
      x_out_mean = x_out_mean,
      x_out_stdev = x_out_stdev,
      site_names = site_names
    )

  } else if (method == 2) {

    #apply method 1 (arithmetic mean for all sites), 1st step
    x_out_method1 <- annual_mean_conc_method1(
      x_in = x_in,
      x_out_mean = x_out_mean,
      x_out_stdev = x_out_stdev,
      site_names = site_names
    )

    #apply method 2 (replace method 1 where correlations exist)

    x_out <- annual_mean_conc_method2(
      x_out_mean = x_out_method1$mean,
      x_out_stdev = x_out_method1$stdev,
      site_names = site_names,
      data.dir = data.dir
    )

  } else {

    stop("method is not established yet")
  }

  x_out
}

# default_statistics -----------------------------------------------------------
#' calculate default statistics
#'
#' calculate default statistics for a grouped data frame (created with
#' dplyr::group_by)
#'
#' @param x data frame with columns \code{DataValue}, \code{CensorCode}
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr summarise
#' @importFrom rlang .data
#' @export
default_statistics <- function(x)
{
  x %>% summarise(
    mean = mean(.data$DataValue),
    N = length(.data$DataValue),
    N_lt = sum(.data$CensorCode == "lt"),
    N_nc = sum(.data$CensorCode == "nc"),
    stdev = stats::sd(.data$DataValue),
    var = stats::var(.data$DataValue),
    se = stats::sd(.data$DataValue) / sqrt(length(.data$DataValue) - 1),
    min = min(.data$DataValue),
    max = max(.data$DataValue),
    median = stats::median(.data$DataValue),
    quantile25 = quant25(.data$DataValue),
    quantile75 = quant75(.data$DataValue),
    quantile95 = quant95(.data$DataValue),
    mean_geom = geom_mean(.data$DataValue)
  )
}

# annual_stats -----------------------------------------------------------------
#' estimates annual mean concentrations per site.
#'
#' Apart from a matrix with mean concentrations for each substance and site (= 0
#' if always below detection limit), matrices with N, 95% confidence, as well as
#' measured min and max are calculated. Result is given as a list, of these
#' matrices.
#'
#' @param x_in name of input data.frame
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom rlang .data
#' @export
annual_stats <- function(x_in)
{
  # Provide vector of SiteCodes ordered by their SiteID
  site_names <- (
    x_in %>%
      dplyr::group_by(
        .data$SiteID,
        .data$SiteCode
      ) %>%
      summarise()
  )$SiteCode

  # Provide statistics grouped by SiteCode and Variable
  x_summary <- x_in %>%
    dplyr::group_by(
      .data$SiteCode,
      .data$VariableID,
      .data$VariableName,
      .data$UnitsAbbreviation
    ) %>%
    default_statistics()

  # Provide statistics grouped by Variable only
  x_total <- x_in %>%
    dplyr::group_by(
      .data$VariableID,
      .data$VariableName,
      .data$UnitsAbbreviation
    ) %>%
    default_statistics()

  # Provide vectors of column names
  columns.variable <- c("VariableID", "VariableName", "UnitsAbbreviation")
  columns.key <- c("SiteCode", columns.variable)
  columns.stat <- setdiff(names(x_summary), columns.key)

  x_out <- list()

  for (column in columns.stat) {

    # Reshape from "long" to "wide"
    result <- stats::reshape(
      data = as.data.frame(x_summary[, c(columns.key, column)]),
      direction = "wide",
      idvar = columns.variable,
      timevar = "SiteCode"
    )

    # Remove attribute given by reshape
    attr(result, "reshapeWide") <- NULL

    # Remove column name prefixes given by rehape
    names(result) <- gsub(paste0("^", column, "\\."), "", names(result))

    # Merge site statistics with totals
    x <- result[, c(columns.variable, site_names)]

    y <- x_total[, c(columns.variable, column)]
    names(y)[ncol(y)] <- "Gesamt"

    result <- merge(x, y) %>% arrange(.data$VariableID)

    x_out[[column]] <- result
  }

  x_out
}

# annual_load_rain -------------------------------------------------------------
#' calculates the load for each substance
#'
#' separates pathways (rain runoff, CSO and WWTP)
#'
#' @param data.dir path of model data (annual mean concentrations
#'   "annual_mean_conc.csv",  rain runoff volumes "Vol_rain.csv", removal at
#'   WWTP "substance_info.csv")
#' @param error_removal_rate relative error in removal at WWTP
#'
#' @return Function returns list with loads and standard deviations,
#' by entry path (sep, cso, wwtp) and by surface water catchment.
#' Concentration in units "mg/L" and "ug/L" is automatically
#' transformed to loads in "kg/yr". Other (unknown) units are
#' left unchanged, resulting in "unit * m3/yr".
#'
#' @export
annual_load_rain <- function(data.dir, error_removal_rate = 0.3)
{
  #load data
  vol_rain <- read_1(file = get_path_or_stop(
    data.dir, "Vol_rain.csv", "rain runoff"
  ))

  error_vol_rain <- read_2(file = get_path_or_stop(
    data.dir, "Vol_rain_relative_error.csv", "rain runoff"
  ))

  x_conc <- read_2(file = get_path_or_stop(
    data.dir, "annual_mean_conc.csv", "annual mean concentrations"
  ))

  error_x_conc <- read_2(file = get_path_or_stop(
    data.dir, "annual_mean_conc_relative_error.csv",
    "annual mean concentrations"
  ))

  removal_rates <- read_1(file = get_path_or_stop(
    data.dir, "substance_info.csv", "removal rates at WWTP"
  ))

  # missing removal rates are set = 0
  removal_rates[, 2] <- as.numeric(removal_rates[, 2])
  removal_rates[which(is.na(removal_rates$Retention_.)), 2] <- 0

  # get removal rates for substances in x_conc only (and in same order)
  removal_rates_red <- x_conc[, 1:2]
  indices <- match(removal_rates_red$VariableName, removal_rates$VariableName)
  removal_rates_red$removal_percent <- removal_rates$Retention_.[indices]

  # order of catchments
  sum_EZG <- stats::aggregate(
    vol_rain$GESAMT, by = list(vol_rain$SUW), FUN = "sum"
  )

  indices <- order(sum_EZG$x, decreasing = TRUE)
  sum_EZG <- sum_EZG[indices, ]
  EZG_names <- sum_EZG[, 1]

  # structure of calculation files
  load_sep <- x_conc[, c(1, 2, 4:9)]
  load_sep[, 3:8] <- 0
  load_sep$TOT <- 0
  error_load_sep <- load_sep
  load_CSO <- load_sep
  error_load_CSO <- load_sep
  load_WWTP <- load_sep
  error_load_WWTP <- load_sep

  # output list
  x_out_by_catchment_kg_yr <- list()
  error_by_catchment_kg_yr <- list()
  x_out_by_pathway_kg_yr <- list()
  error_by_pathway_kg_yr <- list()

  for (EZG in EZG_names) {

    indices <- which(vol_rain$SUW == EZG)

    vol_rain_EZG <- vol_rain[indices, ]
    error_vol_rain_EZG <- error_vol_rain[indices, ]

    # Define vector of site codes
    sites <- c("ALT", "NEU", "STR", "EFH", "GEW", "ANDERE")

    # loads of rain-water based substances via separate sewer system
    for (site in sites) {

      load_sep[[site]] <- x_conc[[site]] * vol_rain_EZG[[site]][1]
    }

    load_sep$TOT <-
      load_sep$ALT +
      load_sep$NEU +
      load_sep$STR +
      load_sep$EFH +
      load_sep$GEW +
      load_sep$ANDERE

    # absolute errors in loads of rain-water based substances via separate sewer system
    for (site in sites) {

      error_load_sep[[site]] <- load_sep[[site]] * sqrt(
        error_x_conc[[site]]^2 + error_vol_rain_EZG[[site]][1]^2
      )
    }

    error_load_sep$TOT <- sqrt(
      error_load_sep$ALT^2 +
      error_load_sep$NEU^2 +
      error_load_sep$STR^2 +
      error_load_sep$EFH^2 +
      error_load_sep$GEW^2 +
      error_load_sep$ANDERE^2
    )

    # loads of rain-water based substances via CSO
    for (site in sites) {

      load_CSO[[site]] <- x_conc[[site]] * vol_rain_EZG[[site]][2]
    }

    load_CSO$TOT <-
      load_CSO$ALT +
      load_CSO$NEU +
      load_CSO$STR +
      load_CSO$EFH +
      load_CSO$GEW +
      load_CSO$ANDERE

    # absolute errors in loads of rain-water based substances via CSO
    for (site in sites) {

      error_load_CSO[[site]] <- load_CSO[[site]] * sqrt(
        error_x_conc[[site]]^2 + error_vol_rain_EZG[[site]][2]^2
      )
    }

    error_load_CSO$TOT <- sqrt(
      error_load_CSO$ALT^2 +
      error_load_CSO$NEU^2 +
      error_load_CSO$STR^2 +
      error_load_CSO$EFH^2 +
      error_load_CSO$GEW^2 +
      error_load_CSO$ANDERE^2
    )

    # loads of rain-water based substances via WWTP
    rate_remaining <- (1 - removal_rates_red$removal_percent / 100)

    for (site in sites) {

      load_WWTP[[site]] <-
        x_conc[[site]] *
        vol_rain_EZG[[site]][3] *
        rate_remaining
    }

    load_WWTP$TOT <-
      load_WWTP$ALT +
      load_WWTP$NEU +
      load_WWTP$STR +
      load_WWTP$EFH +
      load_WWTP$GEW +
      load_WWTP$ANDERE

    # error in loads of rain-water based substances via WWTP
    for (site in sites) {

      error_load_WWTP[[site]] <- load_WWTP[[site]] * sqrt(
        error_x_conc[[site]]^2 +
        error_vol_rain_EZG[[site]][3]^2 +
        error_removal_rate^2
      )
    }

    error_load_WWTP$TOT <- sqrt(
      error_load_WWTP$ALT^2 +
      error_load_WWTP$NEU^2 +
      error_load_WWTP$STR^2 +
      error_load_WWTP$EFH^2 +
      error_load_WWTP$GEW^2 +
      error_load_WWTP$ANDERE^2
    )

    # load unit for all substances in kg/yr (where concentration unit is known)
    indices_mgL <- which(x_conc$UnitsAbbreviation == "mg/L")
    indices_ugL <- which(x_conc$UnitsAbbreviation == "ug/L")

    skip2 <- -(1:2)

    load_sep[indices_mgL, skip2] <- load_sep[indices_mgL, skip2] / 1000
    load_sep[indices_ugL, skip2] <- load_sep[indices_ugL, skip2] / 1e6

    load_CSO[indices_mgL, skip2] <- load_CSO[indices_mgL, skip2] / 1000
    load_CSO[indices_ugL, skip2] <- load_CSO[indices_ugL, skip2] / 1e6

    load_WWTP[indices_mgL, skip2] <- load_WWTP[indices_mgL, skip2] / 1000
    load_WWTP[indices_ugL, skip2] <- load_WWTP[indices_ugL, skip2] / 1e6

    error_load_sep[indices_mgL, skip2] <- error_load_sep[indices_mgL, skip2] / 1000
    error_load_sep[indices_ugL, skip2] <- error_load_sep[indices_ugL, skip2] / 1e6

    error_load_CSO[indices_mgL, skip2] <- error_load_CSO[indices_mgL, skip2] / 1000
    error_load_CSO[indices_ugL, skip2] <- error_load_CSO[indices_ugL, skip2] / 1e6

    error_load_WWTP[indices_mgL, skip2] <- error_load_WWTP[indices_mgL, skip2] / 1000
    error_load_WWTP[indices_ugL, skip2] <- error_load_WWTP[indices_ugL, skip2] / 1e6

    # summary by pathway
    load_summary_path <- x_conc[, 1:2]
    load_summary_path$sep <- load_sep$TOT
    load_summary_path$CSO <- load_CSO$TOT
    load_summary_path$WWTP <- load_WWTP$TOT
    load_summary_path$TOT <- load_sep$TOT + load_CSO$TOT + load_WWTP$TOT

    #absolute error by pathway
    error_load_summary_path <- load_summary_path
    error_load_summary_path$sep <- error_load_sep$TOT
    error_load_summary_path$CSO <- error_load_CSO$TOT
    error_load_summary_path$WWTP <- error_load_WWTP$TOT

    error_load_summary_path$TOT <- sqrt(
      error_load_sep$TOT^2 +
      error_load_CSO$TOT^2 +
      error_load_WWTP$TOT^2
    )

    # summary by source catchment
    load_summary_catch <- load_sep
    j <- 3:9
    load_summary_catch[, j] <- load_sep[, j] + load_CSO[, j] + load_WWTP[, j]

    # absolute error by source catchment
    error_load_summary_catch <- error_load_sep
    error_load_summary_catch[, j] <- sqrt(
      error_load_sep[, j]^2 +
      error_load_CSO[, j]^2 +
      error_load_WWTP[, j]^2
    )

    # organize in list
    x_out_by_catchment_kg_yr[[EZG]] <- load_summary_catch
    x_out_by_pathway_kg_yr[[EZG]] <- load_summary_path
    error_by_catchment_kg_yr[[EZG]] <- error_load_summary_catch
    error_by_pathway_kg_yr[[EZG]] <- error_load_summary_path
  }

  # output
  list(
    by_path = x_out_by_pathway_kg_yr,
    error_by_path = error_by_pathway_kg_yr,
    by_catch = x_out_by_catchment_kg_yr,
    error_by_catch = error_by_catchment_kg_yr
  )
}

# get_path_or_stop -------------------------------------------------------------
get_path_or_stop <- function(data_dir, file_name, subject)
{
  file <- file.path(data_dir, file_name)

  if (! file.exists(file)) stop(
    sprintf("File with %s (%s) not found in data.dir", subject, file_name),
    call. = FALSE
  )

  file
}

# annual_load_sewage -----------------------------------------------------------
#' calculates the load for each substance
#'
#' separates pathways (CSO and WWTP)
#'
#' @param data.dir path of model data (annual mean concentrations
#'   "substance_info.csv",  WWTP runoff volumes "Vol_sewage.csv", removal at
#'   WWTP "substance_info.csv", optional: relative error by substance can be
#'   indicated as additional column "error_conc" in "substance_info.csv")
#' @param error_removal_rate relative error in removal at WWTP
#' @param error_conc constant relative error in concentrations at WWTP outflow
#'   (default = 0.5) or "individual" if relative error by substance is included
#'   in "substance_info.csv"
#'
#' @return Function returns list with loads and standard deviations, by entry
#'   path (cso, wwtp) and by surface water catchment. Concentration in units
#'   "mg/L" and "ug/L" is automatically transformed to loads in "kg/yr". Other
#'   (unknown) units are left unchanged, resulting in "unit * m3/yr".
#'
#' @export
annual_load_sewage <- function(
  data.dir, error_removal_rate = 0.3, error_conc = 0.5
)
{
  #load data
  vol_sewage <- read_2(file = get_path_or_stop(
    data.dir, "Vol_sewage.csv", "sewage runoff"
  ))

  error_vol_sewage <- read_2(file = get_path_or_stop(
    data.dir, "Vol_sewage_relative_error.csv", "sewage runoff"
  ))

  sub_sew_info <- read_1(file = get_path_or_stop(
    data.dir, "substance_info.csv", "substance information WWTP"
  ))

  # read substance information, discard substances with lacking info
  sub_sew_info$CoutWWTP <- as.numeric(sub_sew_info$CoutWWTP)
  sub_sew_info$Retention_. <- as.numeric(sub_sew_info$Retention_.)

  selected <- !is.na(sub_sew_info$CoutWWTP)
  sub_sew_info <- sub_sew_info[selected, ]

  #set retention to zero, where information is lacking
  indices <- which(is.na(sub_sew_info$Retention_.))
  sub_sew_info$Retention_.[indices] <- 0

  #set relative error of substance concentrations (Cout)
  sub_sew_info$error_conc <- if (error_conc == "individual") {
    as.numeric(sub_sew_info$error_conc)
  } else {
    error_conc
  }

  # order of catchments
  sum_EZG <- stats::aggregate(
    vol_sewage$GESAMT, by = list(vol_sewage$SUW), FUN = "sum"
  )

  indices <- order(sum_EZG$x, decreasing = TRUE)
  sum_EZG <- sum_EZG[indices,]
  EZG_names <- sum_EZG[,1]

  # structure of calculation file
  load_sew <- data.frame(VariableName = sub_sew_info$VariableName,
                         CSO = 0, WWTP = 0, TOT = 0, stringsAsFactors=FALSE)
  error_load_sew <- load_sew

  # output list
  x_out_by_pathway_kg_yr <- list()
  error_by_pathway_kg_yr <- list()

  for (EZG in EZG_names) {

    indices <- which(vol_sewage$SUW == EZG)

    vol_sewage_EZG <- vol_sewage[indices, ]
    error_vol_sewage_EZG <- error_vol_sewage[indices, ]

    # loads of sewage based substances per pathway
    load_sew$CSO <- vol_sewage_EZG$GESAMT[1] * sub_sew_info$CoutWWTP /
      (1 - sub_sew_info$Retention_. / 100)

    load_sew$WWTP <- vol_sewage_EZG$GESAMT[2] * sub_sew_info$CoutWWTP

    load_sew$TOT <- load_sew$CSO + load_sew$WWTP

    # absolute errors in loads of sewage based substances per pathway
    error_load_sew$CSO <- load_sew$CSO * sqrt(
      error_vol_sewage_EZG$GESAMT[1]^2 +
      sub_sew_info$error_conc^2 +
      error_removal_rate^2
    )

    error_load_sew$WWTP <- load_sew$WWTP * sqrt(
      error_vol_sewage_EZG$GESAMT[2]^2 + sub_sew_info$error_conc^2
    )

    error_load_sew$TOT <- sqrt(error_load_sew$CSO^2 + error_load_sew$WWTP^2)

    # load unit for all substances in kg/yr (where concentration unit is known)
    indices_mgL <- which(sub_sew_info$UnitsAbbreviation == "mg/L")
    indices_ugL <- which(sub_sew_info$UnitsAbbreviation == "ug/L")

    load_sew[indices_mgL, -1] <- load_sew[indices_mgL, -1] / 1000
    load_sew[indices_ugL, -1] <- load_sew[indices_ugL, -1] / 1e6

    error_load_sew[indices_mgL, -1] <- error_load_sew[indices_mgL, -1] / 1000
    error_load_sew[indices_ugL, -1] <- error_load_sew[indices_ugL, -1] / 1e6

    # organize in list
    x_out_by_pathway_kg_yr[[EZG]] <- load_sew
    error_by_pathway_kg_yr[[EZG]] <- error_load_sew
  }

  # output
  list(
    by_path = x_out_by_pathway_kg_yr,
    error_by_path = error_by_pathway_kg_yr
  )
}

# getNewDetectionLimits --------------------------------------------------------
#' get detection limits for variables
#'
#' that changed (lowered) during the monitoring
#'
#' @return data frame with columns \code{VariableCode}, \code{DetectionLimit}
#' @export
getNewDetectionLimits <- function()
{
  detectionLimits.vector <- c(
    Pb = 0.5,
    Cd = 0.05,
    Cr = 0.2,
    Ni = 0.5,
    V  = 0.1
  )

  data.frame(
    VariableCode = names(detectionLimits.vector),
    DetectionLimit = as.numeric(detectionLimits.vector),
    stringsAsFactors = FALSE
  )
}

# quant25 ----------------------------------------------------------------------
#' function gives the 25 percent quantile
#'
#' @param x vector of numeric values of which to calculate the quantile
#'
#' @export
quant25 <- function(x)
{
  stats::quantile(x, c(0.25))
}

# quant75 ----------------------------------------------------------------------
#' function gives the 75 percent quantile
#'
#' @param x vector of numeric values of which to calculate the quantile
#'
#' @export
quant75<-function(x)
{
  stats::quantile(x, c(0.75))
}

# quant95 ----------------------------------------------------------------------
#' function gives the 95 percent quantile
#'
#' @param x vector of numeric values of which to calculate the quantile
#'
#' @export
quant95<-function(x)
{
  stats::quantile(x, c(0.95))
}

# geom_mean --------------------------------------------------------------------
#' function gives geometric mean
#'
#' @param x vector of numeric values of which to calculate the geometric mean
#'
#' @export
geom_mean <- function (x)
{
  exp(mean(log(x))) # same result as prod(x) ^ (1/length(x))
}

# myFunction_exp ---------------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_exp <- function(x, A, B)
{
  exp(B * x + A)
}

# myErrorFunction_exp ----------------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_exp <- function(x, A, B, RMSE)
{
  exp(B * x + A) * B * RMSE
}

# myFunction_linear ------------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_linear <- function(x, A, B)
{
  B * x + A
}

# myErrorFunction_linear -------------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_linear <- function(B, RMSE)
{
  B * RMSE
}

# myFunction_pot ---------------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_pot <- function(x, A, B)
{
  A * x^B
}

# myErrorFunction_pot ----------------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_pot <- function(x, A, B, RMSE)
{
  A * B * x^(B - 1) * RMSE
}

# myFunction_rcp ---------------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_rcp <- function(x, A, B)
{
  A / x + B
}

# myErrorFunction_rcp ----------------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_rcp <- function(x, A, RMSE)
{
  - A / (x^2) * RMSE
}

# myFunction_log ---------------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_log <- function(x, A, B)
{
  B * ln(x) + A
}

# myErrorFunction_log ----------------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_log <- function(x, B, RMSE)
{
  B / x * RMSE
}

# myFunction_polynom -----------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_polynom <- function(x, A, B)
{
  B * x^2 + A
}

# myErrorFunction_polynom ------------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_polynom <- function(x, B, RMSE)
{
  2 * B * x * RMSE
}

# myFunction_seasonal ----------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_seasonal <- function(x, A, B)
{
  c(A, B)[x]
}

# myErrorFunction_seasonal -----------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_seasonal <- function(x, RMSE_A, RMSE_B)
{
  c(RMSE_A, RMSE_B)[x]
}

# myFunction_quarterly ---------------------------------------------------------
#' @keywords internal
#' @noRd
myFunction_quarterly <- function(x, A, B, C, D)
{
  c(A, B, C, D)[x]
}

# myErrorFunction_quarterly ----------------------------------------------------
#' @keywords internal
#' @noRd
myErrorFunction_quarterly <- function(x, RMSE_A, RMSE_B, RMSE_C, RMSE_D)
{
  c(RMSE_A, RMSE_B, RMSE_C, RMSE_D)[x]
}

# ln ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
ln <- function(x)
{
  log(x, base = exp(1))
}

# annual_mean_conc_method1 -----------------------------------------------------
#' @keywords internal
#' @noRd
annual_mean_conc_method1 <- function(x_in, x_out_mean, x_out_stdev, site_names)
{
  for (site_name in site_names) {

    selected <- (x_in$SiteCode == site_name)

    values <- x_in$DataValue[selected]
    BY <- list(x_in$VariableID[selected])

    x_out_mean[[site_name]] <- stats::aggregate(
      values, by = BY, FUN = "mean"
    )[, 2]

    x_out_stdev[[site_name]] <- stats::aggregate(
      values, by = BY, FUN = "sd"
    )[, 2]
  }

  list(
    mean = x_out_mean,
    stdev = x_out_stdev
  )
}

# annual_mean_conc_method2 -----------------------------------------------------
#' @importFrom kwb.utils renameColumns
#' @keywords internal
#' @noRd
annual_mean_conc_method2 <- function(
  x_out_mean, x_out_stdev, site_names, data.dir
)
{
  #read correlation info and rain series
  x_correlations <- read_2(comment.char = "", file = get_path_or_stop(
    data.dir, "correlations_substances.csv", "correlations"
  ))

  x_rain_events <- read_2(file = get_path_or_stop(
    data.dir, "RainEvents_1961_1990_Dahlem.csv", "rain runoff"
  ))

  #keep only lines active for model calculations
  x_correlations <- x_correlations[which(x_correlations$active == "x"), ]

  #add VariableName
  variables <- kwb.ogre::OGRE_VARIABLES()
  indices <- match(x_correlations$VariableCode, variables$VariableCode)
  x_correlations$VariableName <- variables$VariableName[indices]

  #change column names in rain vents to match names of independents
  x_rain_events <- kwb.utils::renameColumns(
    dframe = x_rain_events, renames = list(
      "mean" = "RImean",
      "max" = "RImax",
      "sum" = "RD",
      "dur" = "DUR",
      "pBefore" = "DWP"
  ))

  #remove rain events with value lower than threshold in mm
  indices <- which(x_rain_events$RD > 0.8)
  x_rain_events <- x_rain_events[indices,]

  #change unit from s to h
  x_rain_events$DUR <- x_rain_events$DUR / 3600
  x_rain_events$DWP <- x_rain_events$DWP / 3600

  #add season and quarters to rain events
  x_rain_events$tBeg <- as.POSIXct(x_rain_events$tBeg)
  x_rain_events$tEnd <- as.POSIXct(x_rain_events$tEnd)

  x_rain_events$quarter <- quarters(x_rain_events$tBeg)

  x_rain_events$quarter[which(x_rain_events$quarter == "Q1")] <- 1
  x_rain_events$quarter[which(x_rain_events$quarter == "Q2")] <- 2
  x_rain_events$quarter[which(x_rain_events$quarter == "Q3")] <- 3
  x_rain_events$quarter[which(x_rain_events$quarter == "Q4")] <- 4

  x_rain_events$quarter <- as.numeric(x_rain_events$quarter)

  indices <- which(x_rain_events$quarter == 2 | x_rain_events$quarter == 3)
  x_rain_events$season <- x_rain_events$quarter
  #summer = 1
  x_rain_events$season[indices] <- 1
  #winter = 2
  x_rain_events$season[- indices] <- 2

  #define functions
  functionName <- list(
    exp = "myFunction_exp",
    linear = "myFunction_linear",
    log = "myFunction_log",
    polynom = "myFunction_polynom",
    pot = "myFunction_pot",
    quarterly = "myFunction_quarterly",
    rcp = "myFunction_rcp",
    seasonal = "myFunction_seasonal"
  )

  #define error functions
  ErrorfunctionName <- list(
    exp = "myErrorFunction_exp",
    linear = "myErrorFunction_linear",
    log = "myErrorFunction_log",
    polynom = "myErrorFunction_polynom",
    pot = "myErrorFunction_pot",
    quarterly = "myErrorFunction_quarterly",
    rcp = "myErrorFunction_rcp",
    seasonal = "myErrorFunction_seasonal"
  )

  # calculate mean and error per site and variable

  for (site_name in site_names) {

    #select correlations for site_name
    indices <- which(x_correlations$SiteCode == site_name)
    x_by_site <- x_correlations[indices,]

    #one correlation at a time for site_name
    for (i in seq_along(x_by_site$VariableName)) {

      #parametrisation for correlation function
      functionCode <- x_by_site$FunctionCode[i]
      independent <- x_rain_events[[x_by_site$x[i]]]
      A <- x_by_site$A[i]
      B <- x_by_site$B[i]
      C <- x_by_site$C[i]
      D <- x_by_site$D[i]

      args <- list(x = independent, A = A, B = B)

      if (! is.na(C)) {
        args <- c(args, C = C, D = D)
      }

      #application correlation i
      x_rain_events$conc <- do.call(
        what = functionName[[functionCode]],
        args = args
      )

      #calculate annual mean
      x_rain_events$C_times_RD <- x_rain_events$conc*x_rain_events$RD

      #select affected VariableName in output file
      indices_out <- match(x_by_site$VariableName[i], x_out_mean$VariableName)

      x_out_mean[[site_name]][indices_out] <-
        sum(x_rain_events$C_times_RD, na.rm = TRUE) /
        sum(x_rain_events$RD, na.rm = TRUE)

      #parametrisation for error function
      RMSE <- x_by_site$RMSE[i]
      RMSE_A <- x_by_site$RMSE_A[i]
      RMSE_B <- x_by_site$RMSE_B[i]
      RMSE_C <- x_by_site$RMSE_C[i]
      RMSE_D <- x_by_site$RMSE_D[i]

      if (functionCode == "exp" | functionCode == "pot") {
        args <- list(x = independent, A = A, B = B, RMSE = RMSE)
      }

      if (functionCode == "log" | functionCode == "polynom") {
        args <- list(x = independent, B = B, RMSE = RMSE)
      }

      if (functionCode == "rcp") {
        args <- list(x = independent, A = A, RMSE = RMSE)
      }

      if (functionCode == "linear") {
        args <- list(B = B, RMSE = RMSE)
      }

      if (functionCode == "seasonal") {
        args <- list(x = independent, RMSE_A = RMSE_A, RMSE_B = RMSE_B)
      }

      if (functionCode == "quarterly") {
        args <- list(
          x = independent, RMSE_A = RMSE_A, RMSE_B = RMSE_B, RMSE_C = RMSE_C,
          RMSE_D = RMSE_D
        )
      }

      #application error function i
      x_rain_events$deltaC <- do.call(
        what = ErrorfunctionName[[functionCode]],
        args = args
      )

      #calculate standard deviation of annual mean
      x_rain_events$dC_times_RD <- x_rain_events$deltaC*x_rain_events$RD

      x_out_stdev[[site_name]][indices_out] <-
        sqrt(sum(x_rain_events$dC_times_RD^2, na.rm = TRUE)) /
        sum(x_rain_events$RD, na.rm = TRUE)
    }
  }

  list(
    mean = x_out_mean,
    stdev = x_out_stdev
  )
}
