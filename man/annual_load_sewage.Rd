% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ogreModelMain.R
\name{annual_load_sewage}
\alias{annual_load_sewage}
\title{calculates the load for each substance}
\usage{
annual_load_sewage(data.dir, error_removal_rate = 0.3, error_conc = 0.5)
}
\arguments{
\item{data.dir}{path of model data (annual mean concentrations
"substance_info.csv",  WWTP runoff volumes "Vol_sewage.csv", removal at
WWTP "substance_info.csv", optional: relative error by substance can be
indicated as additional column "error_conc" in "substance_info.csv")}

\item{error_removal_rate}{relative error in removal at WWTP}

\item{error_conc}{constant relative error in concentrations at WWTP outflow
(default = 0.5) or "individual" if relative error by substance is included
in "substance_info.csv"}
}
\value{
Function returns list with loads and standard deviations, by entry
  path (cso, wwtp) and by surface water catchment. Concentration in units
  "mg/L" and "ug/L" is automatically transformed to loads in "kg/yr". Other
  (unknown) units are left unchanged, resulting in "unit * m3/yr".
}
\description{
separates pathways (CSO and WWTP)
}
