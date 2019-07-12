#' SkewCalc: Estimate reproductive skew
#'
#' Documentatiion here is bare-bones, but will be updated as needed.
#'
#' @name SkewCalc
#' @docType package
NULL

#' @title KipsigisMales
#' @description Reproductive success and age for a census of Kipsigis men.
#' @format A data frame with 848 rows and 3 variables:
#' \describe{
#'   \item{\code{code}}{Personal ID code.}
#'   \item{\code{age}}{Age at death or last census.}
#'   \item{\code{rs}}{Reproductive success.}
#'}
#' @source Monique Borgerhoff Mulder
"KipsigisMales"

#' @title SukumaMales
#' @description Reproductive success and age for a census of Sukuma men.
#' @format A data frame with 63 rows and 3 variables:
#' \describe{
#'   \item{\code{code}}{Personal ID code.}
#'   \item{\code{age}}{Age at death or last census.}
#'   \item{\code{rs}}{Reproductive success.}
#'}
#' @source Monique Borgerhoff Mulder
"SukumaMales"

#' @title ColombiaRS
#' @description Reproductive success, age, sex, and group for a sample of rural Colombians.
#' @format A data frame with 289 rows and 5 variables:
#' \describe{
#'   \item{\code{code}}{Personal ID code.}
#'   \item{\code{age}}{Age at death or last census.}
#'   \item{\code{rs}}{Reproductive success. Children born.}
#'   \item{\code{sex}}{Sex.}
#'   \item{\code{group}}{Afrocolombian or Embera.}
#'}
#' @source Cody T. Ross
"ColombiaRS"
