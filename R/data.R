#' Stocks volatility
#'
#' A dataframe of 2538 stocks' volatility, as well as their corresponding company name. The volatility is estimated by daily high and low price and averaged throughout 2015-2019 (See paper "Controlling the False Split Rate in Tree-Based Aggregation" for more details).
#'
#' @source Center for Research in Security Prices (CRSP)
"stocks_volatility"


#' Hierarchical tree based on NAICS for companies in the stocks dataset
#'
#' A list of length-\code{702} that stores the tree structure (referred to as \code{hc_list} format in the package) for the 2538 stocks. The ith item in the list contains the child nodes of the ith node in the tree. The negative values in the list indicate leaf nodes and the positive values correspond to interior nodes.
#' The tree structure is formed based on North American Industry Classification System (NAICS).
#' @source Information of North American Industry Classification System available at \url{https://www.census.gov/naics/?58967?yearbck=2017}. Company's NAICS code from CRSP/Compustat merged database.
"stocks_tree"


#' Hierarchical tree based on neighborhoods in NYC
#'
#' A list of length-\code{190} that stores the tree structure (referred to as \code{hc_list} format in the package) for the 194 neighborhoods in New York city. The ith item in the list contains the child nodes of the ith node in the tree. The negative values in the list indicate leaf nodes and the positive values correspond to interior nodes.
#' The tree structure is formed based on Neighborhood Tabulation Areas (NTAs).
"taxi_tree"
