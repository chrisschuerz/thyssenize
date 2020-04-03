#' Calculate area weighted values for regions based on station time series data using Thyssen ploygons
#'
#' @param station_data Tibble providing the station time series data of a measured variable. The data must be organized in three columns where the first column provides the date, the second a unique station ID, and the third column the variable value. See examples for example table.
#' @param location sf object that only holds the corresponding station identifiers and the station location.
#' @param shape sf object providing polygons for which the area weighted means should be calculated.
#' @param shape_id Polygon attribute that should be used as identifier for value averaging
#' @param dates Vector of dates for which weighted averages should be calculated.
#' @param cores Parallel computing is implemented to reduce computation time.
#'
#' @return Returns a tibble holding dates, shape_id, and averaged values
#'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom dplyr %>% bind_rows
#' @importFrom foreach foreach %dopar%
#' @importFrom lubridate now
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom sf st_sfc
#'
thyssenize <- function(station_data, location, shape, shape_id, dates, cores = NULL) {

  loc_bbox <- st_sfc(bbox_polygon(location))
  cores <- min(cores,detectCores())

  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  t0 <- now()
  progress <- function(n){
    display_progress(n, n_t, t0,  "Progress:")
  }
  opts <- list(progress = progress)

  n_t <- length(dates)

  thys_list <- foreach(t_i = 1:n_t,
                       .packages = c("tibble", "purrr", "dplyr", "lubridate", "sf"),
                       .options.snow = opts) %dopar% {
                         thys_i <- thyssen_i(station_data = station_data,
                                             location = location,
                                             shape = shape,
                                             bbox_sf = loc_bbox,
                                             date_i = dates[n_t])
                         return(thys_i)
                       }

  stopCluster(cl)
  thys <- bind_rows(thys_list)
  return(thys)
}

#' Calculate area weighted values for regions based on station data for the time step i using Thyssen ploygons
#'
#' @param station_data Tibble providing the station time series data of a measured variable. The data must be organized in three columns where the first column provides the date, the second a unique station ID, and the third column the variable value. See examples for example table.
#' @param location sf object that only holds the corresponding station identifiers and the station location.
#' @param shape sf object providing polygons for which the area weighted means should be calculated.
#' @param shape_id Polygon attribute that should be used as identifier for value averaging.
#' @param bbox_sf bounding box sf object.
#' @param date_i Lubridate ymd format date defining the date for which the average values are calculated.
#'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom dplyr %>% bind_rows
#' @importFrom foreach foreach %dopar%
#' @importFrom lubridate now
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom sf st_sfc
#'
thyssen_i <- function(station_data, location, shape, shape_id, bbox_sf, date_i) {
  var_i <- station_data %>% filter(date == date_i) %>% filter(!is.na(value))
  loc_i <- filter(location, stat_nr %in% var_i$stat_nr)

  center_i <- suppressWarnings(st_centroid(loc_i))
  thys_i <- st_voronoi(st_union(loc_i), bbox_sf) %>%
    st_cast(.)

  stat_order_i <- st_distance(thys_i, loc_i) %>%
    as_tibble(.) %>%
    map_dbl(., which.min) %>%
    set_names(., loc_i$stat_nr)

  thys_i <- st_sf(stat_nr = names(sort(stat_order_i)), geom = thys_i) %>%
    mutate(stat_nr = as.numeric(as.character(stat_nr)))
  sub_thys_i <- suppressWarnings(st_intersection(shape, thys_i))

  sub_val_i <- sub_thys_i %>%
    mutate(area = st_area(geometry) %>% as.numeric(.)) %>%
    group_by(Subbasin) %>%
    mutate(area_frct = area / sum(area)) %>%
    arrange(Subbasin) %>%
    as_data_frame(.) %>%
    select(Subbasin, stat_nr, area_frct) %>%
    left_join(., var_i, by = "stat_nr") %>%
    mutate(value = value*area_frct) %>%
    group_by(Subbasin) %>%
    summarise(value = sum(value)) %>%
    add_column(date = date_i, .before = 1)

  return(sub_val_i)
}

#' bbox
#'
#' @param x sf object for which a bounding rectangle is generated.
#'
#' @importFrom sf st_bbox st_polygon
#'
#' @keywords internal
#'
bbox_polygon <- function(x) {
  bb <- sf::st_bbox(x)

  p <- matrix(
    c(bb["xmin"], bb["ymin"],
      bb["xmin"], bb["ymax"],
      bb["xmax"], bb["ymax"],
      bb["xmax"], bb["ymin"],
      bb["xmin"], bb["ymin"]),
    ncol = 2, byrow = T
  )

  sf::st_polygon(list(p))
}


#' Display the progress if iterative processes
#'
#' @param n Iteration step
#' @param nmax Number of iterations
#' @param t0 initial time step
#' @param word Character string to define printed word
#'
#' @importFrom dplyr %>%
#' @importFrom lubridate as.period interval now seconds
#' @keywords internal
#'
display_progress <- function(n, nmax, t0, word){
  t1 <- now()
  time_elaps  <- interval(t0,t1) %>%
    round(.) %>%
    as.period(.)
  time_remain <- (as.numeric(time_elaps, "seconds")*(nmax-n)/n) %>%
    round(.) %>%
    seconds(.) %>%
    as.period(., unit = "days")
  prgs <- paste0(round(n/nmax*100, digits = 0), "%")

  cat("\r", word, prgs,
      "  Time elapsed:", as.character(time_elaps),
      "  Time remaining:", as.character(time_remain),
      "   ")
}

