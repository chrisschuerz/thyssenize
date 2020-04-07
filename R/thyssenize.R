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
#' @importFrom dplyr %>% bind_rows select
#' @importFrom foreach foreach %dopar%
#' @importFrom lubridate now
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom purrr set_names
#' @importFrom sf st_sfc
#'
#' @export
#'
thyssenize <- function(station_data, location, shape, shape_id, dates, cores = NULL) {
  var_name <- names(station_data)[3]
  names(station_data) <- c("date", "stat_id", "value")
  names(location) <- c("stat_id", "geometry")
  shape <- select(shape, !!shape_id) %>% set_names(c("shape_id", "geometry"))

  loc_bbox <- st_sfc(bbox_polygon(location))

  n_t <- length(dates)
  cores <- min(cores, detectCores(), n_t)
  cl <- makeCluster(cores)
  registerDoSNOW(cl)

  t0 <- now()
  progress <- function(n){
    display_progress(n, n_t, t0,  "Progress:")
  }
  opts <- list(progress = progress)

  thys_list <- foreach(t_i = 1:n_t, .options.snow = opts) %dopar% {
                         thys_i <- thyssen_i(station_data = station_data,
                                             location = location,
                                             shape = shape,
                                             bbox_sf = loc_bbox,
                                             date_i = dates[t_i])
                         return(thys_i)
                       }

  stopCluster(cl)

  thys <- bind_rows(thys_list) %>%
    set_names(c("date",shape_id, var_name))

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
#' @importFrom dplyr %>% arrange filter group_by left_join mutate select summarise
#' @importFrom purrr map_dbl set_names
#' @importFrom sf st_area st_cast st_centroid st_distance st_intersection st_sf st_union st_voronoi
#' @importFrom tibble as_tibble add_column
#'
thyssen_i <- function(station_data, location, shape, bbox_sf, date_i) {
  data_i <- station_data %>% filter(date == date_i) %>% filter(!is.na(value))
  loc_i <- filter(location, stat_id %in% data_i$stat_id)

  center_i <- suppressWarnings(st_centroid(loc_i))
  thys_i <- st_voronoi(st_union(loc_i), bbox_sf) %>%
    st_cast(.)

  stat_order_i <- st_distance(thys_i, loc_i) %>%
    as_tibble(.) %>%
    map_dbl(., which.min) %>%
    set_names(., loc_i$stat_id)

  thys_i <- st_sf(stat_id = names(sort(stat_order_i)), geom = thys_i) %>%
    mutate(stat_id = as.numeric(as.character(stat_id)))
  sub_thys_i <- suppressWarnings(st_intersection(shape, thys_i))

  sub_val_i <- sub_thys_i %>%
    mutate(area = st_area(geometry) %>% as.numeric(.)) %>%
    group_by(shape_id) %>%
    mutate(area_frct = area / sum(area)) %>%
    arrange(shape_id) %>%
    as_tibble(.) %>%
    select(shape_id, stat_id, area_frct) %>%
    left_join(., data_i, by = "stat_id") %>%
    mutate(value = value*area_frct) %>%
    group_by(shape_id) %>%
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

