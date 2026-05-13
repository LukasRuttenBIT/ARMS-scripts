install.packages("webshot", repos = "https://cloud.r-project.org")
webshot::install_phantomjs()
install.packages("leaflet", repos = "https://cloud.r-project.org")
library(leaflet)
library(tidyverse)
install.packages("mapview", repos = "https://cloud.r-project.org")
library(mapview)

# you just need a table with the sampling stations as rows


setwd("/cfs/klemming/projects/supr/naiss2025-23-46/Lukas/ARMS-scripts")

ARMS_coordinates <- read.csv("/cfs/klemming/projects/supr/naiss2025-23-46/Lukas/ARMS-scripts/metadata/generated_meta/arms_coordinates.csv")

# Show a circle at each position
m <- leaflet(data = unique(ARMS_coordinates)) %>%
  addTiles() %>%
  addCircleMarkers(
    ~decimalLongitude, ~decimalLatitude, 
    radius = 5,
    color = "firebrick",
    fill = TRUE,
    fillColor = "firebrick"
  ) %>%
  addLabelOnlyMarkers(
    ~decimalLongitude, ~decimalLatitude, 
    label = ARMS_coordinates$Station,  # Station names as labels
    labelOptions = mapply(function(station) {
      if (station == "RavMarina") {
        labelOptions(
          noHide = TRUE,
          offset = c(-20, 0),  # Move RavMarina 20px to the left
          direction = 'auto',
          textOnly = FALSE,
          style = list(
            "background-color" = "white",
            "border" = "1px solid black",
            "border-radius" = "3px",
            "padding" = "2px",
            "font-size" = "12px"
          )
        )
      } else if (station == "RavHarbour") {
        labelOptions(
          noHide = TRUE,
          offset = c(20, 0),  # Move RavHarbour 20px to the right
          direction = 'auto',
          textOnly = FALSE,
          style = list(
            "background-color" = "white",
            "border" = "1px solid black",
            "border-radius" = "3px",
            "padding" = "2px",
            "font-size" = "12px"
          )
        )
      } else {
        labelOptions(
          noHide = TRUE,
          direction = 'auto',
          textOnly = FALSE,
          style = list(
            "background-color" = "white",
            "border" = "1px solid black",
            "border-radius" = "3px",
            "padding" = "2px",
            "font-size" = "12px"
          )
        )
      }
    },ARMS_coordinates$Station, SIMPLIFY = FALSE)
  ) %>%
  fitBounds(lng1 = -20.266, lat1 = 40, lng2 = 31.0, lat2 = 59.0)

m


## save map
mapshot(m, file = "ARMS_map.png")
