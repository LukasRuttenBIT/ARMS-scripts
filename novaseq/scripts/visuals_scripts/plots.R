library(tidyverse)
library(lubridate)
library(forcats)

# Read metadata
meta <- read.csv(
  "C:/Users/lrutt/OneDrive/Bureaublad/plots_figures/COI_batch3.4.5.6.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract observatory/unit and dates from MaterialSampleID
timeline <- meta %>%
  filter(!is.na(MaterialSampleID)) %>%
  mutate(
    match = str_match(
      MaterialSampleID,
      "^ARMS_(.+)_([0-9]{8})_([0-9]{8})_"
    ),
    observatory_unit = match[, 2],
    deployed = ymd(match[, 3]),
    retrieved = ymd(match[, 4])
  ) %>%
  filter(!is.na(observatory_unit), !is.na(deployed), !is.na(retrieved)) %>%
  
  # Remove sequencing/fraction/sample replicate duplicates
  # Put distinct here
  distinct(observatory_unit, deployed, retrieved, .keep_all = TRUE) %>%
  
  # Optional: keep only deployments/retrievals relevant to 2022-2025
  filter(retrieved >= ymd("2022-01-01"),
         deployed <= ymd("2025-12-31")) %>%
  
  # Create observatory name and unit name
  mutate(
    observatory = str_remove(observatory_unit, "_[^_]+$"),
    unit_id = str_extract(observatory_unit, "[^_]+$")
  )

timeline <- timeline %>%
  mutate(
    region = case_when(
      # Baltic
      str_detect(observatory_unit, "^TZS") ~ "Baltic",
      
      # North Sea
      str_detect(observatory_unit, "^(BelgiumCoast|SWC|Limfjord|Koster)") ~ "North Sea",
      
      # Celtic Seas
      str_detect(observatory_unit, "^(Roscoff|Plymouth|Galway)") ~ "Celtic Seas",
      
      # South European Atlantic Shelf / Bay of Biscay / Iberian Atlantic
      str_detect(observatory_unit, "^(Vigo|Getxo)") ~ "South European Atlantic Shelf",
      
      # Adriatic Sea
      str_detect(observatory_unit, "^(Rav|Rovinj|GulfOfTrieste|Piran|AdriaticEastCoast)") ~ "Adriatic Sea",
      
      # Aegean Sea
      str_detect(observatory_unit, "^Crete") ~ "Aegean Sea",
      
      # Red Sea
      str_detect(observatory_unit, "^Eilat") ~ "Red Sea",
      
      # Arctic / sub-Arctic
      str_detect(observatory_unit, "^(Svalbard|Daneborg|Nuuk)") ~ "Arctic / Sub-Arctic",
      
      # Norwegian Sea
      str_detect(observatory_unit, "^Bodo") ~ "Norwegian Sea",
      
      # Kattegat / Danish straits area
      str_detect(observatory_unit, "^Laeso") ~ "Kattegat",
      
      # Antarctic / Southern Ocean
      str_detect(observatory_unit, "^Rothera") ~ "Southern Ocean / Antarctic",
      
      TRUE ~ "Other"
    )
  )

timeline <- timeline %>%
  mutate(
    region = factor(
      region,
      levels = c(
        "Baltic",
        "North Sea",
        "Celtic Seas",
        "South European Atlantic Shelf",
        "Kattegat",
        "Norwegian Sea",
        "Arctic / Sub-Arctic",
        "Southern Ocean / Antarctic",
        "Adriatic Sea",
        "Aegean Sea",
        "Red Sea",
        "Other"
      )
    )
  )


# Order y-axis by region and deployment date
timeline <- timeline %>%
  arrange(region, deployed, observatory_unit) %>%
  mutate(
    observatory_unit = factor(observatory_unit, levels = rev(unique(observatory_unit)))
  )

p <- ggplot(timeline) +
  geom_segment(
    aes(
      x = deployed,
      xend = retrieved,
      y = observatory_unit,
      yend = observatory_unit
    ),
    linewidth = 0.4,
    colour = "grey30"
  ) +
  geom_text(
    aes(x = deployed, y = observatory_unit),
    label = "◐",
    colour = "red",
    size = 2
  ) +
  geom_text(
    aes(x = retrieved, y = observatory_unit),
    label = "◑",
    colour = "blue",
    size = 2
  ) +
  facet_grid(
    region ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = NULL,
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 11),
    plot.caption = element_text(size = 8)
  )

p




ggsave(
  "C:/Users/lrutt/OneDrive/Bureaublad/plots_figures/ARMS_MBON_sampling_timeline.png",
  p,
  width = 8,
  height = 12,
  dpi = 300
)

