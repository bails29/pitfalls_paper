
# =============================================================================
# DROUGHT-TO-FLOOD EVENT DETECTION: REPRODUCTION SCRIPT
# =============================================================================
#
# DESCRIPTION
#   Reproduces daily drought/flood classifications, event-index CSV files,
#   and time-series plots for six combinations of drought and flood
#   identification methods.
#
#   Event detection is performed using the external functions defined in
#   code/multi_functions_copy.R.
#
# REQUIRED INPUT
#   case_study_data/combined_flow_data.csv
#
#   The input CSV must contain:
#     ID          Catchment identifier
#     Date        Date in YYYY-MM-DD format
#     streamflow  Daily streamflow
#
# DEPENDENCIES
#   dplyr, data.table, evd, zoo, ggplot2
#
#   Install if necessary:
#   install.packages(c("dplyr", "data.table", "evd", "zoo", "ggplot2"))
#
# PROJECT STRUCTURE
#   project/
#   |-- code/
#   |   |-- multi_functions_copy.R
#   |-- case_study_data/
#   |   |-- combined_flow_data.csv
#   |-- reproduction_results/
#   |-- reproduce_events.R
#
# OUTPUTS
#   reproduction_results/
#     Daily classification CSVs (one per method combination)
#
#   reproduction_results/event_indices/
#     Event-index CSVs written by getSeasonal_p()
#
#   reproduction_results/plots/
#     Catchment time-series plots (JPEG)
#
# IMPORTANT
#   The external event-detection functions are unchanged.
#
#   Relative to the original execution, this reproduction uses:
#     - 31-observation smoothing parameter
#     - +/-15-day variable-threshold window
#     - Explicit daily classification derived from merged_condition
#
#   Event-level classification and all other event-index columns are
#   produced directly by the external functions.
#
# =============================================================================


# =============================================================================
# 1. PROJECT CONFIGURATION
# =============================================================================

# Set this to the root directory of the project.
# All subsequent paths are relative to this directory.

project_dir <- "."

setwd(project_dir)

input_file <- "./case_study_data/combined_flow_data.csv"
functions_file <- "./code/multi_functions_copy.R"

output_dir <- "./reproduction_results/git_sample"
event_dir <- file.path(output_dir, "event_indices")
plot_dir <- file.path(output_dir, "plots")


# =============================================================================
# 2. LOAD PACKAGES AND FUNCTIONS
# =============================================================================

library(dplyr)
library(data.table)
library(evd)
library(zoo)
library(ggplot2)

# Load the original event-detection functions without modification.
source(functions_file)


# =============================================================================
# 3. CREATE OUTPUT DIRECTORIES
# =============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(event_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 4. DEFINE EVENT-DETECTION METHODS
# =============================================================================

# Six combinations of drought and flood identification methods.
#
# Dmethod / Fmethod:
#   v = variable threshold
#   f = fixed threshold
#   s = seasonal threshold
#
# Percentiles and snowmelt specifications are retained from the
# original analysis.

methods <- data.table(
  Dmethod = c("v", "v", "f", "f", "s", "s"),
  Fmethod = c("f", "s", "f", "s", "f", "s"),
  approach = rep("raw", 6),
  snowmelt_method = c(
    NA_character_, "Bb", NA_character_, "Bb", "Bb", "Bb"
  ),
  drought_percentile = c(
    0.151, 0.151, 0.128, 0.128, 0.123, 0.123
  ),
  flood_percentile = rep(0.99, 6)
)


# =============================================================================
# 5. DEFINE PLOTTING PERIODS
# =============================================================================

# Catchment-specific plotting periods.
# Catchments not included here are plotted over their full available record.

plot_years <- list(
  "2409" = 2020:2024,
  "8151500" = 2020:2024,
  "11118500" = 2017:2024,
  "7112001" = 1997:2001,
  "27080" = 2008:2013,
  "108002A" = 2016:2020,
  "S. Carlo" = 2021:2023,
  "103.1.0" = 2016:2019
)


# =============================================================================
# 6. LOAD AND PREPARE STREAMFLOW DATA
# =============================================================================

flow <- fread(input_file)

# Check required input columns.
stopifnot(all(c("ID", "Date", "streamflow") %in% names(flow)))

# Standardize identifiers and dates.
flow[, ID := as.character(ID)]
flow[, date := as.Date(Date, format = "%Y-%m-%d")]

# Remove observations without valid streamflow or dates.
flow <- flow[!is.na(streamflow) & !is.na(date)]

# Sort chronologically within each catchment.
setorder(flow, ID, date)

# Check for duplicate observations.
if (flow[, anyDuplicated(paste(ID, date))] > 0L) {
  stop("Duplicate ID/date rows: resolve before running event detection.")
}

# Remove original date column and extract calendar information.
flow[, Date := NULL]

flow[, `:=`(
  year = format(date, "%Y"),
  month = format(date, "%m")
)]

# Calculate monthly mean and maximum streamflow.
flow[, `:=`(
  monthly_mean_flow = mean(streamflow),
  monthly_max_flow = max(streamflow)
), by = .(ID, year, month)]


# =============================================================================
# 7. RUN EVENT DETECTION
# =============================================================================

# Each method combination is applied to all catchments.
#
# getSeasonal_p() returns daily results and writes event-index CSVs
# directly to event_dir.
#
# The original event-index files and their event-level class columns
# are not modified by this script.

for (k in seq_len(nrow(methods))) {
  
  m <- methods[k]
  
  message("Running ", m$Dmethod, "_", m$Fmethod)
  
  pieces <- list()
  
  # ---------------------------------------------------------------------------
  # 7.1 Detect events for each catchment
  # ---------------------------------------------------------------------------
  
  for (site in unique(flow$ID)) {
    
    message("  ", site)
    
    site_flow <- copy(flow[ID == site])
    setorder(site_flow, date)
    
    result <- getSeasonal_p(
      df = site_flow,
      Dmethod = m$Dmethod,
      Fmethod = m$Fmethod,
      approach = m$approach,
      drought_percentile = m$drought_percentile,
      flood_percentile = m$flood_percentile,
      smoothing_param = 31,
      min_length = 30,
      n = 15,
      min_recovery = 15,
      processed_events_path = paste0(event_dir, "/"),
      snowmelt_method = m$snowmelt_method
    )
    
    pieces[[site]] <- as.data.table(result)
  }
  
  # ---------------------------------------------------------------------------
  # 7.2 Combine daily results
  # ---------------------------------------------------------------------------
  
  daily <- rbindlist(
    pieces,
    use.names = TRUE,
    fill = TRUE
  )
  
  setorder(daily, ID, date)
  
  
  # ---------------------------------------------------------------------------
  # 7.3 Add daily classification
  # ---------------------------------------------------------------------------
  
  # The helper returns merged_condition but does not retain the daily
  # class column.
  #
  # Classification:
  #   drought -> drought
  #   flood   -> flood
  #   NA      -> NA
  #   other   -> normal
  #
  # This addition does not modify the event-index CSVs.
  
  if (!"merged_condition" %in% names(daily)) {
    stop(
      "Cannot create class: merged_condition is missing from the daily results."
    )
  }
  
  daily[, class := fifelse(
    merged_condition == "drought", "drought",
    fifelse(
      merged_condition == "flood", "flood",
      fifelse(
        is.na(merged_condition),
        NA_character_,
        "normal"
      )
    )
  )]
  
  stopifnot("class" %in% names(daily))
  
  
  # ---------------------------------------------------------------------------
  # 7.4 Save daily classification CSV
  # ---------------------------------------------------------------------------
  
  daily_file <- file.path(
    output_dir,
    paste0(
      m$Dmethod, "_",
      m$Fmethod,
      "_new_drought_calc.csv"
    )
  )
  
  fwrite(daily, daily_file)
  
  # Verify that the written CSV includes the daily class column.
  written_header <- names(fread(daily_file, nrows = 0L))
  
  if (!"class" %in% written_header) {
    stop("Saved CSV is missing class: ", daily_file)
  }
  
  message("  Saved daily CSV with class: ", daily_file)
  
  
  # ---------------------------------------------------------------------------
  # 7.5 Calculate additional variables for plotting
  # ---------------------------------------------------------------------------
  
  # These variables are added AFTER writing the daily CSV.
  # This preserves the original output schema and processing order.
  
  daily[, water_year := as.integer(year) +
          as.integer(month %in% c("09", "10", "11", "12"))]
  
  if (!is.na(m$snowmelt_method) && "season" %in% names(daily)) {
    
    daily[, `:=`(
      season_start = min(date),
      season_end = max(date)
    ), by = .(ID, water_year, season)]
    
  }
  
  
  # =============================================================================
  # 8. GENERATE CATCHMENT TIME-SERIES PLOTS
  # =============================================================================
  
  # Plot the original streamflow, smoothed streamflow, drought and flood
  # thresholds, and identified drought/flood conditions.
  #
  # Plot appearance, dimensions, and naming are retained.
  
  for (site in unique(daily$ID)) {
    
    plot_data <- daily[ID == site]
    
    # Apply catchment-specific plotting period where defined.
    if (site %in% names(plot_years)) {
      plot_data <- plot_data[
        as.integer(year) %in% plot_years[[site]]
      ]
    }
    
    if (!nrow(plot_data)) next
    
    # Observations classified as drought or flood.
    dots <- plot_data[
      !is.na(merged_condition) &
        merged_condition != "normal"
    ]
    
    # Create time-series plot.
    p <- ggplot(plot_data, aes(x = date)) +
      
      geom_line(
        aes(y = streamflow),
        color = "grey50",
        linewidth = 0.5
      ) +
      
      geom_line(
        aes(y = smoothed, color = "smoothed"),
        linewidth = 0.5
      ) +
      
      geom_line(
        aes(y = threshold, color = "drought threshold"),
        linewidth = 0.5,
        linetype = "dashed"
      ) +
      
      geom_line(
        aes(y = flood_threshold, color = "flood threshold"),
        linewidth = 0.5,
        linetype = "dashed"
      ) +
      
      geom_point(
        data = dots,
        aes(y = 100, color = merged_condition),
        size = 2
      ) +
      
      scale_color_manual(values = c(
        drought = "#7b3294",
        `drought threshold` = "#c2a5cf",
        flood = "#008837",
        `flood threshold` = "#a6dba0",
        smoothed = "black"
      )) +
      
      scale_x_date(
        date_breaks = "3 months",
        date_labels = "%m/%y"
      ) +
      
      labs(
        title = paste(
          site,
          "drought:", m$Dmethod,
          "flood:", m$Fmethod
        ),
        x = "Date",
        y = "Streamflow",
        color = NULL
      ) +
      
      theme_bw() +
      
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5
        ),
        legend.position = "bottom"
      )
    
    
    # -------------------------------------------------------------------------
    # 8.1 Save plot
    # -------------------------------------------------------------------------
    
    # The original naming convention uses "NA" for methods
    # without a seasonal snowmelt specification.
    
    plot_name <- paste0(
      site, "_short_",
      m$Dmethod,
      m$Fmethod,
      m$approach,
      ifelse(
        is.na(m$snowmelt_method),
        "NA",
        m$snowmelt_method
      ),
      ".jpeg"
    )
    
    ggsave(
      file.path(plot_dir, plot_name),
      p,
      dpi = 300,
      units = "mm",
      width = 210,
      height = 150
    )
    
  }
}

message("Finished. Daily CSVs contain the additional class column.")
