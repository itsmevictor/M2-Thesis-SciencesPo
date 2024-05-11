#############################################################################################################################
# SETUP ----
#############################################################################################################################

# Packages
pacman::p_load(tidyverse, arrow, readxl, writexl, janitor, tictoc, stringi, countrycode, rdrobust, rddensity, stargazer, showtext)

# Working directory
wd <- getwd()

# Load convenient functions
lapply(list.files("/Users/victorkreitmann/Desktop/synced/life/academia-pro/useful_scripts/r_functions/", pattern = "*.R", full.names = TRUE), source)

# Load raw data
city_data <- read_parquet(paste0(wd, "/data/aggregating_allyears/city_data_full.parquet")) |> mutate(year_city = paste0(year, "_", code_commune))

font_add("CMU Sans Serif", "/Users/victorkreitmann/Library/Fonts/cmunss.ttf") # Ensure the path points to the correct font file
showtext_auto()

theme_cm_minimal <- function() {
    theme_minimal() +
        theme(text = element_text(family = "CMU Sans Serif"))
}

source(paste0(wd, "/code/analysis/aggregated_analysis/continuity_framework/analysis_city.R"))

#############################################################################################################################
# SPLIT DATASET ----
#############################################################################################################################

median_immi <- city_data |>
    filter(bothgen_UE != 1) |>
    pull(share_immi_city) |>
    median()

data <- city_data |> filter(bothgen_UE != 1)
data_high <- data |> filter(share_immi_city > median_immi)
data_low <- data |> filter(share_immi_city <= median_immi)

data2008 <- data |> filter(year == 2008)
data2014 <- data |> filter(year == 2014)
data2020 <- data |> filter(year == 2020)

bw <- ((rdbwselect(y = data$next_turnout_t1, x = data$vote_diff, c = 0, p = 2))$bws[1])

#############################################################################################################################
# TESTING ESTIMATOR
#############################################################################################################################

# F-test for the variance of the outcome variable

var.test(data_high$next_turnout_t1, data_low$next_turnout_t1)

datahigh_rd_allresults <- map_dfr(polynomials, ~ rd_estimation(df = data_high, polynomial = .x)) |> arrange(Polynomial)
datalow_rd_allresults <- map_dfr(polynomials, ~ rd_estimation(df = data_low, polynomial = .x)) |> arrange(Polynomial)

diff_in_disc_df <- bind_rows(
    datahigh_rd_allresults |> filter(Polynomial == 1 & Type == "Robust"),
    datalow_rd_allresults |> filter(Polynomial == 1 & Type == "Robust")
)

est_datahigh <- diff_in_disc_df |>
    filter(Dataset == "data_high") |>
    pull(Estimate)
est_datalow <- diff_in_disc_df |>
    filter(Dataset == "data_low") |>
    pull(Estimate)

se_datahigh <- diff_in_disc_df |>
    filter(Dataset == "data_high") |>
    pull(SE)
se_datalow <- diff_in_disc_df |>
    filter(Dataset == "data_low") |>
    pull(SE)

delta <- est_datahigh - est_datalow
se_delta <- sqrt(se_datahigh^2 + se_datalow^2)

t_stat <- delta / se_delta

#############################################################################################################################
# PLOTTING DIFFERENT RD ESTIMATES ----
#############################################################################################################################

est_bw <- function(dataset = data, bandwidth_type) {
        # Estimate model at different bandwidths
        tmp <- rdrobust(
            y = dataset$next_turnout_t1, x = dataset$vote_diff,
            bwselect = bandwidth_type, c = 0, p = 1, kernel = "tri", cluster = dataset$year_city,
            all = TRUE
        )

        # Store results
        rd_estimate <- tmp$coef[3]
        rd_confint_low <- tmp$ci[3, 1]
        rd_confint_high <- tmp$ci[3, 2]
        bw <- tmp$bws[1, 1]
        dataset <- deparse(substitute(dataset))

        # Put it in a tibble ; iterative application done later with map::map_dfr
        tibble(
            bw = bw,
            rd_estimate = rd_estimate,
            rd_confint_low = rd_confint_low,
            rd_confint_high = rd_confint_high,
            dataset = dataset
        )
    }

bandwidth_types <- c("mserd", "cerrd")

# bandwidths_high <- get_bandwidths(data_high, bws = c("mserd", "cerrd"))
# bandwidths_low <- get_bandwidths(data_low, bws = c("mserd", "cerrd"))

bandwidth_results_high <- map_dfr(bandwidth_types, ~ est_bw(dataset = data_high, bandwidth_type = .x)) |>
    mutate(bandwidth_type = c("MSE", "CER"), bw = as.numeric(bw))

bandwidth_results_low <- map_dfr(bandwidth_types, ~ est_bw(dataset = data_low, bandwidth_type = .x)) |>
    mutate(bandwidth_type = c("MSE", "CER"), bw = as.numeric(bw))

bandwidth_results_all <- bind_rows(bandwidth_results_high, bandwidth_results_low)

# print the plot for all datasets, where the legend is called "Dataset"
bandwidth_plot_all <- bandwidth_results_all |>
    ggplot(aes(x = bw, y = rd_estimate, ymin = rd_confint_low, ymax = rd_confint_high + 10, color = dataset)) +
    geom_errorbar(aes(ymin = rd_confint_low, ymax = rd_confint_high), width = 0.12) +
    geom_point(size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(
        name = "Bandwidth",
        breaks = bandwidth_results_all$bw,
        labels = paste0(round(bandwidth_results_all$bw, 3), "\n (", bandwidth_results_all$bandwidth_type, ")")
    ) +
    scale_color_discrete(name = "Dataset", labels = c("High", "Low")) +
    labs(
        y = "RD Estimate"
    ) +
    theme_cm_minimal()

bandwidth_plot_all

ggsave(paste0(wd, "/writing_thesis/figures_and_tables/results/bothsets_city_bandwidth.pdf"),
    plot = bandwidth_plot_all, scale = 1.25, width = 7.5, height = 4.5, units = "in", dpi = 450, device = "pdf"
)

#############################################################################################################################
# DIFF-IN-DISC AT ALT IMMI CUTOFFS ----
#############################################################################################################################

diff_in_disc_function <- function(cutoff) {
    cutoff_val <- city_data |>
        filter(bothgen_UE != 1) |>
        pull(share_immi_city) |>
        quantile(cutoff)

    data_high <- data |> filter(share_immi_city > cutoff_val)
    data_low <- data |> filter(share_immi_city <= cutoff_val)

    data_high_rd_results <- rd_estimation(df = data_high) |> mutate(cutoff = cutoff)
    data_low_rd_results <- rd_estimation(df = data_low) |> mutate(cutoff = cutoff)

    all_results <- bind_rows(data_high_rd_results, data_low_rd_results)

    return(all_results)
}

cutoffs <- c(0.3, 0.4, 0.6, 0.7)

diff_in_disc_results <- map_df(cutoffs, ~ diff_in_disc_function(.x)) |>
    filter(Polynomial == 1 & Type == "Robust") |>
    arrange(cutoff) |>
    mutate(across(c("Bandwidth", "Estimate", "SE", "Z", "p-value", "CI_low", "CI_high"), ~ round(as.numeric(.), 3)))
