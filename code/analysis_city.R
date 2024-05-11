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
# ESTIMATION PER SE ----
#############################################################################################################################

# Creating year & city fixed effects
data1 <- data |> mutate(year_nextelection = as.factor(year_nextelection), code_commune = as.factor(code_commune))

year_nextelection_dummies <- model.matrix(~ year_nextelection - 1, data = data1)
city_dummies <- model.matrix(~ code_commune - 1, data = data1)

covariates <- cbind(
    year_nextelection_dummies,
    city_dummies
)

rd_estimation <- function(df = data, polynomial = 1, bandwidth = NULL, covariates = NULL) {
    if (exists("bandwidth")) {
        rd_data1 <- rdrobust(y = df$next_turnout_t1, x = df$vote_diff, c = 0, p = polynomial, kernel = "tri", h = bandwidth, cluster = df$year_city, all = TRUE)
    } else {
        rd_data1 <- rdrobust(y = df$next_turnout_t1, x = df$vote_diff, c = 0, p = polynomial, kernel = "tri", bwselect = "mserd", cluster = df$year_city, all = TRUE)
    }

    results_conv <- cbind(
        rd_data1$coef[1],
        rd_data1$se[1],
        rd_data1$z[1],
        rd_data1$pv[1],
        rd_data1$ci[1, 1],
        rd_data1$ci[1, 2]
    )

    results_biasco <- cbind(
        rd_data1$coef[2],
        rd_data1$se[2],
        rd_data1$z[1],
        rd_data1$pv[2],
        rd_data1$ci[2, 1],
        rd_data1$ci[2, 2]
    )

    results_robust <- cbind(
        rd_data1$coef[3],
        rd_data1$se[3],
        rd_data1$z[1],
        rd_data1$pv[3],
        rd_data1$ci[3, 1],
        rd_data1$ci[3, 2]
    )

    if (exists("covariates_prep")) {
        results_info <- cbind(
            deparse(substitute(df)),
            rd_data1$bws[1, 1],
            paste(names(covariates_prep), collapse = ", "),
            rd_data1$N[1],
            rd_data1$N[2],
            rd_data1$N_h[1],
            rd_data1$N_h[2],
            rd_data1$p
        ) |>
            as_tibble() |>
            set_names(c("Dataset", "Bandwidth", "Covariates", "Left observations", "Right observations", "Effective obs. left", "Effective obs. right", "Polynomial"))
    } else {
        results_info <- cbind(
            deparse(substitute(df)),
            rd_data1$bws[1, 1],
            "None",
            rd_data1$N[1],
            rd_data1$N[2],
            rd_data1$N_h[1],
            rd_data1$N_h[2],
            rd_data1$p
        ) |>
            as_tibble() |>
            set_names(c("Dataset", "Bandwidth", "Covariates", "Left observations", "Right observations", "Effective obs. left", "Effective obs. right", "Polynomial"))
    }

    all_results <- rbind(results_conv, results_biasco, results_robust) |>
        as_tibble() |>
        set_names(c("Estimate", "SE", "Z", "p-value", "CI_low", "CI_high"))


    main_results <- cbind(results_info, all_results) |>
        mutate(
            Type = c("Conventional", "Bias-corrected", "Robust"),
            .before = Dataset,
        )

    return(main_results)
}

polynomials <- c(1, 2, 3, 4)

# This estimates the model with all possible polynomials, for each of the three datasets
batch_estimate <- function(dataset, dataset_name) {
    batch_results <- map_dfr(polynomials, ~ rd_estimation(df = dataset, polynomial = .x)) |> arrange(Dataset, Polynomial)

    batch_results$Dataset <- dataset_name

    return(batch_results)
}

full_results <- bind_rows(
    batch_estimate(data_high, "data_high"),
    batch_estimate(data_low, "data_low")
) |>
    filter(Polynomial == 1 & Type == "Robust")
    # select(Dataset, Type, Estimate, SE, Z, 'p-value', CI_low, CI_high)

#############################################################################################################################
# ESTIMATION PLOTS ----
#############################################################################################################################

results_plot <- function(dataset = data, polynomial) {
    covariates <- model.matrix(~ year_nextelection + 0, data = dataset)

    tmp <- rdplot(
        y = data$next_turnout_t1,
        x = data$vote_diff,
        c = 0,
        p = polynomial,
        x.label = "Difference in vote share (X)",
        y.label = "Turnout in the next election (Y)",
        title = "",
        hide = TRUE
    )

    plot <- tmp$rdplot + theme_cm_minimal()

    return(plot)
}

polynomials <- c(1, 2, 3, 4)
# apply the function to print four plots, based on each value of the polynomials
list_of_plots <- map(polynomials, ~ results_plot(data, .x))
# save the plots separately
for (i in seq_along(list_of_plots)) {
    plot_name <- paste0(wd, "/writing_thesis/figures_and_tables/results/continuity_results_", sprintf("plot%d.pdf", i)) # Creating a filename for each plot
    ggsave(plot_name, plot = list_of_plots[[i]], width = 5, height = 5)
}


tmp <- rdplot(
    y = data$next_turnout_t1,
    x = data$vote_diff,
    c = 0,
    p = 1,
    x.label = "Difference in vote share (X)",
    y.label = "Turnout in the next election (Y)",
    title = "",
    # hide = TRUE
)$rdplot +
    theme_cm_minimal()

tmp


#############################################################################################################################
# DENSITY TEST & PLOTS ----
#############################################################################################################################

# Density test
density_test <- rddensity(
    X = data$vote_diff,
    c = 0,
)

# Density plot
density_plot <- rdplotdensity(
    rdd = rddensity(X = data$vote_diff, c = 0),
    X = data$vote_diff,
    lwd = 1.10,
    alpha = 0.05,
    CIshade = 0.23,
    # CItype = "region",
    CIcol = c("#EA3324", "#00b7ff"),
    histFillShade = 0.20,
    histFillCol = "#355b67",
    # title = "Density plot for the running variable",
    xlabel = "Difference in vote share (X)",
    ylabel = "Density",
    noPlot = TRUE
)$Estplot +
    # add a black "full" (not dotted) vertical line at 0
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    # change the color of the lines on each side of the cutoff (for the density) to orange and blue
    scale_color_manual(values = c("red", "blue")) +
    theme_cm_minimal() +
    # remove legend
    theme(legend.position = "none")

density_plot

ggsave(paste0(wd, "/writing_thesis/figures_and_tables/results/city_density_plot.pdf"),
    plot = density_plot,
    scale = 1.25, width = 7.5, height = 4.5, units = "in", dpi = 450, device = "pdf"
)

#############################################################################################################################
# PLACEBO CUTOFFS ----
#############################################################################################################################

cutoff_table <- function(cutoff) {
    rd_alt_cutoff <- rdrobust(y = data$next_turnout_t1, x = data$vote_diff, c = cutoff, p = 2, kernel = "tri", bwselect = "mserd", cluster = data$year_city, all = TRUE)

    print(paste0("Alternative cutoff: ", cutoff))
    print(paste0("Bandwidth: ", round(rd_alt_cutoff$bws[1, 1], 3)))
    print(paste0("RD Estimator: ", round(rd_alt_cutoff$coef[3], 3)))
    print(paste0("p-value: ", round(rd_alt_cutoff$pv[3], 3)))
    print(paste0("Confidence intervals: ", paste(rd_alt_cutoff$ci[3, ], collapse = " ")))

    print(paste0("Left observations: ", (rd_alt_cutoff$N_h |> str_split(" ") |> unlist())[1]))
    print(paste0("Right observations: ", (rd_alt_cutoff$N_h |> str_split(" ") |> unlist())[2]))
}

cutoff_table(-3)


#############################################################################################################################
# SENSITIVITY TO BANDWIDTH CHOICE ----
#############################################################################################################################

# Estimation
get_bandwidths <- function(dataset = data, bws = c("cerrd", "cerrd2", "mserd", "mserd2")) {
    bandwidths <- list()

    # Calculate "cerrd" bandwidth if requested
    if ("cerrd" %in% bws) {
        cer_bw <- rdbwselect(
            y = dataset$next_turnout_t1, x = dataset$vote_diff,
            kernel = "tri", cluster = dataset$year_city, 
            bwselect = "cerrd"
        )$bws[1]
        bandwidths[["cerrd"]] <- cer_bw
    }

    # Calculate "cerrd2" bandwidth if requested
    if ("cerrd2" %in% bws) {
        cer_bw <- rdbwselect(
            y = dataset$next_turnout_t1, x = dataset$vote_diff,
            kernel = "tri", cluster = dataset$year_city, 
            bwselect = "cerrd"
        )$bws[1]
        cer_bw2 <- cer_bw * 2
        bandwidths[["cerrd2"]] <- cer_bw2
    }

    # Calculate "mserd" bandwidth if requested
    if ("mserd" %in% bws) {
        mserd_bw <- rdbwselect(
            y = dataset$next_turnout_t1, x = dataset$vote_diff,
            kernel = "tri", cluster = dataset$year_city,
            bwselect = "mserd"
        )$bws[1]
        bandwidths[["mserd"]] <- mserd_bw
    }

    # Calculate "mserd2" bandwidth if requested
    if ("mserd2" %in% bws) {
        mserd_bw <- rdbwselect(
            y = dataset$next_turnout_t1, x = dataset$vote_diff,
            kernel = "tri", cluster = dataset$year_city, 
            c = 0, bwselect = "mserd"
        )$bws[1]
        mserd_bw2 <- mserd_bw * 2
        bandwidths[["mserd2"]] <- mserd_bw2
    }

    return(bandwidths)
}

bandwidth_estimate <- function(dataset = data, bandwidth) {
    # Estimate model at different bandwidths
    tmp <- rdrobust(
        y = dataset$next_turnout_t1, x = dataset$vote_diff,
        h = bandwidth, c = 0, p = 1, kernel = "tri", cluster = dataset$year_city, 
        all = TRUE
    )
    # rd_data1 <- rdrobust(y = df$next_turnout_t1, x = df$vote_diff, c = 0, p = polynomial, kernel = "tri", bwselect = "mserd", cluster = df$year_city, all = TRUE)

    # Store results
    rd_estimate <- tmp$coef[3]
    rd_confint_low <- tmp$ci[3, 1]
    rd_confint_high <- tmp$ci[3, 2]
    bw <- bandwidth
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

bandwidths <- get_bandwidths(data)
bandwidth_results <- map_dfr(bandwidths, ~ bandwidth_estimate(dataset = data, bandwidth = .x)) |>
    mutate(bandwidth_type = c("CER", "CER 2", "MSE", "MSE 2")) 

# Create a plot based on these results
bandwidth_plot <- bandwidth_results |>
    ggplot(aes(x = bw, y = rd_estimate, ymin = rd_confint_low, ymax = rd_confint_high)) +
    geom_errorbar(aes(ymin = rd_confint_low, ymax = rd_confint_high), width = 0.12, color = "#000085") +
    geom_point(color = "#EA3323", size = 4) + # Separate layer for the red points
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(
        name = "Bandwidth",
        breaks = bandwidth_results$bw,
        labels = paste0(round(bandwidth_results$bw, 3), "\n (", bandwidth_results$bandwidth_type, ")")
    ) +
    labs(
        # title = paste0("Sensitivity to bandwidth choice (", bandwidth_results$dataset[1], ")"),
        y = "RD Estimate"
    ) +
    theme_cm_minimal()

# Print the plot
bandwidth_plot

# Save the plot
ggsave(paste0(wd, "/writing_thesis/figures_and_tables/results/city_bandwidths.pdf"),
    plot = bandwidth_plot, scale = 1.25, width = 7.5, height = 4.5, units = "in", dpi = 450, device = "pdf"
)


