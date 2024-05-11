# Code to perform the local-randomization analysis.

#############################################################################################################################
# SETUP ----
#############################################################################################################################

# Packagesw
pacman::p_load(tidyverse, arrow, readxl, writexl, janitor, tictoc, stringi, countrycode, rdlocrand, rddensity, stringdist, showtext)

# Working directory
wd <- getwd()

# Load convenient functions
lapply(list.files("/Users/victorkreitmann/Desktop/synced/life/academia-pro/useful_scripts/r_functions/", pattern = "*.R", full.names = TRUE), source)

# Font stuff
font_add("CMU Sans Serif", "/Users/victorkreitmann/Library/Fonts/cmunss.ttf") # Ensure the path points to the correct font file
showtext_auto()

theme_cm_minimal <- function() {
    theme_minimal() +
        theme(text = element_text(family = "CMU Sans Serif"))
}

# Load raw data
bv_data <- read_parquet(paste0(wd, "/data/aggregating_allyears/bv_data_full.parquet")) |> mutate(year_bv = paste0(year, "_", id_bv_miom))

#############################################################################################################################
# SPLIT DATASET ----
#############################################################################################################################

data <- bv_data |>
    filter(bothgen_UE != 1) |>
    mutate(share_immi_bv_cat = as.integer(share_immi_bv_cat))

data2008 <- data |> filter(year == 2008)
data2014 <- data |> filter(year == 2014)
data2020 <- data |> filter(year == 2020)

median_immi <- median(data$share_immi_bv)

data_high <- data |> filter(share_immi_bv > median_immi)
data_low <- data |> filter(share_immi_bv <= median_immi)

#############################################################################################################################
# SIMPLE SCATTERPLOT ----
#############################################################################################################################

sampled_data <- data |>
    filter(year == 2020) |>
    distinct(code_commune, year) |>
    sample_n(80)

data |>
    inner_join(sampled_data, by = c("year", "code_commune")) |>
    ggplot(aes(x = vote_diff, y = next_turnout_t1, color = code_commune)) +
    geom_point(alpha = 0.40) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    labs(
        # title = "Scatterplot of vote share difference and turnout",
        x = "Difference in vote share (X)",
        y = "Turnout at the next election (Y)"
    ) +
    theme_cm_minimal() +
    theme(legend.position = "none") +
    # cut the plot for values below y = 55
    coord_cartesian(ylim = c(0, 55))

ggsave(
    paste0(wd, "/writing_thesis/figures_and_tables/empirical_framework/bv_scatterplot.pdf"),
    plot = last_plot(), scale = 0.90, width = 9, height = 6, units = "in", dpi = 450, device = "pdf"
)


#############################################################################################################################
# WINDOW SELECTION ----
########################### W##################################################################################################

rdlocrand::rdwinselect(
    R = data$vote_diff,
    X = data |> select(
        # turnout_t1_city
        n_voters_bv, 
        share_immi_bv),
    c = 0,
    wobs = 2,
)




#############################################################################################################################
# DENSITY TESTING ----
#############################################################################################################################

binom_table <- function(w_params) {
    tmp <- rdlocrand::rdrandinf(
        Y = data$next_turnout_t1,
        R = data$vote_diff,
        cutoff = 0,
        wl = w_params[1],
        wr = w_params[2],
        statistic = "diffmeans",
        seed = "999",
        quietly = TRUE,
        wmasspoints = TRUE
    )

    control_obs <- (tmp$sumstats[2, ] |> as.character() |> strsplit(" ") |> unlist())[1] |> as.numeric()
    treated_obs <- (tmp$sumstats[2, ] |> as.character() |> strsplit(" ") |> unlist())[2] |> as.numeric()
    total_obs <- as.numeric(control_obs) + as.numeric(treated_obs)

    test <- binom.test(control_obs, total_obs, 1 / 2)
    print(test$p.value)

    w_left <- w_params[1]
    w_right <- w_params[2]
    n_left <- control_obs
    n_right <- treated_obs
    p_value <- test$p.value

    tbl <- tibble(
        dataset = deparse(substitute(data)),
        w_left = w_left,
        w_right = w_right,
        n_left = n_left,
        n_right = n_right,
        p_value = p_value
    )

    return(tbl)
}

tmp <- rdlocrand::rdrandinf(
    Y = data$next_turnout_t1,
    R = data$vote_diff,
    cutoff = 0,
    wl = -0.75,
    wr = 0.75,
    statistic = "diffmeans",
    seed = "999"
)

tmp <- rdlocrand::rdrandinf(
    Y = data$next_turnout_t1,
    R = data$vote_diff,
    cutoff = 0,
    wl = -1,
    wr = 1,
    statistic = "diffmeans",
    seed = "999",
    # quietly = TRUE,
    wmasspoints = TRUE
)



#############################################################################################################################
# ESTIMATION PER SE ----
#############################################################################################################################

rd_inference <- function(df = data, w_params) {
    ci_vec <- c(0.05, seq(from = -20, to = 20, by = 0.1))
    bernoulli_probs <- rep(0.5, nrow(df))

    tmp <- rdlocrand::rdrandinf(
        Y = df$next_turnout_t1,
        R = df$vote_diff,
        wl = w_params[1],
        wr = w_params[2],
        seed = 999,
        bernoulli = bernoulli_probs,
        wmasspoints = TRUE,
        ci = ci_vec,
        d = -10,
        quietly = TRUE
    )

    t_test <- tmp$obs.stat[1]
    p_value <- tmp$p.value
    confint_left <- tmp$ci[1]
    confint_right <- tmp$ci[2]
    wleft <- w_params[1]
    wright <- w_params[2]
    nleft <- tmp$sumstats[2, 1]
    nright <- tmp$sumstats[2, 2]
    binom_test <- binom_table(w_params) |> select(p_value) |> pull()

    tbl <- tibble(
        dataset = deparse(substitute(df)),
        w_left = wleft,
        w_right = wright,
        n_left = nleft,
        n_right = nright,
        t_test = t_test,
        p_value = p_value,
        confint_left = confint_left,
        confint_right = confint_right,
        binom_test = binom_test
    )

    return(tbl)
}

params_list <- list(
    c(-2, 2),
    c(-1.75, 1.75),
    c(-1.5, 1.5),
    c(-1.25, 1.25),
    c(-1, 1),
    c(-0.75, 0.75),
    c(-0.5, 0.5),
    c(-0.25, 0.25)
)

rd_results <- lapply(params_list, rd_inference, df = data) |> data.table::rbindlist()

binom_results <- lapply(params_list, binom_table) |> data.table::rbindlist()




df = data
ci_vec <- c(0.05, seq(from = -20, to = 20, by = 0.1))
w_params = c(-2, 2)
# create a vector equal to the length of the dataframe that specifices the probabilities of treatment for each unit (0.5, 0.5)
bernoulli_probs <- rep(0.5, nrow(df))

tmp1 <- rdlocrand::rdrandinf(
    Y = df$next_turnout_t1,
    R = df$vote_diff,
    wl = w_params[1],
    wr = w_params[2],
    bernoulli = bernoulli_probs,
    seed = 999,
    wmasspoints = TRUE,
    ci = ci_vec,
    d = -10,
    # quietly = TRUE
)



#############################################################################################################################
# DENSITY PLOT ----
#############################################################################################################################

density <- rddensity(X = data$vote_diff, c = 0, bino = TRUE, binoW = 0.25, binoNStep = 5, massPoints = TRUE)

density_plot <- rdplotdensity(
    rdd = density,
    X = data$vote_diff,
    plotRange = c(-10, 10),
    lwd = 1.10,
    alpha = 0.05,
    CIshade = 0.23,
    histFillShade = 0.20,
    CIcol = c("#EA3324", "#00b7ff"),
    histFillCol = "#355b67",
    # title = "Density plot for the running variable",
    xlabel = "Difference in vote share (X)",
    ylabel = "Density",
    noPlot = TRUE
)$Estplot +
    # add a black "full" (not dotted) vertical line at 0
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    # change the color of the lines on each side of the cutoff (for the density) to orange and blue
    scale_color_manual(values = c("red", "blue"))

density_plot

ggsave(
    paste0(wd, "/writing_thesis/figures_and_tables/results/bv_density_plot.pdf"),
    plot = density_plot,
    scale = 1.25, width = 7.5, height = 4.5, units = "in", dpi = 450, device = "pdf"
)
