using RCall, JLD2, Distributions

# ==============================================================================
# TPC
# ==============================================================================

R"""
library("tidyverse")

true_cov <- function(mat) {
    n <- ncol(mat)
    M <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:n) {
            M[i,j] <- mean((mat[,i] - mean(mat[,i])) * (mat[,j] - mean(mat[,j])))
        }
    }
    rownames(M) <- colnames(M) <- colnames(mat)
    return(M)
}

percs <- c(0.01, 0.99)

dat <- read.csv("data/rezende2019-supp.csv") %>%
    dplyr::filter(Type != "run") %>%
    dplyr::select(Type, Q10:d) %>%
    # these two appear to weird so we log now to exp later
    mutate(C = log(C), d = log(d)) %>%
    rename("y" = C, "q" = Q10, "z" = Tth, "w" = d)

# ------------------------------------------------------------------------------
# plants
# ------------------------------------------------------------------------------

plant_dat <- dat %>%
    dplyr::filter(Type != "run") %>%
    dplyr::select(-Type)

plant_summary <- plant_dat %>%
    summarise(across(everything(), ~ list(quantile(., percs))))

plant_vars <- plant_dat %>%
    filter(
        between(q, plant_summary$q[[1]][1], plant_summary$q[[1]][2]),
        between(y, plant_summary$y[[1]][1], plant_summary$y[[1]][2]),
        between(z, plant_summary$z[[1]][1], plant_summary$z[[1]][2]),
        between(w, plant_summary$w[[1]][1], plant_summary$w[[1]][2])
        )

# ------------------------------------------------------------------------------
# insects
# ------------------------------------------------------------------------------

insect_dat <- dat %>%
    dplyr::filter(Type != "fit") %>%
    dplyr::select(-Type)

insect_summary <- insect_dat %>%
    summarise(across(everything(), ~ list(quantile(., percs))))

insect_vars <- insect_dat %>%
    filter(
        between(q, insect_summary$q[[1]][1], insect_summary$q[[1]][2]),
        between(y, insect_summary$y[[1]][1], insect_summary$y[[1]][2]),
        between(z, insect_summary$z[[1]][1], insect_summary$z[[1]][2]),
        between(w, insect_summary$w[[1]][1], insect_summary$w[[1]][2])
        )

# ------------------------------------------------------------------------------
# export
# ------------------------------------------------------------------------------

ret <- list(
    #
    plant_names = colnames(plant_vars),
    plant_means = colMeans(plant_vars),
    plant_covs = true_cov(plant_vars),
    #
    insect_names = colnames(insect_vars),
    insect_means = colMeans(insect_vars),
    insect_covs = true_cov(insect_vars)
    )
"""

ret = @rget ret

save("data/T_dist.jld2",
    Dict(
        "plant_T_dist" => MvNormal(ret[:plant_means], ret[:plant_covs]),
        "plant_names" => Tuple(Symbol.(ret[:plant_names])),
        "insect_T_dist" => MvNormal(ret[:insect_means], ret[:insect_covs]),
        "insect_names" => Tuple(Symbol.(ret[:insect_names]))
        )
    )

# ==============================================================================
# temperature distribution
# ==============================================================================

# thank you to https://askubuntu.com/questions/1211782/curl-openssl-4-not-found-required-by-curl for resolving the CURL issue that was arising here
R"""
library("tidyverse")
library("raster")
library("terra")

compute_average <- function(stack_summary, name) {
    lats <- seq(ymin(stack_summary), ymax(stack_summary), length = nrow(stack_summary))
    sum_by_lat <- rowSums(stack_summary, na.rm = TRUE)
    n_obs_by_lat <- rowSums(!is.na(stack_summary), na.rm = TRUE)
    average_by_lat <- sum_by_lat / n_obs_by_lat
    return (data.frame(type = name, latitude = lats, value = average_by_lat))
}

# data from https://worldclim.org/data/worldclim21.html
dir_min <- "data/wc2.1_10m_tmin"
dir_avg <- "data/wc2.1_10m_tavg"
dir_max <- "data/wc2.1_10m_tmax"

# we need the yearly average as well as the average minimum in the winter and the average maximum in the summer
min_files <- list.files(dir_min, full.names = TRUE)
avg_files <- list.files(dir_avg, full.names = TRUE)
max_files <- list.files(dir_max, full.names = TRUE)

min_stack <- stack(lapply(min_files, function(i) {raster(i)}))
avg_stack <- stack(lapply(avg_files, function(i) {raster(i)}))
max_stack <- stack(lapply(max_files, function(i) {raster(i)}))

min_stack_summary <- min(min_stack, na.rm = TRUE)
avg_stack_summary <- mean(avg_stack, na.rm = TRUE)
max_stack_summary <- max(max_stack, na.rm = TRUE)

Tmin_lats <- compute_average(min_stack_summary, "tmin")
Tavg_lats <- compute_average(avg_stack_summary, "tavg")
Tmax_lats <- compute_average(max_stack_summary, "tmax")

temperature <- bind_rows(Tmin_lats, Tmax_lats, Tavg_lats) %>%
    mutate(type = factor(type,
        levels = c("tmin", "tavg", "tmax"),
        labels = c("Minimum", "Average", "Maximum")))

t_colours <- setNames(hcl.colors(3, "Berlin"), levels(temperature$type))

mods <- t(data.frame(
    min = coef(lm(value ~ poly(latitude, 2, raw = TRUE),
        subset(temperature, type == "Minimum"))),
    max = coef(lm(value ~ poly(latitude, 2, raw = TRUE),
        subset(temperature, type == "Maximum"))))) %>% as.data.frame %>%
    rownames_to_column("type")

colnames(mods) <- c("type", "intercept", "linear_term", "quadratic_term")

ret <- list(
    min = c(subset(mods, type == "min"))[-1],
    max = c(subset(mods, type == "max"))[-1])
"""

const coefs = @rget ret

save("data/latitude_T.jld2",
    Dict(
        "min" => coefs[:min],
        "max" => coefs[:max]
        )
    )
