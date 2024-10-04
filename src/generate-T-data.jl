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

percs <- c(0.15, 0.85)

rez_dat <- read.csv("data/rezende2019-supp.csv") %>%
    filter(Type == "biochem") %>%
    dplyr::select(Q10:d)

res_summary <- rez_dat %>%
    summarise(across(everything(), ~ list(quantile(., percs))))

rez_vars <- rez_dat %>%
    filter(between(Q10, res_summary$Q10[[1]][1], res_summary$Q10[[1]][2]),
        between(C, res_summary$C[[1]][1], res_summary$C[[1]][2]),
        between(Tth, res_summary$Tth[[1]][1], res_summary$Tth[[1]][2]),
        between(d, res_summary$d[[1]][1], res_summary$d[[1]][2]) )

ret <- list(names = colnames(rez_vars),
    means = colMeans(rez_vars),
    covs = true_cov(rez_vars))
"""

ret = @rget ret

save("data/T_dist.jld2",
    Dict(
        "T_dist" => MvNormal(ret[:means], ret[:covs]),
        "names" => Tuple(Symbol.(ret[:names]))
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
dir_min <- "/home/mark_holmes/computron/data/wc2.1_10m_tmin"
dir_max <- "/home/mark_holmes/computron/data/wc2.1_10m_tmax"
dir_avg <- "/home/mark_holmes/computron/data/wc2.1_10m_tavg"

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
    # avg = coef(lm(value ~ poly(latitude, 2, raw = TRUE),
    #     subset(temperature, type == "Average"))),
    max = coef(lm(value ~ poly(latitude, 2, raw = TRUE),
        subset(temperature, type == "Maximum"))))) %>% as.data.frame %>%
    rownames_to_column("type")

colnames(mods) <- c("type", "intercept", "linear_term", "quadratic_term")

ret <- list(
    min = c(subset(mods, type == "min"))[-1],
    max = c(subset(mods, type == "max"))[-1])
"""

const coefs = @rget ret

function Trange_from_latitude(x)
    cmin = coefs[:min]
    cmax = coefs[:max]
    Tmin = cmin[:intercept] + x * cmin[:linear_term] + x^2 * cmin[:quadratic_term]
    Tmax = cmax[:intercept] + x * cmax[:linear_term] + x^2 * cmax[:quadratic_term]
    return Tmin, Tmax
end

save("data/latitude_T.jld2",
    Dict(
        "min" => coefs[:min],
        "max" => coefs[:max],
        "Trange_from_latitude" => Trange_from_latitude
        )
    )
