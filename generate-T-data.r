library("raster")
library("tidyverse")
library("geodata")

# ==============================================================================
# TPC
# ==============================================================================

true_cov <- function (mat) {
    n <- ncol(mat)
    M <- matrix(0, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            M[i,j] <- mean((mat[,i] - mean(mat[,i])) * (mat[,j] - mean(mat[,j])))
        }
    }
    return(M)
}

dat <- read.csv("data/rezende2019-supp.csv") %>%
    filter(Type != "run") %>%
    select(Type, Q10, C, Tth, d) %>%
    mutate(C = log(C), d = log(d)) %>%
    rename("y" = C, "q" = Q10, "z" = Tth, "w" = d)

percs <- c(0.1, 0.9)
qpercs <- function(x) quantile(x, percs)

# ≬(x, y) = y[1] < x < y[2]

# plants
plant_dat <- dat %>%
    filter(Type == "biochem") %>%
    select(q, y, z, w)

plant_percs <- lapply(plant_dat, qpercs)

plant_vars <- plant_dat %>%
    filter(
        q > plant_percs$q[1], q < plant_percs$q[2],
        y > plant_percs$y[1], y < plant_percs$y[2],
        z > plant_percs$z[1], z < plant_percs$z[2],
        w > plant_percs$w[1], w < plant_percs$w[2])

# insects
insect_dat <- dat %>%
    filter(Type == "fit") %>%
    select(q, y, z, w)

insect_percs <- lapply(insect_dat, qpercs)

insect_vars <- insect_dat %>%
    filter(
        q > insect_percs$q[1], q < insect_percs$q[2],
        y > insect_percs$y[1], y < insect_percs$y[2],
        z > insect_percs$z[1], z < insect_percs$z[2],
        w > insect_percs$w[1], w < insect_percs$w[2])

# compute

plant_means <- data.frame(lapply(plant_vars, mean))

plant_covs <- true_cov(plant_vars)
rownames(plant_covs) <- colnames(plant_covs) <- names(plant_vars)
plant_covs <- data.frame(plant_covs)

write.csv(plant_means, "data/plant-means.csv", row.names = FALSE)
write.csv(plant_covs, "data/plant-covs.csv")

#
insect_means <- data.frame(lapply(insect_vars, mean))

insect_covs <- true_cov(insect_vars)
rownames(insect_covs) <- colnames(insect_covs) <- names(insect_vars)
insect_covs <- data.frame(insect_covs)

write.csv(insect_means, "data/insect-means.csv", row.names = FALSE)
write.csv(insect_covs, "data/insect-covs.csv")

# ==============================================================================
# temperature distribution
# ======================0========================================================

# data from https://worldclim.org/data/worldclim21.html
dir_min <- "data/worldclim-tmin"
dir_avg <- "data/worldclim-tavg"
dir_max <- "data/worldclim-tmax"

# 
worldclim_global("tmin", res = 10, path = dir_min)
worldclim_global("tavg", res = 10, path = dir_avg)
worldclim_global("tmax", res = 10, path = dir_max)

# we need the yearly average as well as the average minimum in the winter and the average maximum in the summer

min_files <- list.files(dir_min, recursive = TRUE, full.names = TRUE)
avg_files <- list.files(dir_avg, recursive = TRUE, full.names = TRUE)
max_files <- list.files(dir_max, recursive = TRUE, full.names = TRUE)

min_stack <- stack(lapply(min_files, raster))
avg_stack <- stack(sapply(avg_files, raster))
max_stack <- stack(sapply(max_files, raster))

min_min <- min(min_stack, na.rm = TRUE)
avg_avg <- mean(avg_stack, na.rm = TRUE)
max_max <- max(max_stack, na.rm = TRUE)

min_lats <- rowSums(min_min, na.rm = TRUE) / rowSums(!is.na(min_min), na.rm = TRUE)
avg_lats <- rowSums(avg_avg, na.rm = TRUE) / rowSums(!is.na(avg_avg), na.rm = TRUE)
max_lats <- rowSums(max_max, na.rm = TRUE) / rowSums(!is.na(max_max), na.rm = TRUE)
lats <- unique(coordinates(min_min)[,2])

dat_by_lat <- data.frame(
    latitude = lats,
    Tmin = min_lats,
    Tavg = avg_lats,
    Tmax = max_lats) %>%
    na.omit

write.csv(dat_by_lat, "data/dat-by-lat.csv", row.names = FALSE)

ggplot(dat_by_lat) +
    aes(x = latitude) +
    geom_path(aes(y = Tmin), color = "blue") +
    geom_smooth(aes(y = Tmin), method = "glm", formula = y ~ poly(x, 2),
        lty = 2, color = "blue") +
    geom_path(aes(y = Tavg), color = "grey") +
    geom_smooth(aes(y = Tavg), method = "glm", formula = y ~ poly(x, 2),
        lty = 2, color = "grey") +
    geom_path(aes(y = Tmax), color = "red") + 
    geom_smooth(aes(y = Tmax), method = "glm", formula = y ~ poly(x, 2),
        lty = 2, color = "red") +
    theme_bw() +
    labs(x = "Latitude", y = "Temperature")

min_mod <- lm(Tmin ~ poly(latitude, 2, raw = TRUE), dat_by_lat)
avg_mod <- lm(Tavg ~ poly(latitude, 2, raw = TRUE), dat_by_lat)
max_mod <- lm(Tmax ~ poly(latitude, 2, raw = TRUE), dat_by_lat)

term_names <- c("intercept", "linear_term", "quadratic_term")

dat_coef <- data.frame(rbind(coef(min_mod), coef(avg_mod), coef(max_mod))) %>%
    mutate(type = c("min", "avg", "max"))
colnames(dat_coef) <- c(term_names, "type")
dat_coef <- relocate(dat_coef, type)

write.csv(dat_coef, "data/model-coef-temperature.csv", row.names = FALSE)
