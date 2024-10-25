using JLD2, Distributions, CSV, DataFrames

# ==============================================================================
# TPC
# ==============================================================================

function true_cov(mat)
    n = size(mat)[2]
    M = zeros(n, n)
    for i in 1:n
        for j in 1:n
            M[i,j] = mean((mat[:,i] .- mean(mat[:,i])) .* (mat[:,j] .- mean(mat[:,j])))
        end
    end
    # rownames(M) <- colnames(M) <- colnames(mat)
    return M
end

dat = CSV.read("data/rezende2019-supp.csv", DataFrame)
filter!(:Type => x -> x != "run", dat)
select!(dat, :Type, :Q10, :C, :Tth, :d)
transform!(dat, :C => ByRow(C -> log(C)) => :C)
transform!(dat, :d => ByRow(d -> log(d)) => :d)
rename!(dat, :C => "y", :Q10 => "q", :Tth => "z", :d => "w")

percs = (0.01, 0.99)
qpercs(x) = quantile(x, percs)
≬(x, y) = y[1] < x < y[2]

# plants
plant_dat = filter(:Type => x -> x == "biochem", dat)
select!(plant_dat, :q, :y, :z, :w)

par_percs = NamedTuple{Tuple(Symbol.(names(plant_dat)))}(
    (qpercs(plant_dat[!, par]) for par in Symbol.(names(plant_dat)))
    )

plant_vars = plant_dat[
    .≬(plant_dat.q, Ref(par_percs.q)) .& .≬(plant_dat.y, Ref(par_percs.y)) .&
    .≬(plant_dat.z, Ref(par_percs.z)) .& .≬(plant_dat.w, Ref(par_percs.w)),:
]

# insects
insect_dat = filter(:Type => x -> x == "fit", dat)
select!(insect_dat, :q, :y, :z, :w)

par_percs = NamedTuple{Tuple(Symbol.(names(insect_dat)))}(
    (qpercs(insect_dat[!, par]) for par in Symbol.(names(insect_dat)))
    )

insect_vars = insect_dat[
    .≬(insect_dat.q, Ref(par_percs.q)) .& .≬(insect_dat.y, Ref(par_percs.y)) .&
    .≬(insect_dat.z, Ref(par_percs.z)) .& .≬(insect_dat.w, Ref(par_percs.w)),:
]


# compute

ret =  (
    #
    plant_names = names(plant_vars),
    plant_means = mean.(eachcol(plant_vars)),
    plant_covs = true_cov(plant_vars),
    #
    insect_names = names(insect_vars),
    insect_means = mean.(eachcol(insect_vars)),
    insect_covs = true_cov(insect_vars)
)

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
