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
println("Saved TPC distributions")

# ==============================================================================
# temperature distribution
# ==============================================================================

using ArchGDAL, Rasters

# data from https://worldclim.org/data/worldclim21.html
dir_min = "data/wc2.1_10m_tmin"
dir_avg = "data/wc2.1_10m_tavg"
dir_max = "data/wc2.1_10m_tmax"

# we need the yearly average as well as the average minimum in the winter and the average maximum in the summer
min_files = readdir(dir_min, join = true)
avg_files = readdir(dir_avg, join = true)
max_files = readdir(dir_max, join = true)

min_stack = Raster.(min_files)
avg_stack = Raster.(avg_files)
max_stack = Raster.(max_files)

min_stack_summary = mosaic(minimum, min_stack)
avg_stack_summary = mosaic(mean, avg_stack)
max_stack_summary = mosaic(maximum, max_stack)

using GLM, CairoMakie

function lat_means(x)
    nan = x.missingval
    lats = eachcol(x.data)
    fixed_lats = filter.(x -> x != nan, lats)
    return mean.(fixed_lats)
end

lats = dims(min_stack_summary)[2].val.data
Tmin_lats = lat_means(min_stack_summary)
Tavg_lats = lat_means(avg_stack_summary)
Tmax_lats = lat_means(max_stack_summary)

temperature = DataFrame(
    latitude = lats,
    Tmin = Float64.(Tmin_lats),
    Tavg = Float64.(Tavg_lats),
    Tmax = Float64.(Tmax_lats))

subset!(temperature, (names(temperature) .=> ByRow(x -> !isnan(x)))...)

min_mod = lm(@formula(Tmin ~ latitude + latitude^2), temperature)
avg_mod = lm(@formula(Tavg ~ latitude + latitude^2), temperature)
max_mod = lm(@formula(Tmax ~ latitude + latitude^2), temperature)

begin
    fig = Figure(size = (1000, 400))
    ax = Axis(fig[1,1], ylabel = "Latitude", xlabel = "Longitude")
    plot!(avg_stack_summary)
    ax2 = Axis(fig[1,2], xlabel = "Temperature", yticklabelsvisible = false)
    ylims!(extrema(lats)...)
    lmin = lines!(Tmin_lats, lats, color = :blue)
    lavg = lines!(Tavg_lats, lats, color = :grey)
    lmax = lines!(Tmax_lats, lats, color = :orange)
    lines!(Float64.(predict(min_mod, DataFrame(latitude = lats))),
        lats, color = :blue, linestyle = :dash)
    lines!(Float64.(predict(avg_mod, DataFrame(latitude = lats))),
        lats, color = :grey, linestyle = :dash)
    lines!(Float64.(predict(max_mod, DataFrame(latitude = lats))),
        lats, color = :orange, linestyle = :dash)
    axislegend(ax2, [lmin, lavg, lmax], ["Min", "Average", "Max"]; position = :lt)
    fig
end

save("plots/temperature-latitude.png", fig)

term_names = (:intercept, :linear_term, :quadratic_term)
mods = Dict("min" => NamedTuple{term_names}(coef(min_mod)),
    "max" => NamedTuple{term_names}(coef(max_mod)))

save("data/latitude_T.jld2", mods)
println("Saved latitude coefficients")
