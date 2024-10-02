# todo stop this being a csv maybe? this is gross

mean_dat = CSV.read("data/rezende2019-plant-means.csv", DataFrame)
means = mean_dat.x
covs = Matrix(CSV.read("data/rezende2019-plant-covs.csv", DataFrame)[1:4,2:5])

save("data/T_dist.jld2",
    Dict(
        "T_dist" => MvNormal(means, covs),
        "names" => Tuple(Symbol.(mean_dat.Column1))
        )
    )

# todo include r script from playground that shows where this all came from
