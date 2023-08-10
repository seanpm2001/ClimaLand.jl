function moving_average(x_ma, x_new, window)

end

function distribution_mean(T_ma::FT, θ_ma::FT, parameters) where {FT}
    (; a, b) = parameters
    # mean must be positive and smaller than some bound
    return FT(1)
end

function distribution_std(T_ma::FT, θ_ma::FT, parameters) where {FT}
    # std must be positive and smaller than some limit
    return FT(1)
end
