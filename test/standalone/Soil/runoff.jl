using ClimaLand
using Test
using ClimaCore
FT = Float32
@testset "Base runoff functionality, FT = $FT" begin
    runoff = ClimaLand.Soil.Runoff.NoRunoff()
    precip = 5.0
    @test ClimaLand.Soil.Runoff.soil_surface_infiltration(runoff, precip) == precip
    @test ClimaLand.Soil.Runoff.subsurface_runoff_source(runoff) == nothing
    struct Foo{FT} <: ClimaLand.Soil.Runoff.AbstractSoilSource{FT} end
    srcs = (1, 2, 3)
    @test ClimaLand.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLand.Soil.append_source(Foo{FT}(), srcs) == (srcs..., Foo{FT}())
end

@testset "TOPMODEL runoff FT =$FT" begin

    domain = ClimaLand.Domains.SphericalShell(;
                                              radius = FT(100.0),
                                              depth = FT(10.0),
                                              nelements = (10, 10),
                                              npolynomial = 1,
                                              dz_tuple = FT.((1.0, 0.05)),
                                              );
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface
    f_max = ClimaCore.Fields.ones(surface_space).*FT(0.5)
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4/1000) # m/s

    subsurface_source = ClimaLand.Soil.Runoff.TOPMODELSubsurfaceRunoff{FT}(R_sb)
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT, typeof(f_max)}(f_over, f_max, subsurface_source)
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    vg_α = ClimaCore.Fields.ones(subsurface_space) .* FT(0.2)
    vg_n = ClimaCore.Fields.ones(subsurface_space) .* FT(2.2)
    hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(;α = vg_α, n = vg_n)
    θ_r = ClimaCore.Fields.zeros(subsurface_space)
    
    soil_parameters = (; ν = ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.6), Ksat = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-6), θ_r = θ_r, hydrology_cm = hydrology_cm)
    precip = ClimaCore.Fields.zeros(domain.space.subsurface) .- FT(2e-6)
    ϑ_l = FT(0.5) .- FT(0.3/10) .*(z.+FT(10))
    θ_i = FT(0.1) .+ClimaCore.Fields.zeros(axes(z))
    T = FT(273.5) .+ClimaCore.Fields.zeros(axes(z))
    Y = (; soil = (;ϑ_l = ϑ_l, θ_i = θ_i))
    p = (; soil = (;T = T,))

    
    
end

