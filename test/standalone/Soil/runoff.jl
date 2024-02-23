using ClimaLand
using Test

FT = Float32
@testset "Base runoff functionality, FT = $FT" begin
    runoff = ClimaLand.Soil.NoRunoff()
    precip = 5.0
    @test ClimaLand.Soil.soil_surface_infiltration(runoff, precip) == precip
    @test ClimaLand.Soil.subsurface_runoff_source(runoff) == nothing
    struct Foo{FT} <: ClimaLand.Soil.AbstractSoilSource{FT} end
    srcs = (1, 2, 3)
    @test ClimaLand.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLand.Soil.append_source(Foo{FT}(), srcs) == (srcs..., Foo{FT}())
end

@testset "TOPMODEL runoff FT =$FT" begin

    domain = ClimaLand.Domains.SphericalShell(;
                                              radius = FT(2.0),
                                              depth = FT(1.0),
                                              nelements = (10, 10),
                                              npolynomial = 1,
                                              dz_tuple = FT.((1.0, 0.05)),
                                              );
    surface_space = domain.space.surface
    f_max = ClimaCore.Fields.ones(surface_space)
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4/1000) # m/s

    subsurface_source = ClimaLand.Soil.Runoff.TOPMODELSubsurfaceRunoff{FT}(R_sb)
    runoff_model = TOPMODELRunoff{FT}(f_over, f_max, subsurface_source)
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    soil_parameters = (; ν = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 0.5, Ksat = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-6)
    net_water_flux = ClimaCore.Fields.zeros(domain.space.subsurface) .- 2e-6
    Y = (; soil = )

    #etc
    
    
end

