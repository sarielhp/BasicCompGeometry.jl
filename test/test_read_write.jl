using Test
using BasicCompGeometry
using BasicCompGeometry.ReadWrite
using LinearAlgebra

@testset "ReadWrite" begin
    @testset "fbin" begin
        filename = tempname() * ".fbin"
        data = rand(Float32, 5, 3) # 5 vectors of dimension 3 (if we follow write_fbin logic)
        # Actually write_fbin uses size(vecs) -> nvecs, dim
        # So it writes vecs' (dim x nvecs)

        write_fbin(filename, data)

        # read_fbin returns (dim, nvecs)
        read_data = read_fbin(filename)
        @test size(read_data) == (3, 5)
        @test read_data ≈ data'

        rm(filename)
    end

    @testset "ibin" begin
        filename = tempname() * ".ibin"
        data = rand(Int32(1):Int32(100), 5, 3)
        write_ibin(filename, data)

        # read_ibin returns matrix' -> (nvecs, dim)
        read_data = read_ibin(filename)
        @test size(read_data) == (5, 3)
        @test read_data == data

        rm(filename)
    end

    @testset "fvecs" begin
        filename = tempname() * ".fvecs"
        dim = Int32(4)
        nvecs = 3
        data = rand(Float32, dim, nvecs)

        open(filename, "w") do io
            for i = 1:nvecs
                write(io, dim)
                write(io, data[:, i])
            end
        end

        read_data = read_fvecs(filename)
        @test size(read_data) == (dim, nvecs)
        @test read_data ≈ data

        rm(filename)
    end
end

@testset "Binary subsets and malformed fvecs" begin
    mktempdir(mkpath(joinpath(@__DIR__, "..", "output", "tmp"))) do dir
        floats = Float32[1 2; 3 4; 5 6]
        fbin = joinpath(dir, "subset.fbin")
        write_fbin(fbin, floats)
        @test read_fbin_n(fbin) == floats'
        @test read_fbin_n(fbin, 2) == floats[1:2, :]'
        @test read_fbin_n(fbin, 20) == floats'
        @test read_fbin(fbin, 1, 1) == floats[2:2, :]'
        @test read_fbin(fbin, 1) == floats[2:3, :]'
        ints = Int32.(floats)
        ibin = joinpath(dir, "subset.ibin")
        write_ibin(ibin, ints)
        @test read_ibin(ibin, 1, 1) == ints[2:2, :]

        malformed = joinpath(dir, "malformed.fvecs")
        write(malformed, UInt8[])
        @test_throws ErrorException read_fvecs(malformed)
        write(malformed, Int32(0))
        @test_throws ErrorException read_fvecs(malformed)
        write(malformed, Int32[2, 0]) # Incomplete record.
        @test_throws ErrorException read_fvecs(malformed)
        open(malformed, "w") do io
            write(io, Int32(2), Float32[1, 2], Int32(3), Float32[3, 4])
        end
        @test_throws ErrorException read_fvecs(malformed)
        @test_throws ErrorException read_fvecs(joinpath(dir, "missing.fvecs"))
    end
end
