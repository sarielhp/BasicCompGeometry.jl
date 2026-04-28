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
