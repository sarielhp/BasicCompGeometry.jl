using BasicCompGeometry
using BasicCompGeometry: push_distinct!, split_edges, write_plt
using BasicCompGeometry.ReadWrite
using Test

@testset "Curve operations and text readers" begin
    @testset "push_distinct!" begin
        P = PntSeq2F()
        @test push_distinct!(P, point(1.0, 2.0))
        @test !push_distinct!(P, point(1.0, 2.0))
        @test push_distinct!(P, point(1.0, 2.5))
        @test !push_distinct!(P, point(1.0, 2.5))
        @test push_distinct!(P, point(1.0, 2.0))   # only the last vertex is compared
        @test length(P) == 3
    end

    @testset "split_edges" begin
        P = PntSeq([point(0.0, 0.0), point(2.0, 0.0), point(2.0, 2.0)])
        Q = split_edges(P)
        @test Q.pnts == [point(0.0, 0.0), point(1.0, 0.0), point(2.0, 0.0), point(2.0, 1.0), point(2.0, 2.0)]
        @test length(split_edges(PntSeq2F())) == 0
        @test split_edges(PntSeq([point(1.0, 1.0)])).pnts == [point(1.0, 1.0)]
        # Repeated points collapse.
        R = split_edges(PntSeq([point(0.0, 0.0), point(0.0, 0.0), point(1.0, 0.0)]))
        @test R.pnts == [point(0.0, 0.0), point(0.5, 0.0), point(1.0, 0.0)]
    end

    @testset "reverse" begin
        P = PntSeq([point(0.0, 0.0), point(1.0, 1.0), point(2.0, 0.0)])
        @test reverse(P).pnts == reverse(P.pnts)
        @test P.pnts[1] == point(0.0, 0.0)          # input untouched
    end

    @testset "read_plt" begin
        mktempdir() do dir
            f = joinpath(dir, "a.plt")
            P = PntSeq([point(0.0, 0.0), point(1.5, 2.5), point(3.0, 4.0)])
            write_plt(f, P)
            Q = read_plt(f)
            @test Q.pnts == P.pnts

            g = joinpath(dir, "b.txt")
            write(g, "x y k tid\n1 2 7 8\n1 2 7 8\n3 4 7 8\n")
            @test read_plt(g).pnts == [point(1.0, 2.0), point(3.0, 4.0)]

            h = joinpath(dir, "c.plt")
            write_plt(h, PntSeq([point(1.0, 2.0, 3.0), point(4.0, 5.0, 6.0)]))
            @test read_plt(h; dim = 3).pnts == [point(1.0, 2.0, 3.0), point(4.0, 5.0, 6.0)]
            @test_throws ErrorException read_plt(h; dim = 4)
        end
    end

    @testset "read_geolife_plt" begin
        mktempdir() do dir
            f = joinpath(dir, "t.plt")
            write(f, """Geolife trajectory
WGS 84
Altitude is in Feet
Reserved 3
0,2,255,My Track,0,0,2,8421376
0
39.9,116.3,0,492,39744.12,2008-10-23,02:53:04
39.9,116.3,0,492,39744.12,2008-10-23,02:53:05
39.95,116.4,0,492,39744.13,2008-10-23,02:53:06
""")
            @test read_geolife_plt(f).pnts == [point(39.9, 116.3), point(39.95, 116.4)]
        end
    end
end
