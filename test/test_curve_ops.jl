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
        @test Q.pnts == [
            point(0.0, 0.0),
            point(1.0, 0.0),
            point(2.0, 0.0),
            point(2.0, 1.0),
            point(2.0, 2.0),
        ]
        @test length(split_edges(PntSeq2F())) == 0
        @test_throws ArgumentError split_edges(
            PntSeq([Point{2,Int}(1, 2), Point{2,Int}(3, 4)]),
        )
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

            # Mixed separators, tabs, blank lines, and comments.
            m = joinpath(dir, "d.txt")
            write(m, "# comment\n1.0, 2.0\n\n3.0\t4.0\n  5.0 ,6.0  \n\n")
            @test read_plt(m).pnts == [point(1.0, 2.0), point(3.0, 4.0), point(5.0, 6.0)]
            @test length(read_plt(joinpath(dir, "empty.txt") |> f -> (write(f, ""); f))) ==
                  0

            # Type stability across dimensions, UTF-8 whitespace, trailing commas, and dim validation.
            @test (@inferred read_plt(f)).pnts == Q.pnts
            read_plt_3d(p) = read_plt(p; dim = 3)
            @test (@inferred read_plt_3d(h)).pnts == [point(1.0, 2.0, 3.0), point(4.0, 5.0, 6.0)]
            @test (@inferred BasicCompGeometry.ReadWrite._read_plt(h, Val(3))).pnts ==
                  [point(1.0, 2.0, 3.0), point(4.0, 5.0, 6.0)]
            @test_throws ArgumentError read_plt(f; dim = 0)

            extra = joinpath(dir, "extra.txt")
            write(extra, "1.0, 2.0,\n3.0,\u00a04.0 tag\n")
            @test read_plt(extra).pnts == [point(1.0, 2.0), point(3.0, 4.0)]

            # Zero-allocation per-line point parsing.
            @test BasicCompGeometry.ReadWrite._parse_point("1.0, 2.0", Val(2), "1.0, 2.0") ==
                  point(1.0, 2.0)
            @test (@allocated BasicCompGeometry.ReadWrite._parse_point(
                "1.0, 2.0",
                Val(2),
                "1.0, 2.0",
            )) == 0
        end
    end

    @testset "read_geolife_plt" begin
        mktempdir() do dir
            f = joinpath(dir, "t.plt")
            write(
                f,
                """Geolife trajectory
WGS 84
Altitude is in Feet
Reserved 3
0,2,255,My Track,0,0,2,8421376
0
39.9,116.3,0,492,39744.12,2008-10-23,02:53:04
39.9,116.3,0,492,39744.12,2008-10-23,02:53:05
39.95,116.4,0,492,39744.13,2008-10-23,02:53:06
""",
            )
            @test read_geolife_plt(f).pnts == [point(39.9, 116.3), point(39.95, 116.4)]
            write(f, readchomp(f) * "\n\n")      # trailing blank lines are ignored
            @test length(read_geolife_plt(f)) == 2

            two_col = joinpath(dir, "two_col.plt")
            write(two_col, "h1\nh2\nh3\nh4\nh5\nh6\n10.5,20.5\n30.5,40.5\n")
            @test read_geolife_plt(two_col).pnts == [point(10.5, 20.5), point(30.5, 40.5)]

            bad = joinpath(dir, "bad.plt")
            write(bad, "h1\nh2\nh3\nh4\nh5\nh6\n10.5\n")
            @test_throws ErrorException read_geolife_plt(bad)
        end
    end
end
