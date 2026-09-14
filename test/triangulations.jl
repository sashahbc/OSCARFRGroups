@testset "triangulations: exact regression counts (img/tst/chapter-9-b.tst)" begin
    for field = [AcbField(128), QQBarField()]
        @info "Field" field
        Fr = field isa QQBarField ? field : ArbField(128)
        s = sqrt(Fr(1)//3)
        cube_pts = [(s,s,s),(s,s,-s),(s,-s,s),(-s,s,s),(s,-s,-s),(-s,s,-s),(-s,-s,s),(-s,-s,-s)]
        cube = [P1Point_R3(p) for p in cube_pts]

        # DelaunayTriangulation(cube{[1,5]}) -> bare octahedron
        t0 = delaunay_triangulation(cube[[1,5]])
        @test length(t0.v)==6 && length(t0.e)==24 && length(t0.f)==8
        @test check_triangulation(t0) === true

        # DelaunayTriangulation(cube) -> 11 vertices, 54 edges, 18 faces
        t1 = delaunay_triangulation(cube)
        @test length(t1.v)==11 && length(t1.e)==54 && length(t1.f)==18
        @test check_triangulation(t1) === true
    end

    # DelaunayTriangulation(p,100.) on 4 near-degenerate points -> 32/180/60
    field = AcbField(128)
    Fr = ArbField(128)
    p1 = P1Point_R3((Fr(0),Fr(0),Fr(1)))
    p2 = P1Point_R3((Fr(0),Fr(0),Fr(-1)))
    p3 = P1Point(field(1//10^4))
    p4 = P1Point(field(1//10^4) + field(10^4)*field(im))
    t = delaunay_triangulation([p1,p2,p3,p4]; quality=100.0)
    @test length(t.v)==32 && length(t.e)==180 && length(t.f)==60
    @test check_triangulation(t) === true
end

# The broader mutation smoke test below uses AcbField only. QQBarField's exact
# algebraic arithmetic can blow up unpredictably (very slow, though not
# incorrect) through the accumulated sqrt/circumcircle operations of many
# incremental insertions and flips -- a known general characteristic of exact
# algebraic-number computation, not a correctness issue (confirmed by the
# fixed-input regressions above, which do use QQBarField).
@testset "triangulations: incremental operations" begin
    field = AcbField(128)

    pts = equidistributed_p1_points(20, field)
    @test length(pts) == 20
    t = delaunay_triangulation(pts)
    @test check_triangulation(t) === true

    v_before = length(t.v)
    p_new = centre(t.f[1])
    add_to_triangulation!(t, p_new; delaunay=true)
    @test length(t.v) == v_before+1
    @test check_triangulation(t) === true

    # removal only guarantees a structurally valid triangulation, not global
    # Delaunay-ness (it locally minimizes circumcircle power while flipping
    # the vertex down to trivalent, matching GAP's RemoveFromTriangulation) --
    # explicit fix_delaunay! restores full Delaunay-ness.
    remove_from_triangulation!(t, t.v[end])
    @test length(t.v) == v_before
    @test check_triangulation(t; nodelaunay=true) === true
    fix_delaunay!(t)
    @test check_triangulation(t) === true

    # a single non-multi flip breaks Delaunay-ness in general, but the DCEL
    # itself must stay structurally valid; fix_delaunay! restores it.
    flip_edge!(t.e[1], false)
    @test check_triangulation(t; nodelaunay=true) === true
    fix_delaunay!(t)
    @test check_triangulation(t) === true

    # wiggled_triangulation via a Möbius map, and a bare shallow_copy
    m = moebius_map(pos(t.v[1]), pos(t.v[2]))
    t2 = wiggled_triangulation(t, m)
    @test check_triangulation(t2) === true
    @test length(t2.v)==length(t.v) && length(t2.e)==length(t.e) && length(t2.f)==length(t.f)

    t3 = shallow_copy(t)
    @test check_triangulation(t3) === true

    # closest_* dispatches
    v1 = t.v[1]
    @test closest_vertex(v1) === v1
    @test v1 in closest_vertices(closest_face(v1))
    @test closest_face(closest_face(v1)) === closest_face(v1)
end
