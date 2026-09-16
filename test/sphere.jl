@testset "sphere groups" begin
    # (2,3,7) triangle group orbifold: generators x1,x2,x3 with x3*x2*x1=1
    # (ordering [3,2,1], i.e. sphere_group(3)'s default ordering),
    # x1^2=x2^3=x3^7=1. Classic hyperbolic triangle group.
    g = sphere_group(3, [2,3,7])
    @test sphere_rank(g) == 3
    @test sphere_ordering(g) == [3,2,1]
    @test sphere_exponents(g) == [2,3,7]
    @test euler_characteristic(g) == -1//42 # -1 + 1/2 + 1/3 + 1/7, halved

    x1,x2,x3 = gens(g)
    @test x1^2 == one(g)
    @test x2^3 == one(g)
    @test x3^7 == one(g)
    @test x3*x2*x1 == one(g)
    @test order(x1) == 2 && order(x2) == 3 && order(x3) == 7

    peri = peripheral_classes(g)
    @test length(peri) == 3
    @test all(is_peripheral, peri)
    @test is_peripheral(x1)

    # (2,4,4) is a Euclidean (Euler characteristic 0) triangle group
    g0 = sphere_group(3, [2,4,4])
    @test euler_characteristic(g0) == 0

    # explicit ordering vector, no exponents (n=4 generators, isomorphic to a
    # free group of rank n-1=3, Euler characteristic 1-3=-2)
    gf = sphere_group([2,1,4,3])
    @test sphere_rank(gf) == 4
    @test sphere_exponents(gf) == [0,0,0,0]
    @test euler_characteristic(gf) == 1 - 3

    iso = isomorphism_free_group(gf)
    @test iso !== nothing

    # intersection numbers: two distinct peripheral (generator) classes of a
    # sphere group are disjoint simple curves around different punctures
    cx1, cx2 = conjugacy_class(g, x1), conjugacy_class(g, x2)
    @test intersection_number(cx1, cx2) == 0
    @test self_intersection_number(cx1) == 0

    # as_sphere_group / isomorphism_sphere_group round trip via a plain f.p. group
    F = free_group(3)
    f1,f2,f3 = gens(F)
    G, _ = quo(F, [f1^2, f2^3, f3^7, f3*f2*f1])
    g2 = as_sphere_group(G)
    @test g2 !== nothing
    @test sphere_exponents(g2) == [2,3,7]

    iso2 = isomorphism_sphere_group(G)
    @test iso2 !== nothing

    # automorphism group / epimorphism to Out: smoke test only (complex GAP machinery)
    a = sphere_automorphism_group(g)
    @test a !== nothing
    eo = epimorphism_to_out(a)
    @test eo !== nothing

    # amalgamated free product: needs infinite-order generators, so use gf
    gf2 = sphere_group([2,1,4,3])
    amal = amalgamated_free_product(gf, gf2, gens(gf)[1], gens(gf2)[1])
    @test amal isa SphereGroup
    embs = embeddings_of_amalgamated_free_product(amal)
    @test length(embs) == 2
end
