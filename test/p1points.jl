@testset "P1 points" begin
    for field = [ComplexF64,AcbField(),QQBarField()]
        @info "Field" field
        a = Dict((i,j)=>P1Point(field(i+j*im)) for i=-2:3,j=-3:4)
        b = P1Inf(a[0,0])
        @test proj_equiv(b, antipode(a[0,0]))
        @test proj_equiv(a[1,0], antipode(a[-1,0]))
        @test proj_equiv(distance(a[0,1],a[1,0]), pi/2)
        @test 0.456874 < distance(a[1,2],a[3,4]) < 0.456875
        for p=[values(a)...; b]
            @test proj_equiv(P1Point_R3(R3_P1Point(p)), p)
        end
        for p=values(a), q=[a[0,-1],b,a[1,2]]
            p==antipode(q) && continue # midpoint impossible for antipodes
            m = midpoint(p,q)
            @test proj_equiv(distance(m,p),distance(q,m)) && proj_equiv(distance(q,m),distance(p,q)/2)
        end
        for p=values(a)
            @test x_ratio(b,a[0,0],a[1,0],p) == p
            @test P1Point(cross_ratio(b,a[0,0],a[1,0],p)) == p
        end
        for (p,q,r)=[(a[0,0],a[1,0],a[0,1]),(a[1,0],a[-1,-3],a[-1,3])]
            m,rad = circumcircle(p,q,r)
            @test proj_equiv(distance(m,p),distance(m,q)) && proj_equiv(distance(m,q),distance(m,r)) && proj_equiv(distance(m,r),rad)
        end
        @test distance(barycentre(values(a) |> collect), P1Point(field((12043 + 8578im)//1000))) < 0.001

        @test proj_equiv(cleaned_p1point(P1Point(field(1//10^9))), P1zero(field))
        @test proj_equiv(cleaned_p1point(b), b)
        @test proj_equiv(cleaned_p1point(a[1,0]), a[1,0])

        near = P1Point(field(1) + field(1//10^9))
        cl = collected_p1_points([a[1,0],near,a[2,0]])
        @test length(cl) == 2
        @test any(x -> x[2]==2, cl)

        cands = [a[2,0],a[1,0],a[-1,0]]
        @test match_p1_points([a[1,0]],[cands]) == [2]
        @test proj_equiv(closest_p1_point(cands,a[1,0]), a[1,0])
        @test match_p1_points([a[0,0]],[[a[1,0],a[-1,0]]]) === nothing # equidistant, ambiguous
    end
end

@testset "P1 maps" begin
    for field = [AcbField(),QQBarField()]
        @info "Field" field
        R,z = rational_function_field(field,:z)
        zn = numerator(z)

        @test p1_map(zn) == p1_map(z) == p1_map(moebius_map(zn)) == p1_map(moebius_map(z)) == gen(p1_map(z))

        p = P1Point(field(2+3im))
        p0 = P1zero(field)
        p1 = P1one(field)
        p∞ = P1Inf(field)
        m = moebius_map(p)
        @test proj_equiv(image(m,P1Inf(field)),p) && proj_equiv(p, antipode(image(m,P1zero(field))))
        s = moebius_map(1/(1-z))
        t = moebius_map(1/z)
        @test s == moebius_map(P1zero(field),P1one(field),P1Inf(field))
        @test t == moebius_map(field(0),field(1),field(1),field(0))
        for p=[p0,p1,p∞,p]
            @test proj_equiv((s∘t)(p), s(t(p)))
            @test proj_equiv(preimage(s,s(p)), s(preimage(s,p))) && proj_equiv(s(preimage(s,p)), p)
        end

        f = p1_map(s)^2 / p1_map(t)^3 # z^3 / (1-z)^2
        @test preimages(f,p0) == [p0 => 3]
        for (r,i)=preimages(f,p∞)
            @test (proj_equiv(r,p1) && i==2) || (proj_equiv(r,p∞) && i==1)
        end
        for (r,i)=preimages(f,p1)
            @test proj_equiv(f(r),p1) && i==1
        end
        for (r,i)=critical_points(f)
            @test (proj_equiv(r,p0) && i==2) || (proj_equiv(r,p1) && i==1) || (proj_equiv(r,P1Point(field(3))) && i==1)
        end

        ff = p1_map([p∞,p1=>2],[p0=>3],P1Point(field(2))=>P1Point(field(8)))
        for z=[p∞,p0,p1,P1Point(field(3))]
            @test proj_equiv(f(z), ff(z))
        end

        @test proj_equiv(derivative(f,p∞), field(1))
        @test derivative(f,p1) == derivative(f,p0) == field(0)
        let (z,dz) = image_df(f,P1Point(field(2)))
            @test proj_equiv(z, P1Point(field(8)))
            @test proj_equiv(dz, field(-4//13))
        end
        o = p1_map(z)
        df = derivative(f) # z^2*(3-z)/(1-z)^3

        @test proj_equiv(df(P1Point(field(3))), p0)
        @test proj_equiv(df(p∞), p1)
        @test integral(o^3+1) == o^4/4 + o

        # p1_map_by_coefficients / coefficients_of_p1_map
        m2 = p1_map_by_coefficients([field(1),field(2),field(3)],[field(1),field(0),field(1)]) # (1+2z+3z^2)/(1+z^2)
        cn,cd,d = coefficients_of_p1_map(m2)
        @test cn == [field(1),field(2),field(3)] && cd == [field(1),field(0),field(1)] && d == 2

        # p1_monomial
        @test proj_equiv(p1_monomial(field,3)(P1Point(field(2))), P1Point(field(8)))
        @test proj_equiv(p1_monomial(field,-2)(P1Point(field(2))), P1Point(field(1//4)))

        # sl2_p1_map / moebius_map(::Matrix) round trip
        mm = moebius_map(field(2),field(3),field(1),field(5))
        @test moebius_map(sl2_p1_map(mm)) == mm

        # moebius_map(Vector) forms
        pts3 = [p0,p1,p∞]
        @test moebius_map(pts3) == moebius_map(pts3...)
        dst3 = [P1Point(field(5)),P1Point(field(6)),P1Point(field(7))]
        msd = moebius_map(pts3,dst3)
        for i=1:3
            @test proj_equiv(msd(pts3[i]), dst3[i])
        end

        # Base.conj(MoebiusP1Map)
        mc = conj(mm)
        @test mc.a==conj(mm.a) && mc.b==conj(mm.b) && mc.c==conj(mm.c) && mc.d==conj(mm.d)

        # p1_map_scaling agrees with image_df's second (derivative-like) output
        @test proj_equiv(p1_map_scaling(f,P1Point(field(2))), abs(image_df(f,P1Point(field(2)))[2]))

        # cleaned_p1map
        noisy = p1_map_by_coefficients([field(1)+field(1//10^12),field(2)],[field(1)])
        @test cleaned_p1map(noisy) == p1_map_by_coefficients([field(1),field(2)],[field(1)])

        # p1_intersect: two Möbius paths meeting at a known point (t,u) = (1, 0.37)
        delta = moebius_path(P1Point(field(1-2im)), P1Point(field(-3+im)))
        u0 = field(37//100)
        target = f(delta(P1Point(u0)))
        gamma = moebius_path(P1Point(field(7-3im)), target) # γ(1) = target = δ(u0)
        res = p1_intersect(gamma,f,delta)
        @test any(r -> proj_equiv(field(r.t),field(1)) && proj_equiv(field(r.u),u0), res)
        for r=res
            @test proj_equiv(gamma(P1Point(field(r.t))), f(delta(P1Point(field(r.u)))))
            @test proj_equiv(r.γt, gamma(P1Point(field(r.t))))
            @test proj_equiv(r.δu, delta(P1Point(field(r.u))))
        end

        # rotating_moebius_map: sends the last point to ∞
        pts = [p0,p1,P1Point(field(2+3im)),P1Point(field(-1+im)),P1Point(field(5-2im))]
        rot = rotating_moebius_map(pts)
        @test proj_equiv(image(rot,pts[end]), p∞)

        # normalizing_moebius_map: also sends the last point to ∞
        norm_map = normalizing_moebius_map(pts)
        @test proj_equiv(image(norm_map,pts[end]), p∞)
        newpts = norm_map.(pts)
        @test proj_equiv(newpts[end], p∞)
    end
end
