@testset "P1 points" begin
    for field = [ComplexF64,AcbField(),QQBarField()]
        @info "Field" field
        a = Dict((i,j)=>P1Point(field(i+j*im)) for i=-2:3,j=-3:4)
        b = P1Inf(a[0,0])
        @test b ≎ antipode(a[0,0])
        @test a[1,0] ≎ antipode(a[-1,0])
        @test distance(a[0,1],a[1,0]) ≎ pi/2
        @test 0.456874 < distance(a[1,2],a[3,4]) < 0.456875
        for p=[values(a)...; b]
            @test P1Point_R3(R3_P1Point(p)) ≎ p
        end
        for p=values(a), q=[a[0,-1],b,a[1,2]]
            p==antipode(q) && continue # midpoint impossible for antipodes
            m = midpoint(p,q)
            @test distance(m,p) ≎ distance(q,m) ≎ distance(p,q)/2
        end
        for p=values(a)
            @test x_ratio(b,a[0,0],a[1,0],p) == p
        end
        for (p,q,r)=[(a[0,0],a[1,0],a[0,1]),(a[1,0],a[-1,-3],a[-1,3])]
            m,rad = circumcircle(p,q,r)
            @test distance(m,p) ≎ distance(m,q) ≎ distance(m,r) ≎ rad
        end
        @test distance(barycentre(values(a) |> collect), P1Point(field((12043 + 8578im)//1000))) < 0.001
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
        @test image(m,P1Inf(field)) ≎ p ≎ antipode(image(m,P1zero(field)))
        s = moebius_map(1/(1-z))
        t = moebius_map(1/z)
        @test s == moebius_map(P1zero(field),P1one(field),P1Inf(field))
        @test t == moebius_map(field(0),field(1),field(1),field(0))
        for p=[p0,p1,p∞,p]
            @test (s∘t)(p) ≎ s(t(p))
            @test preimage(s,s(p)) ≎ s(preimage(s,p)) ≎ p
        end

        f = p1_map(s)^2 / p1_map(t)^3 # z^3 / (1-z)^2
        @test preimages(f,p0) == [p0 => 3]
        for (r,i)=preimages(f,p∞)
            @test (r ≎ p1 && i==2) || (r ≎ p∞ && i==1)
        end
        for (r,i)=preimages(f,p1)
            @test f(r) ≎ p1 && i==1
        end
        for (r,i)=critical_points(f)
            @test (r ≎ p0 && i==2) || (r ≎ p1 && i==1) || (r ≎ P1Point(field(3)) && i==1)
        end

        ff = p1_map([p∞,p1=>2],[p0=>3],P1Point(field(2))=>P1Point(field(8)))
        for z=[p∞,p0,p1,P1Point(field(3))]
            @test f(z) ≎ ff(z)
        end

        @test derivative(f,p∞) ≎ field(1)
        @test derivative(f,p1) == derivative(f,p0) == field(0)
        let (z,dz) = image_df(f,P1Point(field(2)))
            @test z ≎ P1Point(field(8))
            @test dz ≎ field(-4//13)
        end
        o = p1_map(z)
        df = derivative(f) # z^2*(3-z)/(1-z)^3
        
        @test df(P1Point(field(3))) ≎ p0
        @test df(p∞) ≎ p1
        @test integral(o^3+1) == o^4/4 + o
        
        # p1_intersect

        # rotating_moebius_map
        # normalizing_moebius_map
    end
end
