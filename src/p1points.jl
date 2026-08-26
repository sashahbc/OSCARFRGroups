################################################################
# points on the sphere

export P1Point, P1Inf, P1zero, P1one, ≎, P1Point_R3, R3_P1Point
export distance, antipode, barycentre, midpoint, x_ratio, circumcircle

export p1_map, moebius_map, p1_path, ℂ, ℂ⁽ᶻ⁾, 𝓏
export preimages, critical_points, periodic_points, fixed_points, image_df

frexp_2(x::Complex) = max(frexp(real(x))[2],frexp(imag(x))[2])
frexp_2(x::AcbFieldElem) = max(x.real_mid_exp,x.imag_mid_exp)
frexp_2(x::QQBarFieldElem) = 1 # no scaling required
Base.ldexp(x::Complex,i) = typeof(x)(ldexp(real(x),i),ldexp(imag(x),i))
Base.ldexp(x::QQBarFieldElem,i) = x*2^i
Base.abs2(x::AcbFieldElem)::ArbFieldElem = real(x)*real(x)+imag(x)*imag(x)
≎(x,y) = x ≈ y # default
≎(x::Union{ArbFieldElem,AcbFieldElem},y::Union{ArbFieldElem,AcbFieldElem}) = overlaps(x,y)
≎(x::Tuple,y::Tuple...) = all(splat(≎),zip(x,y...))
≎(x::AbstractVector,y::AbstractVector...) = all(splat(≎),zip(x,y...))

# P1 points over type T, represented as quotients num/denom
# with the convention that either |num|<=1 and denom=1, or num=1 and |denom|<=1
struct P1Point{T}
    num::T
    den::T
    function P1Point{T}(num::T,den::T=one(num)) where T
        e = 1-max(frexp_2(num),frexp_2(den))
        new(ldexp(num,e),ldexp(den,e))
    end
end
field(x::P1Point{T}) where T <: Complex = T
field(x::P1Point{T}) where T <: FieldElem = parent(x.num)
real_field(x::P1Point{T}) where T <: Complex = real(T)
real_field(x::P1Point{T}) where T <: FieldElem = parent(real(x.num))
Base.real(f::Field) = parent(real(zero(f)))

real_field(x::T) where T = real(T)
real_field(x::AcbFieldElem) = parent(real(x))
real_field(x::QQBarFieldElem) = parent(x)
complex_field(x::T) where T = complex(T)
complex_field(x::ArbFieldElem) = AcbField(precision(parent(x)))
complex_field(x::QQBarFieldElem) = parent(x) # same real as complex field
(f::QQBarField)(r::QQBarFieldElem,i::QQBarFieldElem) = r + f(im)*i
                
P1Point(z::T,w::T=one(z)) where T = P1Point{T}(z,w)
P1Point(f::Field,z,w=one(z)) = P1Point{typeof(zero(f))}(f(z),f(w))

P1Inf(z::P1Point{T}) where T = P1Point{T}(one(z.num),zero(z.num))
P1Inf(T::Type = ComplexF64) = P1Point{T}(T(1),T(0))
P1Inf(f::T) where T <: Field = P1Point{typeof(zero(f))}(f(1),f(0))
Base.zero(z::P1Point{T}) where T = P1Point{T}(zero(z.num),one(z.num))
P1zero(T::Type = ComplexF64) = P1Point{T}(T(0),T(1))
P1zero(f::T) where T <: Field = P1Point{typeof(zero(f))}(f(0),f(1))
Base.one(z::P1Point{T}) where T = P1Point{T}(one(z.num),one(z.num))
P1one(T::Type = ComplexF64) = P1Point{T}(T(1),T(1))
P1one(f::T) where T <: Field = P1Point{typeof(zero(f))}(f(1),f(1))

Base.iszero(z::P1Point) = iszero(z.num)
Base.isone(z::P1Point) = z.num==z.den
Base.isinf(z::P1Point) = iszero(z.den)
maybezero(z::P1Point) = overlaps(z,P1zero(z))
maybeone(z::P1Point) = overlaps(z,P1one(z))
maybeinf(z::P1Point) = overlaps(z,P1Inf(z))

zcoord(z::P1Point) = z.num/z.den # north pole → 0
wcoord(z::P1Point) = z.den/z.num # north pole → ∞
xycoord(z::P1Point) = 2z.num*conj(z.den) / (abs2(z.num) + abs2(z.den)) # projection to xy plane
ycoord(z::P1Point) = imag(xycoord(z))

Oscar.overlaps(z::P1Point{AcbFieldElem},w::P1Point{AcbFieldElem}) = overlaps(z.num*w.den, z.den*w.num)
Oscar.overlaps(z::P1Point{T},w::P1Point{T}) where T = z == w

function Base.show(io::IO,z::P1Point)
    if maybeinf(z)
        print(io, "ℙ¹(∞)")
    else
        print(io, "ℙ¹(",zcoord(z),")")
    end
end

Base.:(==)(z::P1Point,w::P1Point) = z.num*w.den == z.den*w.num
≎(z::P1Point,w::P1Point) = z.num*w.den ≎ z.den*w.num
Base.hash(z::P1Point,h) = hash(z.num,hash(z.den,h))

function P1Point_R3(v::NTuple{3,T}) where T
    U = complex_field(v[1])
    P1Point(U(v[1],v[2]),U(sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])+v[3]))
end
    
function R3_P1Point(z::P1Point{T}) where T
    p = 2z.num*conj(z.den)
    n = abs2(z.num)
    d = abs2(z.den)
    (real(p),imag(p),d-n) ./ (d+n)
end

function antipode(z::P1Point{T}) where T
    P1Point{T}(-conj(z.den),conj(z.num))
end

barycentre(z::P1Point) = z
barycentre(z::P1Point,w::P1Point) = midpoint(z,w)
function barycentre(z::Vector{P1Point{T}}) where T
    @assert length(z) ≥ 1
    U = real_field(z[1])
    b = (U(0),U(0),U(0))
    for p=z
        b = b .+ R3_P1Point(p)
    end
    P1Point_R3(b)
end
barycentre(z::P1Point,w::P1Point...) = barycentre([z;w...])

function Oscar.midpoint(z::P1Point{T},w::P1Point{T}) where T
    z==w && return z

    d = w.num*conj(z.num) + w.den*conj(z.den)
    iszero(d) && error("cannot find midpoint of antipodes")
    a = sqrt((abs2(z.num)+abs2(z.den))/(abs2(w.num)+abs2(w.den))*abs2(d))
    P1Point{T}(a*w.num+d*z.num,a*w.den+d*z.den)
end

"""Compute the spherical distance between points z,w. Returns a Float64."""
function distance(z::P1Point{T},w::P1Point{T}) where T
    z==w && return 0.0
    
    d = w.num*conj(z.num) + w.den*conj(z.den)
    if iszero(d)
        Float64(pi)
    else
        v = abs((w.num*z.den-w.den*z.num)/d) # image of w under map z↦0
        2atan(Float64(v))
    end
end

# collected p1 points, gather and cluster

# match p1 points

# closest p1 point

function x_ratio(p₁::P1Point{T},p₂::P1Point{T},p₃::P1Point{T},p₄::P1Point{T}) where T
    P1Point{T}((p₁.num*p₃.den-p₃.num*p₁.den) * (p₂.num*p₄.den-p₄.num*p₂.den),
               (p₂.num*p₃.den-p₃.num*p₂.den) * (p₁.num*p₄.den-p₄.num*p₁.den))
end

function circumcircle(p₁::P1Point{T},p₂::P1Point{T},p₃::P1Point{T}) where T
    F = field(p₁)
    p = zero(F)
    q = zero(F)
    for (a,b,c)=[(p₁,p₂,p₃),(p₂,p₃,p₁),(p₃,p₁,p₂)]
        p += abs2(a.num)*b.den*c.den*conj(b.num*c.den-c.num*b.den)
        q += a.num*b.den*conj(b.num*a.den)*(abs2(c.num)+abs2(c.den))
    end

    qimag = imag(q)
    
    centre = P1Point{T}(F(im)*(-qimag+sqrt(abs2(p)+qimag*qimag)),p)

    d = abs(centre.num*p₁.den - p₁.num*centre.den) / abs(conj(centre.num)*p₁.num + p₁.den*conj(centre.den));
    
    if d > 1
        d = inv(d)
        centre = antipode(centre)
    end
    
    (centre, 2atan(Float64(d)))
end    

################################################################
# rational maps on spheres

abstract type P1Map{T} end

struct RationalP1Map{T,R} <: P1Map{T}
    map::AbstractAlgebra.Generic.RationalFunctionFieldElem{T,R}
end

Oscar.degree(map::RationalP1Map) = max(degree(numerator(map.map)),degree(denominator(map.map)))
Oscar.function_field(map::RationalP1Map) = parent(map.map)
Oscar.polynomial_ring(map::RationalP1Map) = parent(numerator(map.map))
Oscar.coefficient_field(map::RationalP1Map) = base_ring(function_field(map))

Base.show(io::IO, map::RationalP1Map) = print(io, "ℙ¹(z ↦ ",map.map,")")
Base.:(==)(m::RationalP1Map{T},n::RationalP1Map{T}) where T = m.map == n.map
Base.hash(m::RationalP1Map{T},h) where T = hash(m.map,h)

struct MoebiusP1Map{T} <: P1Map{T}
    a::T
    b::T
    c::T
    d::T
    function MoebiusP1Map{T}(a::T,b::T,c::T,d::T) where T
        e = 1-max(frexp_2(a),frexp_2(b),frexp_2(c),frexp_2(d))
        new{T}(ldexp(a,e),ldexp(b,e),ldexp(c,e),ldexp(d,e))
    end
end

Oscar.degree(::MoebiusP1Map) = 1

Base.show(io::IO, map::MoebiusP1Map) = print(io, "ℙ¹(z ↦ (",map.a,"*z + ",map.b,")//(",map.c,"*z + ",map.d,")")
Base.:(==)(m::MoebiusP1Map{T},n::MoebiusP1Map{T}) where T = m.a==n.a && m.b==n.b && m.c==n.c && m.d==n.d
Base.hash(m::MoebiusP1Map{T},h) where T = hash((m.a,m.b,m.c,m.d),h)
          
moebius_map(a::T, b::T, c::T, d::T) where T = MoebiusP1Map{T}(a,b,c,d)

function as_rational_p1_map(m::RationalP1Map{T},n::MoebiusP1Map{T}) where T
    F = function_field(m)
    R = polynomial_ring(m)
    F(R([n.b,n.a]))//F(R([n.d,n.d]))
end

function p1_map(p::Union{AcbPolyRingElem,AbstractAlgebra.Generic.Poly{T}}) where T
    r,z = rational_function_field(coefficient_ring(p),var(parent(p)))
    RationalP1Map{eltype(typeof(p)),typeof(p)}(r(p))
end
function moebius_map(p::Union{AcbPolyRingElem,AbstractAlgebra.Generic.Poly{T}}) where T
    degree(p)==1 || error("only degree-1 maps can be made Moebius")

    c = coefficients(p)
    F = parent(c[0])
    MoebiusP1Map{eltype(typeof(p))}(c[1],c[0],F(0),F(1))
end
            
p1_map(r::AbstractAlgebra.Generic.RationalFunctionFieldElem{T,U}) where {T,U} = RationalP1Map{T,U}(r)
moebius_map(r::AbstractAlgebra.Generic.RationalFunctionFieldElem{T,U}) where {T,U} = moebius_map(p1_map(r))

function p1_map(m::MoebiusP1Map{T}) where T
    r,z = rational_function_field(parent(m.a),:z)
    RationalP1Map{T,typeof(numerator(z))}((m.a*z+m.b)//(m.c*z+m.d))
end

function moebius_map(r::RationalP1Map)
    degree(r)==1 || error("only degree-1 maps can be made Moebius")

    n = coefficients(numerator(r.map))
    d = coefficients(denominator(r.map))
    
    moebius_map(n[1],n[0],d[1],d[0])
end

# Möbius map w/ images of ∞,0,1
function moebius_map(p∞::P1Point{T},p0::Union{P1Point{T},Nothing} = nothing, p1::Union{P1Point{T},Nothing} = nothing) where T
    if p0==p1==nothing
        moebius_map(p∞.num,-conj(p∞.den),p∞.den,conj(p∞.num))
    elseif p1==nothing
        moebius_map(p∞.num,p0.num,p∞.den,p0.den)
    else
        p01 = p1.num*p0.den - p0.num*p1.den
        p1∞ = p∞.num*p1.den - p1.num*p∞.den
        moebius_map(p∞.num*p01,p0.num*p1∞,p∞.den*p01,p0.den*p1∞)
    end
end

# Möbius map w/ sources and images
function moebius_map(m::Pair{P1Point{T},P1Point{T}}...) where T
    length(m)==3 || error("Need three source and three range points")

    inv(moebius_map(m[1].first,m[2].first,m[3].first))∘moebius_map(m[1].second,m[2].second,m[3].second)
end

function moebius_path(p::P1Point{T},q::P1Point{T}) where T # Möbius transformation 0↦p, 1↦q, ∞↦antipode(p)
    moebius_map(antipode(p),p,q)
end

function moebius_map(m::Union{Matrix{T},MatrixElem{T}}) where T
    size(m)==(2,2) || error("Matrix should be 2×2")

    moebius_map(m[1,1],m[1,2],m[2,1],m[2,2])
end

Base.inv(m::MoebiusP1Map{T}) where T = moebius_map(m.d,-m.b,-m.c,m.a)

Base.one(m::MoebiusP1Map{T}) where T = (z = zero(m.a); o = one(m.a); moebius_map(o,z,z,o))
Base.one(::Type{MoebiusP1Map{T}}) where {T <: Complex} = moebius_map(T(1),T(0),T(0),T(1))
Base.zero(r::RationalP1Map{T}) where T = p1_map(zero(r.map))
Base.one(r::RationalP1Map{T}) where T = p1_map(one(r.map))
Oscar.gen(r::RationalP1Map{T}) where T = p1_map(gen(function_field(r)))
Oscar.identity_map(z::P1Point{T}) where T = (F = field(z); moebius_map(F(1),F(0),F(0),F(1)))

function Base.:∘(m::MoebiusP1Map{T},n::MoebiusP1Map{T}) where T
    moebius_map(m.a*n.a + m.b*n.c, m.a*n.b + m.b*n.d, m.c*n.a + m.d*n.c, m.c*n.b + m.d*n.d)
end

function Base.:∘(m::RationalP1Map{T},n::MoebiusP1Map{T}) where T
    rn = as_rational_p1_map(m,n)
    subst(numerator(m),rn)//subst(denominator(m),rn)
end

function Base.:∘(m::MoebiusP1Map{T},n::RationalP1Map{T}) where T
    (m.a*n + m.b) // (m.c*n + m.d)
end

function Base.:∘(m::RationalP1Map{T},n::RationalP1Map{T}) where T
    my_quo(parent(m.map),subst(numerator(m),rn),subst(denominator(m),rn))
end

function Base.:^(m::P1Map{T},n::MoebiusP1Map{T}) where T
    inv(n)∘m∘n
end

is_polynomial(m::MoebiusP1Map) = iszero(m.c)
is_polynomial(m::RationalP1Map) = isone(denominator(m.map))
Base.isone(m::MoebiusP1Map) = m.a == m.d && iszero(m.b) && iszero(m.c)
Base.iszero(m::RationalP1Map) = iszero(m.map)
Base.isone(m::RationalP1Map) = isone(m.map)

function Oscar.derivative(m::RationalP1Map) # seems broken in Oscar
    n, d = numerator(m.map), denominator(m.map)
    if is_constant(d)
        dm = derivative(n) // d
    else
        dm = my_quo(parent(m.map),derivative(n)*d - n*derivative(d),d^2)
    end
    p1_map(dm)
end

function Oscar.integral(m::RationalP1Map)
    if is_polynomial(m)
        p1_map(function_field(m)(integral(numerator(m.map))))
    else
        error("integral of rational map not (yet) implemented")
    end
end

Base.:+(m::RationalP1Map{T}, n::RationalP1Map{T}) where T = p1_map(m.map+n.map)
Base.:-(m::RationalP1Map{T}, n::RationalP1Map{T}) where T = p1_map(m.map-n.map)
Base.:*(m::RationalP1Map{T}, n::RationalP1Map{T}) where T = p1_map(m.map*n.map)
Base.:/(m::RationalP1Map{T}, n::RationalP1Map{T}) where T = p1_map(my_quo(parent(m.map),numerator(m.map)*denominator(n.map),numerator(n.map)*denominator(m.map)))

Base.:+(m::RationalP1Map{T}, n) where T = p1_map(m.map+function_field(m)(n))
Base.:-(m::RationalP1Map{T}, n) where T = p1_map(m.map-function_field(m)(n))
Base.:*(m::RationalP1Map{T}, n) where T = p1_map(m.map*function_field(m)(n))
Base.:/(m::RationalP1Map{T}, n) where T = p1_map(m.map//function_field(m)(n))

Base.:+(m, n::RationalP1Map{T}) where T = p1_map(function_field(n)(m)+n.map)
Base.:-(m, n::RationalP1Map{T}) where T = p1_map(function_field(n)(m)-n.map)
Base.:*(m, n::RationalP1Map{T}) where T = p1_map(function_field(n)(m)*n.map)
Base.:/(m, n::RationalP1Map{T}) where T = p1_map(function_field(n)(m)//n.map)

Base.:^(m::RationalP1Map{T}, n::Integer) where T = p1_map(m.map^n)
Base.inv(m::RationalP1Map{T}) where T = p1_map(inv(m.map))
Base.:-(m::RationalP1Map{T}) where T = p1_map(-m.map)
Base.conj(m::RationalP1Map{T}) where T = p1_map(conj(m.map))

function Oscar.image(m::MoebiusP1Map{T}, p::P1Point{T}) where T
    P1Point{T}(m.a*p.num + m.b*p.den, m.c*p.num + m.d*p.den)
end
(m::P1Map{T})(p::P1Point{T}) where T = image(m,p)

function Oscar.preimage(m::MoebiusP1Map{T}, p::P1Point{T}) where T
    P1Point{T}(m.d*p.num - m.b*p.den, -m.c*p.num + m.a*p.den)
end

#!!! should we do it in reverse if one of num, den is very small?
function homogeneous_poly_eval(p::PolyRingElem{T},num::T,den::T,d = degree(p)) where T
    c = coefficients(p)
    v = c[0]
    inum = one(v)
    for i=1:d
        inum *= num
        v = den*v + inum*c[i]
    end
    v
end

# computing the derivative on-the-spot
function homogeneous_dpoly_eval(p::PolyRingElem{T},num::T,den::T,d = degree(p)) where T
    c = coefficients(p)
    v = c[1]
    inum = one(v)
    for i=2:d
        inum *= num
        v = den*v + inum*i*c[i]
    end
    v
end

function Oscar.image(m::RationalP1Map{T}, p::P1Point{T}) where T
    n = homogeneous_poly_eval(numerator(m.map),p.num,p.den,degree(m))
    d = homogeneous_poly_eval(denominator(m.map),p.num,p.den,degree(m))
    P1Point{T}(n,d)
end

"""returns the gcd g of polynomials a,b, as well as the quotients a/g, b/g.

Allows for approximate calculations when the coefficients of a,b are compatible with a common factor
"""
gcd_quo(a::PolyRingElem{T},b::PolyRingElem{T}) where T = let g = gcd(a,b)
    g,divexact(a,g),divexact(b,g)
end

function gcd_quo(a::PolyRingElem{AcbFieldElem},b::PolyRingElem{AcbFieldElem})
    # we won't use a smart half-gcd algorithm
    va, vb = collect(coefficients(a)), collect(coefficients(b))
    while length(va)>0 && contains(va[end],zero(coefficient_ring(a)))
        pop!(va)
    end
    while length(vb)>0 && contains(vb[end],zero(coefficient_ring(a)))
        pop!(vb)
    end
    while length(va)>0 && length(vb)>0
        if length(va)>length(vb) # make vb longer
            va, vb = vb, va
        end
        c = vb[end] / va[end]
        for i=1:length(va)
            vb[end+1-i] -= c*va[end+1-i]
        end
        pop!(vb)
        while length(vb)>0 && contains(vb[end],zero(coefficient_ring(a)))
            pop!(vb)
        end
    end
    if length(va)>length(vb) # make vb longer
        va, vb = vb, va
    end
    for i=1:length(vb)
        vb[i] /= vb[end]
    end

    vb[end] = one(vb[end])
    g = parent(a)(vb)

    return g, divexact(a,g), divexact(b,g)
end

function my_quo(r::Field,a::PolyRingElem{AcbFieldElem},b::PolyRingElem{AcbFieldElem})
    _, qa, qb = gcd_quo(a,b)

    # major hack! force quotient without executing broken gcd
    q = AbstractAlgebra.Generic.FracFieldElem{AcbPolyRingElem}(qa,qb)
    q.parent = AbstractAlgebra.Generic.FracDict[parent(a)]
    r(q)
end

my_quo(r::Field,a::PolyRingElem{T},b::PolyRingElem{T}) where T = r(a//b)

function yun_algorithm(callback::Function,f)
    g, c, d = gcd_quo(f,derivative(f))
    d -= derivative(c)
    i = 1
    while !is_constant(c)
        g, c, d = gcd_quo(c,d)
        d -= derivative(c)
        callback(g,i)
        i += 1
    end
end

function yun_algorithm(f)
    p = typeof(f)[]
    yun_algorithm(f) do g, _
        push!(p,g)
    end
    p
end

# try to implement https://arxiv.org/abs/2301.07880 ?
function roots_with_multiplicity(callback::Function, p::PolyRingElem{T}) where T
    yun_algorithm(p) do g, i
        degree(g)>0 && for w=roots(g) # bug: roots on 1 fails
            callback(w,i)
        end
    end
end

function roots_with_multiplicity(p::PolyRingElem{T}) where T
    r = Pair{T,Int}[]
    roots_with_multiplicity(p) do w, i
        push!(r,w=>i)
    end
    r
end


"""Return the preimages of p under the map m, with multiplicity, as a list of pairs P1Point=>Int"""
function preimages(m::RationalP1Map{T}, p::P1Point{T}) where T
    r = Pair{P1Point{T},Int}[]
    d = degree(m)
    roots_with_multiplicity(numerator(m.map)*p.den - denominator(m.map)*p.num) do w, i
        push!(r,P1Point{T}(w)=>i)
        d -= i
    end
    d > 0 && push!(r,P1Inf(p)=>d)
    r
end

"""Return the fixed points of the map m, with multiplicity, as a list of pairs P1Point=>Int"""
function fixed_points(m::RationalP1Map{T}) where T
    r = Pair{P1Point{T},Int}[]
    d = degree(m)+1
    roots_with_multiplicity(numerator(m.map) - inflate(denominator(m.map),1,1)) do w, i
        push!(r,P1Point{T}(w)=>i)
        d -= i
    end
    d > 0 && push!(r,P1Inf(coefficient_field(m))=>d)
    r    
end

function periodic_points(m::RationalP1Map, p::Int)
    error("!!! periodic points not yet implemented")
end

"""Return the critical points of the map m, with multiplicity, as a list of pairs P1Point=>Int"""
function critical_points(m::RationalP1Map{T}) where T
    num, den = numerator(m.map), denominator(m.map)
    r = Pair{P1Point{T},Int}[]
    d = 2degree(m)-2
    roots_with_multiplicity(derivative(num)*den - num*derivative(den)) do w, i
        push!(r,P1Point{T}(w)=>i)
        d -= i
    end
    d > 0 && push!(r,P1Inf(coefficient_field(m))=>d)
    r
end

function image_df(m::RationalP1Map{T}, p::P1Point{T}) where T
    num, den = numerator(m.map), denominator(m.map)
    d = degree(m)
    
    qnum = homogeneous_poly_eval(num,p.num,p.den,d)
    qden = homogeneous_poly_eval(den,p.num,p.den,d)

    #!!! make more robust, with homogeneous coordinates
    if maybeinf(p) # p = ∞
        cnum, cden = coefficients(num), coefficients(den)
        qdiff = cnum[d]*cden[d-1] - cnum[d-1]*cden[d]
    else
        qdiff = homogeneous_dpoly_eval(num,p.num,p.den,d)*qden - homogeneous_dpoly_eval(den,p.num,p.den,d)*qnum
        qdiff /= p.den
    end
    qdiff *= (abs2(p.den) + abs2(p.num)) / (abs2(qnum) + abs2(qden))

    P1Point{T}(qnum,qden), qdiff
end

Oscar.derivative(m::RationalP1Map{T}, p::P1Point{T}) where T = image_df(m,p)[2]

function p1_map(poles::Vector, zeros::Vector, img::Pair{P1Point{T},P1Point{T}}) where T
    r, z = rational_function_field(field(img.first),:z)
    zp = numerator(z)
    
    function make_poly(v)
        all(a->isa(a,P1Point{T}) || isa(a,Pair{P1Point{T},U} where U <: Integer),v) || error("zeros/poles should be T or T=>Integer")

        p = one(zp)
        for a=v
            if isa(a,P1Point{T})
                p *= a.den*zp - a.num
            else
                p *= (a.first.den*zp - a.first.num)^a.second
            end
        end
        p
    end
    
    m = p1_map(my_quo(r,make_poly(zeros),make_poly(poles)))
    m * (zcoord(img.second) / zcoord(m(img.first)))
end

"""compute the (t,u) in [0,1]x[0,1] such that γ(t) = f(δ(u)).
returns a list of (t,u,Im(γ^-1*f*δ)'(u),γ(t),δ(u))
γ, δ are Möbius transformations, and f is a rational map.
"""
function p1_intersect(γ::MoebiusP1Map{T},f::RationalP1Map{T},δ::MoebiusP1Map{T}) where T
    U = real(T)
    f₁ = inv(γ)∘f∘δ
    p = imag(numerator(f₁)*conj(denominator(f₁)))
    
    intersections = @NamedTuple{t::U,u::U,direction::Int,γt::P1Point{T},δu::P1Point{T}}[]
    for u=roots(p,isolate_real=true)
        is_real(r) || continue
        pu = P1Point(u)
        z = f₁(pu)
        isinf(z) && continue
        t = real(coord(z))
        0 ∈ imag(coord(z)) || continue
        0 ≤ t ≤ 1 || continue #!!! interval arithmetic test
        slope = homogeneous_dpoly_eval(poly,u.num,u.den,2degree(f))
        slopei = imag(slope)
        direction = (0 ∈ slopei ? 0 : (slopei > 0 ? 1 : -1))
        push!(intersections, (t = t, u, direction, γt = γ(z), δu = δ(pu)))
    end
    intersections
end    

"""find a Möbius transformation that sends the last of points to P1Inf, and either
- matches points and extra as well as possible, if oldpoints is a list;
- is only a rotation, otherwise."""
function rotating_moebius_map(points::Vector{P1Point{T}}, oldpoints = nothing) where T
    m = inv(moebius_map(points[end])) # sends points[end] to ∞

    if extra≠nothing
        @assert length(points) == length(oldpoints)
        @assert oldpoints[end] ≎ P1Inf(oldpoints[end])
        points = m.(points)
        n = length(points)

        theta = zero(field(points[end]))
        norm = theta
        for i=1:n
            proj = xycoord(points[i])
            oldproj = xycoord(oldpoints[i])
            theta += conj(proj)*oldproj
            norm += abs2(proj)
        end
        if 0∈norm || abs2(theta) < 0.7norm
            theta = zero(theta)
        end

        if theta == 0
            # as last resort, just force the point of largest distance to be on the positive real axis
            theta = one(theta) # default
            norm = zero(theta)
            for i=1:n
                proj = xycoord(points[i])
                newnorm = abs2(proj)
                if newnorm > norm
                    norm = newnorm
                    theta = conj(proj)
                end
            end
        end
        theta /= abs(theta) # make it of norm 1
        
        m = moebius_map(theta,zero(theta),zero(theta),one(theta))∘m
    end
    m
end

#= given a set of points on S^2 \subset R^3, there is, up to rotations
   of the sphere, a unique M\"obius transformation that centers these
   points, i.e. such that their barycentre is (0,0,0). This follows
   from GIT, as Burt Totaro told me:

   Dear Laurent,

   Geometric invariant theory (GIT) gives a complete answer to your
   question. More concretely, the answer follows from the Kempf-Ness
   theorem in GIT, as I think Frances Kirwan first observed.

   Namely, given a sequence of N points p_1,...,p_N on the 2-sphere,
   there is a Mobius transformation that moves these points
   to have center of mass at the origin of R^3 if and only if either
   (1) fewer than N/2 of the points are equal to any given point
   in the 2-sphere; or
   (2) N is even, N/2 of the points are equal to one point
   in the 2-sphere, and the other N/2 points are equal
   to a different point in the 2-sphere.
   (In GIT terminology, condition (1) describes which
   N-tuples of points in S^2 = CP^1 are "stable",
   and (2) describes which N-tuples are "polystable"
   but not stable.) In particular, if p_1,...,p_N are all distinct
   and N is at least 2, then they can be centered
   by some Mobius transformation.

   This result was the beginning of many developments in GIT,
   such as Donaldson's notion of "balanced" metrics.
   Here is a good survey (where Theorem 4.13 is the statement
   above).

   R. P. Thomas. Notes on GIT and symplectic reduction for bundles
   and varieties. arXiv:math/0512411

   Burt Totaro

   This is also proven in <Cite Ref="MR2121737">.
=#

"""x is a "shifting" parameter; it is a vector in R^3, and
describes a Möbius transformation with north-south dynamics.
More precisely, let t=|x|. in R^3, the transformation sends
     P to (2(1-t)P+(2-t+(v*P))v)/(1+(1-t)^2+(2-t)(v*P)).
In particular, for t=0 it sends everything to v, and for t=1 it fixes P.
"""
function __solve_barycenter(x, p::Vector{NTuple{3,T}}) where T
    t = sqrt(sum(x.*x))
    n = length(p)
    
    sum = zero(x)
    for i=1:n
        z = sum(p[i].*x)
        d = 1 + (1-t)*(1-t) + (2-t)*z

        for j=1:3
            sum[j] += (2(1-t)*p[i][j] + (2-t+z)*x[j]) / d
        end
    end
    sum ./ n
end


function __find_barycenter(r3points::Vector{NTuple{3,T}}) where T
    x0 = zeros(T,3)
    problem = NonlinearSolve.NonlinearProblem(__solve_barycenter,x0,r3points)
    sol = NonlinearSolve.NonlinearSolve(problem)

    @info sol

    sol.u
end

"""Compute a Möbius transformation that sends points[end] to ∞,
the barycenter to 0, and makes the new points as close as possible
to oldpoints by a rotation fixing 0-∞.
oldpoints is allowed to be 'nothing', in which case we just return a good Möbius transformation.
"""
function normalizing_moebius_map(points, oldpoints = nothing)
    n = length(points)
    F = field(points[1])
    
    if n==2
        return inv(moebius_map(points[1],points[2]))
    elseif n==3
        r3 = 1/sqrt(F(3))
        return moebius_map(points[1]=>P1Point(r3),points[2]=>P1Point(-r3),points[3]=>P1Inf(points[3]))
    end

    barycenter = __find_barycenter(R3_P1Point.(points))
    dilate = sqrt(barycenter.*barycenter)
    if iszero(dilate)
        map = identity_map(points[1])
    else
        map = inv(moebius_map(P1Point_R3(-barycenter/dilate)))*(1-dilate)
    end
    
    if oldpoints==nothing
        rotating_moebius_map(map(points[end]))∘map
    else
        newpoints = map.(points)
        rotating_mobius_map(newpoints,oldpoints)∘map
    end
end

