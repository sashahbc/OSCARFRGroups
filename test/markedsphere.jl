@testset "marked spheres" begin
  field = AcbField(128)
  n = 6
  pts = equidistributed_p1_points(n, field)

  s = new_marked_sphere(pts)
  @test check_triangulation(triangulation(s)) === true
  @test length(vertices_of_marked_sphere(s)) == n
  for p in vertices_of_marked_sphere(s)
    @test any(q -> proj_equiv(p, q), pts)
  end

  g = Oscar.group(s)
  @test Int(GAP.Globals.RankOfFreeGroup(g)) == n - 1
  @test length(tree_edges(s)) == n - 1

  mk = marking(s)
  @test mk !== GAP.Globals.fail
  modelfam = GAP.Globals.FamilyObj(one(GapObj(sphere_model(s))))
  for x in GAP.Globals.GeneratorsOfGroup(g)
    img = GAP.Globals.Image(mk, x)
    @test GAP.Globals.FamilyObj(img) == modelfam
  end

  # explicit model with a non-default ordering
  model2 = sphere_group([1, 3, 5, 2, 4, 6])
  s2 = new_marked_sphere(pts, model2)
  @test sphere_model(s2) === model2
  @test check_triangulation(triangulation(s2)) === true

  # wiggled_marked_sphere / shallow_copy. A gentle perturbation, not one
  # sending a data point to infinity: that extreme a wiggle can blow up a
  # face past a hemisphere and break the "barycentre is inside its own
  # face" assumption behind check_triangulation's Locate check -- a
  # genuine geometric effect, not a bug, but not what this test is for.
  m = moebius_map(field(1), field(1//10), field(-1//10), field(1))
  sw = wiggled_marked_sphere(s, m)
  @test check_triangulation(triangulation(sw)) === true
  @test length(vertices(triangulation(sw))) == length(vertices(triangulation(s))) &&
    length(edges(triangulation(sw))) == length(edges(triangulation(s))) &&
    length(faces(triangulation(sw))) == length(faces(triangulation(s)))
  @test tree_cost(sw) == tree_cost(s)
  @test Oscar.group(sw) === Oscar.group(s)
  @test marking(sw) === marking(s)

  sc = shallow_copy(s)
  @test check_triangulation(triangulation(sc)) === true
  @test tree_cost(sc) == tree_cost(s)
  @test marking(sc) === marking(s)
end
