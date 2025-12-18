function import_ring(i::Int) 
  js = JSON.parsefile( joinpath(@__DIR__, "data","FusionRingsJSON", "ring"*string(i)*".json") )
  fc = [ js["formal_code"][i] for i in 1:4 ]
  r = fc[1]
  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r 
      mt[i,j,k] = Int.(js["mult_tab"][i][j][k])
  end
  FusionRings.fusion_ring( mt, formal_code = fc)
end

export fusion_ring_list

fusion_ring_list =  [ import_ring(i) for i in 1:28451 ]

export frl

frl = fusion_ring_list