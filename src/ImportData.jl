# Just some test code 

function frl(i::Int) 
  js = JSON.parsefile("/home/gert/Tests/JSONExport/ring_"*string(i)*".json")
  fc = [ js["formal_code"][i] for i in 1:4 ]
  r = fc[1]
  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r 
      mt[i,j,k] = Int.(js["mt"][i][j][k])
  end
  FusionRings.fusion_ring( mt, formal_code = fc)
end
