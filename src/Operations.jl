function permute(r::FusionRing,perm::Array{Int,1})::FusionRing
end

function which_permutation(r1::FusionRing,r2::FusionRing)::Array{Int,1}
  
end

function sort(r::FusionRing;sortby="fpdims")::FusionRing

end

function perm_vec_qd(r::FusionRing)::Array{Int,1}

end

function perm_vec_sd_conj(r::FusionRing)::Array{Int,1}

end

function tensor_product(r1::FusionRing,r2::FusionRing)::FusionRing

end

⊗(r1::FusionRing,r2::FusionRing) = tensor_product(r1,r2)

function replace_by_known(r::FusionRing)::FusionRing
end

