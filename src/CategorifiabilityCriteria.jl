# ========================================================================
#   COMMUTATIVE SCHUR PRODUCT CRITERION
# ========================================================================

# TODO: check whether s is correct

# csp_criterion( fusion_ring ) returns true if fusion_ring does not have a unitary categorification due to the commutative schur product criterion

# Compat helpers (ported from Mathematica naming)
ISREAL(x) = isreal(x)
lesserthan(x, y) = x < y

function csp_criterion( ring::FusionRing )
    is_commutative( ring ) || error("Commutative Schur Product Criterion only applies to commutative fusion rings")

    chars = characters( ring )
    r = rank( ring )

    for j1 in 1:r, j2 in 1:r, j3 in 1:r
        s = sum( chars[ i, j1 ] * chars[ i, j2 ] * chars[ i, j3 ] / chars[ i, 1 ] for i in 1:r )
        if ISREAL( s ) && lesserthan( s, 0 )
            return true
        else 
            continue
        end
    end

    return false
end


# ========================================================================
#   PIVOTAL DRINFELD CENTER CRITERION
# ========================================================================

# pdc_criterion returns true if ring has no complex pivotal categorification due to the pivotal Drinfeld center criterion
function pdc_criterion( r::FusionRing )
#
# PDCCriterion[ ring_FusionRing?CommutativeQ ] :=
#   Module[{chars,c},
#     chars =
#       FusionRingCharacters[ring];
#     c =
#       #.ConjugateTranspose[#]& /@ chars;
#
#     Catch[
#       Do[
#         If[
#           And @@ Flatten @ AlgebraicIntegerQ[ c[[j]] / c ],
#           Throw[ False ]
#         ],
#         { j, Length[c] }
#       ];
#       True
#     ]
#
#   ];
end



#  ========================================================================
#    PSEUDO-UNITARY DRINFELD CENTER CRITERION
#  ========================================================================
#  Returns True if ring has no complex pseudo-unitary categorification
#
#
#PackageExport["PUDCC"]
#
#PUDCC[ ring_FusionRing?CommutativeQ ] :=
#  With[{ chars = FusionRingCharacters[ring] },
#    And @@ Flatten @
#    AlgebraicIntegerQ[ chars[[1]].ConjugateTranspose[chars[[1]]] / chars ]
#  ];
#
#
#(*========================================================================
#    D-NUMBER CRITERION
#  ========================================================================
#  The function returns true if the fusion ring cannot be categorified.
#*)
#
#PackageExport["DNCriterion"]
#
#DNCriterion[ ring_FusionRing?CommutativeQ ] :=
#  Module[{ chars, c, DNumberQ, n },
#    chars = FusionRingCharacters[ring];
#    c = RootReduce[ #.ConjugateTranspose[#]& /@ chars ];
#
#    DNumberQ[ x_ ] :=
#      Module[{ p, a, y },
#        If[ !AlgebraicIntegerQ[x], Return[ True ] ];
#
#        p =
#          MinimalPolynomial[x][y];
#
#        a =
#          ( Rest @ MonomialList[p] /. y -> 1 );
#
#        n =
#          Exponent[ p, y ];
#
#        And @@
#        Table[
#          Mod[ a[[i]]^n, a[[-1]]^i ] == 0,
#          { i, Length[a] }
#        ]
#      ];
#
#    Not[ And @@ ( DNumberQ /@ c ) ]
#
#  ];
#
#(*========================================================================
#    EXTENDED CYCLOTOMIC CRITERION
#  ========================================================================
#  The function returns True if the fusion ring cannot be categorified.
#
#
#
#
#*)
#
#
#
#
#(*========================================================================
#    LAGRANGE CRITERION
#  ========================================================================
#  The function returns False if the fusion ring cannot be categorified.
#
#*)
#
#PackageExport["LCriterion"]
#
#LCriterion[ ring_FusionRing ] :=
#  Module[{ subringDims },
#    subringDims =
#      FPDim /@
#      DeleteDuplicates @
#      SubFusionRings[ ring ][[;;,2]];
#
#    Not[And @@ AlgebraicIntegerQ @ RootReduce[ FPDim[ring]/subringDims ]]
#
#  ];
#
#
#(*========================================================================
#    ZERO SPECTRUM CRITERION
#  ========================================================================
#  The function returns True if the fusion ring cannot be categorified.
#*)
#
#PackageExport["ZSCriterion"]
#
#ZSCriterion::usage =
#  "ZSCriterion[ fusionRing ] returns True if the fusion ring fusionRing cannot be categorified due to " <>
#  "the Zero Spectrum criterion. (arXiv:2203.06522v1)";
#
#(*See arXiv:2203.06522v1 for more info.";*)
#
#SetAttributes[ ZSC, Listable ];
#
#ZSCriterion[ ring_FusionRing ] :=
#  Module[{ mt, non0Cons, ones, d, i1, i2, i3, i4, i5, i6, i7, i8, i9, matches1, matches2, matches3, matches4 },
#    mt =
#      MT[ring];
#    non0Cons =
#      NonZeroStructureConstants[ring];
#    ones =
#      Position[ mt, x_Integer/; x == 1 ];
#    d =
#      CC[ring] /@ Range[Rank[ring]];
#
#    Catch[
#      Do[
#        { i2, i1, i3 } = ind;
#        matches1 = Cases[ non0Cons, { i4_, i1, i6_ } ];
#        Do[
#          { i4, i6 } = ind1[[{1,3}]];
#          matches2 = Cases[ non0Cons, { i5_, i4, i2 } ];
#          Do[
#            i5 = ind2[[1]];
#            If[
#              mt[[ i5, i6, i3 ]] != 0 &&
#              MemberQ[ crit1[ {i1,i2,i3,i4,i5,i6}, d, mt ], 1 ],
#              matches3 = Cases[ non0Cons, { i7_, i9_, i1 } ];
#              Do[
#                { i7, i9 } = ind3[[{1,2}]];
#                matches4 = Cases[ non0Cons, { i2, i7, i8_ } ];
#                Do[
#                  i8 = ind4[[3]];
#                  If[
#                    mt[[ i8, i9, i3 ]] != 0 &&
#                    crit3[ { i4, i5, i6, i7, i8, i9 }, d, mt ] == 0 &&
#                    MemberQ[ crit2[ { i1, i2, i3, i7, i8, i9 }, d, mt ], 1 ],
#                    Throw[ True ]
#                  ]
#                  ,{ ind4, matches4 }]
#                ,{ ind3, matches3 }]
#            ]
#            ,{ ind2, matches2 }]
#          ,{ ind1, matches1 }]
#        , { ind, ones } ];
#      False
#    ]
#
#  ];
#
function crit1( i::Vector{Int64}, d::Vector{Int64}, mt::Array{Int64,3} )::Bool
    r = size( mt, 1 )
    sum( mt[ i[5], i[4], k ] * mt[ i[3], d[ i[1] ], k ] for  k in 1:r ) == 1 ||
    sum( mt[ i[2], d[ i[4] ], k ] * mt[ i[3], d[ i[6] ], k ] for k in 1:r ) == 1 ||
    sum( mt[ d[ i[5] ], i[2], k ] * mt[ i[6], d[ i[1] ], k ] for k in 1:r ) == 1
end

function crit2( i::Vector{Int64}, d::Vector{Int64}, mt::Array{Int64,3} )::Bool
    r = size( mt, 1 )
    sum( mt[ i[2], i[4], k ] * mt[ i[3], d[ i[6] ], k ] for  k in 1:r ) == 1 ||
    sum( mt[ i[5], d[ i[4] ], k ] * mt[ i[3], d[ i[1] ], k ] for k in 1:r ) == 1 ||
    sum( mt[ d[ i[2] ], i[5], k ] * mt[ i[1], d[ i[6] ], k ] for k in 1:r ) == 1
end

function crit3( i::Vector{Int64}, d::Vector{Int64}, mt::Array{Int64,3} )::Bool
    r = size( mt, 1 )
    sum(
        mt[ i[1], i[4], k ] *
        mt[ d[ i[2] ], i[5], k ] *
        mt[ i[3], d[ i[6] ], k ]
        for k in 1:r
    ) == 0
end

#PackageExport["OSCriterion"]
#
#OSCriterion::usage =
#  "OSCriterion[ fusionRing ] returns True if the fusion ring fusionRing cannot be categorified due to " <>
#  "the One Spectrum criterion. (arXiv:2203.06522v1)";
#
#OSCriterion[ ring_FusionRing ] :=
#  Module[{ mt, non0Cons, zeros, d, i1, i2, i3, i4, i5, i6, i7, i8, i9, matches1, matches2, matches3, matches4 },
#    mt = MT[ring];
#    non0Cons = NZSC[ring];
#    zeros = Position[ mt, x_Integer/; x === 0 ];
#    d = CC[ring] /@ Range[Rank[ring]];
#
#    Catch[
#      Do[
#        { i2, i1, i3 } = ind;
#        matches1 = Cases[ non0Cons, { i4_, i1, i6_ } ];
#        Do[
#          { i4, i6 } = ind1[[ { 1, 3 } ]];
#          matches2 = Cases[ non0Cons, { i5_, i4, i2 } ];
#          Do[
#            i5 = ind2[[ 1 ]];
#            If[
#              mt[[ i5, i6, i3 ]] =!= 0,
#              matches3 = Cases[ non0Cons, { i7_, i9_, i1 } ];
#              Do[
#                { i7, i9 } = ind3[[ { 1, 2 } ]];
#                matches4 =
#                  Cases[
#                    Range @ Rank @ ring, i0_ /;
#                    mt[[i4,i7,i0]] === 1 &&
#                    mt[[i6,d[[i9]],i0]] == 1 &&
#                    MemberQ[ crit1[ { i9, i0, i6, i7, i4, i1  }, d, mt ], 1 ]
#                  ];
#                Do[
#                  i0 = ind4;
#                  If[
#                    !MissingQ[
#                      FirstCase[
#                        Range @ Rank @ ring, i8_ /;
#                        mt[[i2,i7,i8]] =!= 0 && mt[[i8,i9,i3]] =!= 0 &&
#                        mt[[d[[i5]],i8,i0]] === 1 &&
#                        MemberQ[ crit1[ { i7, i2, i8, i4, i5, i0  }, d, mt ], 1 ] &&
#                        MemberQ[ crit1[ { i9, i8, i3, i0, i5, i6  }, d, mt ], 1 ] &&
#                        crit3[ { i4, i5, i6, i7, i8, i9 }, d, mt ] == 1
#                      ]
#                    ],
#                    Throw @ True
#                  ]
#                ,{ ind4, matches4 }]
#              ,{ ind3, matches3 }]
#            ]
#            ,{ ind2, matches2 }]
#          ,{ ind1, matches1 }]
#        ,{ ind, zeros }
#      ];
#      False
#    ]
#
#  ];
