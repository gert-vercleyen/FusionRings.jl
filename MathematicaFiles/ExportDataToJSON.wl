(* ::Package:: *)

(* ::Input:: *)
(**)


PacletDirectoryLoad["~/Projects"];
<<Anyonica`


ringtojson[ r_FusionRing ]:= 
	Module[
		{ fc, asso, cats, catProps, categorifiableQ, sl2zreps, tpdecomps, subfusrings },
		fc = FC @ r;
		
		cats = FusionCategories @ r;
		
		categorifiableQ = 
			Which[
                MissingQ[cats] && ( Rank[r] > 7 || Mult[r] > 1 ),
                    Null,
				MissingQ[cats] && Rank[r] < 8 && Mult[r] == 1 && FC[r] =!= {7,1,2,11},
                    False,
                FC[r] === {7,1,2,11},
                    Null,
				True,
					True
			];
		
		If[
			TrueQ[categorifiableQ],
			catProps = MemberQ[True] /@ 
				Comap[ Map /@ { BraidedQ, UnitaryQ, SphericalQ, RibbonQ, ModularQ }, cats ],
			catProps = ConstantArray[ Null, 5 ]
		];

        tpdecomps = WhichDecompositions @ r;

        subfusrings = SubFusionRings @ r;

		sl2zreps = 
			With[{md = ModularData[r]},
				If[ 
					MissingQ[md],
					Null,
					Table[ 
						"reps_"<>ToString[i] -> Evaluate[Chop /@ ReIm @ N[  md[[i]], 64 ]],
						{i,Length @ md}
					]
				]
			];
			
		asso =<|
			"references" -> {""}
			,
			"software" -> {"https://doi.org/10.5281/zenodo.10686859"}
			,
			"mult_tab" -> MT @ r
			,
			"texnames" -> TeXNames @ r
			,
			"barcode"  -> ToString @ Barcode @ r
			,
			"formal_code" -> FC @ r
			,
			"tensor_product_decompositions" ->
                If[
                    tpdecomps === {},
                    {},
                    <|
                        "identifiers" -> {"anyonwiki_code"},
                        "value" -> Map[ FC, tpdecomps, {2} ]
                    |>
                ]
			,
			"non_trivial_sub_fusion_rings" ->
                If[
                    subfusrings === {},
                    {},
                    <|
                        "identifiers" -> {{"injection","anyonwiki_code"}},
                        "value" -> MapAt[FC, subfusrings, {All,2}]
                    |>
                ]
			,
			"numeric_characters" ->
				Map[Chop@* ReIm, N[ FusionRingCharacters @ r, 64 ], {2} ]
			,
			"numeric_frobenius_perron_dimensions" -> N[ ReIm @ FPDims @ r, 64 ]
			,
			"numeric_frobenius_perron_dimension" -> N[ ReIm @ FPDim @ r, 64 ]
			,
			"numeric_projective_SL2Z_reps" -> sl2zreps
			,
			"has_categories_with_props" ->
                If[
                    TrueQ @ categorifiableQ,
                    Transpose @ {
                        { "Braided", "Unitary", "Spherical", "Ribbon", "Modular"  },
                        catProps
                    },
                    Null
                ]
			,
			"categorifiable" -> categorifiableQ
			,
			"categorifications" ->
				If[ 
					TrueQ @ categorifiableQ, 
					<|
						"identifiers" -> {"anyonwiki_code"},
						"categories" -> FC /@ cats
					|>,
					Null
			  ]			
			,
			"comments" -> {""}
            ,
            "info" ->
                "Fusion ring. mult_tab: structure constants." <>
                " barcode & formal_code: unique identifiers see (DOI: 10.1063/5.0148848)." <>
                " non_trivial_sub_fusion_rings: tuples where the first element = elements of ring that" <>
                " form subring isomorphic to subring identified by second element of the tuple." <>
                " software: doi of original software used to represent fusion ring." <>
                " references: doi of paper from which data was obtained." <>
                " categorifiable: false=not categorifiable, null= unknown." <>
                " categorifications: if categorifiable then anyonwiki codes of" <>
                " pivotal (braided) fusion cats that categorify ring." <>
                " numeric_projective_SL2Z_reps: each rep consists of a generalized S-matrix and" <>
                " a vector of vectors representing the ln(diag(T))/(2 pi i) of a generalized T-matrix."
		|>

	];
	
exportringstojson[dir_String] :=
    With[{ l = Length @ FRL },
        Do[
            If[ Mod[i,1000] == 0, Print["First "<> ToString[i] <> " rings exported. "] ];
            exportringtojson[dir] @ FRL[[i]],
            { i, l }

        ]
    ]

(* NOTE: only works for rings of which the formal code is known *)
exportringtojson[dir_String][r_FusionRing] :=
	Module[
		{asso, codeToString, info, codestring},
		codestring =
			StringReplace[ ToString @ FC @ r, { " "->"", ","->"_","{"->"","}"->"" } ];

		Export[
			FileNameJoin[{dir,"fr_" <> codestring <>".json"}],
			ringtojson @ r,
			"JSON"
		];
	]
	
TeXNames[r_FusionRing] :=
	Module[{names},
		If[
			Length[ names = Names @ r ] == 0,
			{ FRToTexString[FC @r] },
			If[ 
				SongNameQ[ # ], 
				If[ 
					MemberQ[ FC @ r ] @ probSongCodes, 
					fixProbSongName @ FC @ r, 
					FixSongName @ # 
				], 
				FixLatex @ ToString[ TeXForm @ # ] 	
			]& /@ names
		]	
	];


SongNameQ[ s_String ] := StringMatchQ[ s, "\!\(\*TemplateBox[List["~~x__];

FixSongName[ songName_String ] :=
	If[
		TrueQ @ SongNameQ[songName],
		songName //
		StringReplace[
		"\!\(\*TemplateBox[List["~~Shortest[x__]~~"\"\","~~Shortest[y__]~~"\"\","~~z__~~"], \"Subsuperscript\", Rule[SyntaxForm, SubsuperscriptBox]]\)":>x<>"@"<>y<>"@"<>z]//
		StringReplace[{"\\\\("->"(","\\\\)"->")","\\\\!"->""}] //
		StringReplace["(\\\\*SubscriptBox[" ~~ Shortest[x__] ~~ ", (" ~~Shortest[ y__] ~~ ")])" :>
		    x <> "_{" <> y <> "}"]//
		StringReplace["*StyleBox["~~x__~~y_?DigitQ~~Shortest[z__]~~"False]]" :> "\\mathbf{"<>y<>"}"]//
		StringReplace[{"[DoubleStruckCapitalZ]"->"\\mathbb{Z}", "[LeftTriangleEqual]"->"\\trianglelefteq","[Cross]"->"\\times"}]//
		StringReplace[{"\\"->"","\""->""}]//
		StringReplace[{"("->"",")"->""}]//
		StringReplace[{"["->"\\left[","]"->"\\right]"}]//
		StringReplace[{"mathbb"->"\\mathbb","tria"->"\\tria","times"->"\\times","mathbf"->"\\mathbf","Id"->"\\text{Id}"}]//
		StringReplace[x__~~"@"~~y__~~"@"~~z__:>x<>"_{"<>y<>"}^{"<>z<>"}"]//
		StringReplace["|"->"\\|"]//
		StringReplace["\\left[CircleTimes\\right]"->"\\otimes"],
		songName
	];

FixLatex[str_String] :=
	str //
	StringReplace["\\unicode{22ca}" -> "\\rtimes "] //
	StringReplace[Shortest["\\text{$\\times $"~~x__~~"}"] :> "\\times \\text{"<>x<>"}" ] //
	StringReplace[Shortest["\\left.\\text{"~~x__~~"}"~~y__~~"\\right)"] :>"\\text{"<>x<>"} "<>y<>")" ] //
	StringReplace[Shortest["\\text{"~~x__~~"$\\times $"~~y__~~"}"] :> "\\text{"<>x<>"} \\otimes \\text{" <> y <> "}"] //
	StringReplace[Shortest["\\text{"~~P___~~"SU(2}"] :> "\\text{"<>P<>"SU}(2" ]//
	StringReplace[ "\\text{}" ->""] //
	StringReplace["Rep(}" -> "Rep}("] //
	StringReplace["TY(}" -> "TY}("] //
	StringReplace["HI(}" -> "HI}("] //
	StringReplace["Adj(}" -> "Adj}("]//
	StringReplace["\\text{Fib$\\otimes $Fib}"->"\\text{Fib}\\otimes\\text{Fib}"]//
	StringReplace["\\text{Potts}"->"\\text{TY}(\\mathbb{Z}_3)"]//
	StringReplace["\\mathbb{Z}_2\\text{$\\otimes $Ising}"->"\\mathbb{Z}_2\\otimes \\text{Ising}"]//
	StringReplace["\\left.\\mathbb{Z}_2\\text{$\\otimes $Rep}(D_3\\right)"->"\\mathbb{Z}_2\\otimes \\text{Rep}(D_3)"] // 
	StringReplace["\\text{Fib$\\otimes $Ising}"->"\\text{Fib}\\otimes\\text{Ising}"]//
	StringReplace["\\times \\text{Rep}"->"\\otimes \\text{Rep}"]//
	StringReplace["\\mathbb{Z}_2\\times \\mathbb{Z}_2]"->"\\mathbb{Z}_2\\otimes \\mathbb{Z}_2]"]//
	StringReplace["\\mathbb{Z}_2\\times \\text{Ising}"->"\\mathbb{Z}_2\\otimes \\text{Ising}"]//
	StringReplace["\\trianglelefteq \\mathbb{Z}_{2}\\otimes\\mathbb{Z}_{2}"->"\\trianglelefteq \\mathbb{Z}_{2}\\times\\mathbb{Z}_{2}"]
	
(* Some names need to be fixed manually *)
probSongCodes = {{2,1,0,1},{2,1,0,2},{4,1,0,1},{4,1,0,4},{4,1,2,1},{4,1,2,4},{6,1,4,1},{6,1,4,8},{8,1,0,1},{8,1,0,23},{8,1,4,1},{8,1,4,16},{8,1,4,17},{8,1,6,2},{8,1,6,10},{9,1,4,3},{9,1,4,11}};
fixProbSongName = 
	Association @ 
	Thread[ 
		{{2,1,0,1},{2,1,0,2},{4,1,0,1},{4,1,0,4},{4,1,2,1},{4,1,2,4},{6,1,4,1},{6,1,4,8},{8,1,0,1},{8,1,0,23},{8,1,4,1},{8,1,4,16},{8,1,4,17},{8,1,6,2},{8,1,6,10},{9,1,4,3},{9,1,4,11}}
		->
		{"\\left[I\\trianglelefteq I\\right]_{\\mathbf{1}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq I\\right]_{\\mathbf{1}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\right]_{\\mathbf{1}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\right]_{\\mathbf{1}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\right]_{\\mathbf{2}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\right]_{\\mathbf{2}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_3\\right]_{\\mathbf{1}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_3\\right]_{\\mathbf{1}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\times\\mathbb{Z}_2\\right]_{\\mathbf{1}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\times\\mathbb{Z}_2\\right]_{\\mathbf{1}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_4\\right]_{\\mathbf{1}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_4\\right]_{\\mathbf{1}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_2\\times\\mathbb{Z}_2\\right]_{\\mathbf{2}\\|1}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_4\\right]_{\\mathbf{3}\\|0}^{\\text{Id}}","\\left[I\\trianglelefteq \\mathbb{Z}_4\\right]_{\\mathbf{3}\\|1}^{\\text{Id}}","\\left[\\mathbb{Z}_2\\trianglelefteq \\mathbb{Z}_6\\right]_{\\mathbf{1}\\|1}^{(\\mathbf{2}\\ \\mathbf{3})}","\\left[\\mathbb{Z}_2\\trianglelefteq \\mathbb{Z}_6\\right]_{\\mathbf{1}\\|1}^{(\\mathbf{2}\\ \\mathbf{3})}"}	
	];
	
FRToTexString[ {a_,b_,c_,d_} ]:=
	"\\text{FR}^{" <> ToString[a] <> "," <> ToString[b] <> "," <> ToString[c] <> "}_{" <> ToString[d] <> "}";
	


["~/Projects/FusionRings.jl/src/data/FusionRingsJSON/"]


(* ::Input:: *)
(*ringtojson[FRL[[505]]]*)



