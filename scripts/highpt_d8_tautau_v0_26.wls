(* v0.26 exact HighPT Run-II di-tau recast for active d=8 UHE operators.
   WC labels:
     charm   WC["l2u2D22",{3,3,2,2}]
     strange WC["l2d2D22",{3,3,2,2}]
   Positivity: C <= 0.
   EFTorder->4 keeps SM x d8 interference, the consistent O(Lambda^-4) truncation.
*)

loadHighPT[] := Module[{ok},
  ok = Quiet[Check[Needs["HighPT`"]; True, False]];
  If[!ok,
    Print["HighPT not found; installing official package..."];
    Quiet[Check[
      Import["https://github.com/HighPT/HighPT/raw/master/install.m"];
      Needs["HighPT`"]; ok = True,
      ok = False
    ]];
  ];
  If[!ok, Print["FATAL: HighPT could not be loaded."]; Exit[2]];
];
loadHighPT[];

$PrintingProcessInfo = False;
$ParallelHighPT = True;
outdir = FileNameJoin[{Directory[], "results", "v0.26"}];
If[!DirectoryQ[outdir], CreateDirectory[outdir, CreateIntermediateDirectories -> True]];

scanOne[label_, lam_] := Module[{wc, chiExpr, f, grid, chiMin, allowed, lower, best, csv},
  InitializeModel["SMEFT", EFTorder -> 4, OperatorDimension -> 8, EFTscale -> lam];
  wc = WC[label, {3,3,2,2}];
  chiExpr = Total@ChiSquareLHC["di-tau-ATLAS",
      Coefficients -> {wc}, EFTorder -> 4, OperatorDimension -> 8,
      EFTscale -> lam, RescaleError -> True];
  f[c_?NumericQ] := N[chiExpr /. {wc -> c, Conjugate[wc] -> c}];
  (* Dense scan around the strict unitarity interval plus a wide tail for the collider bound. *)
  grid = Join[
    Table[{c,f[c]}, {c,-60.,-12.,0.5}],
    Table[{c,f[c]}, {c,-12.,0.,0.05}]
  ];
  grid = DeleteDuplicatesBy[grid, First];
  chiMin = Min[grid[[All,2]]];
  best = First@MinimalBy[grid, Last];
  allowed = Select[grid, Last[#]-chiMin <= 3.84 &];
  lower = If[Length[allowed]>0, Min[allowed[[All,1]]], Missing["NoAllowedGridPoint"]];
  csv = FileNameJoin[{outdir,"highpt_"<>label<>"_L"<>ToString[Round[lam]]<>"GeV_v0_26.csv"}];
  Export[csv, Prepend[grid,{"C","chi2"}], "CSV"];
  {label,lam/1000.,best[[1]],chiMin,lower,3.84,csv}
];

rows = Flatten[Table[
  {scanOne["l2u2D22",lam],scanOne["l2d2D22",lam]},
  {lam,{3000.,4000.,5000.}}],1];
summary = Prepend[rows,{"label","Lambda_TeV","C_best","chi2_min","C95_lower_Cle0","DeltaChi2","grid_csv"}];
Export[FileNameJoin[{outdir,"highpt_tautau_d8_v0_26.csv"}],summary,"CSV"];
Print[Grid[summary,Frame->All]];
Print["STRICT_UNITARITY_POINT C = ",-8 Pi/3,"; compare to each C95_lower."];
Print["WROTE ",FileNameJoin[{outdir,"highpt_tautau_d8_v0_26.csv"}]];
