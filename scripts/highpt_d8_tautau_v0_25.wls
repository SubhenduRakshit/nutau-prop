(* HighPT Run-II di-tau recast for the active d=8 UHE operators.
   Exact HighPT WC labels:
     charm:   WC["l2u2D22",{3,3,2,2}]
     strange: WC["l2d2D22",{3,3,2,2}]
   Positivity branch used in the UHE project: C <= 0.

   IMPORTANT: EFTorder->4 gives the consistent O(Lambda^-4) SMEFT
   prediction, i.e. SM x d=8 interference. It intentionally does NOT
   include (d=8)^2 = O(Lambda^-8).
*)

Quiet[Needs["HighPT`"]];
$PrintingProcessInfo = False;
$ParallelHighPT = True;

outdir = FileNameJoin[{Directory[], "results", "v0.25"}];
If[!DirectoryQ[outdir], CreateDirectory[outdir, CreateIntermediateDirectories -> True]];

scanOne[label_, lam_] := Module[
  {wc, chiExpr, f, grid, chiMin, allowed, bound, best, csv},
  InitializeModel["SMEFT", EFTorder -> 4, OperatorDimension -> 8, EFTscale -> lam];
  wc = WC[label, {3, 3, 2, 2}];
  chiExpr = Total@ChiSquareLHC["di-tau-ATLAS",
      Coefficients -> {wc}, EFTorder -> 4, OperatorDimension -> 8,
      EFTscale -> lam, RescaleError -> True];
  f[c_?NumericQ] := N[chiExpr /. {wc -> c, Conjugate[wc] -> c}];
  grid = Table[{c, f[c]}, {c, -500., 0., 0.25}];
  chiMin = Min[grid[[All, 2]]];
  best = First@MinimalBy[grid, Last];
  allowed = Select[grid, Last[#] - chiMin <= 3.84 &];
  bound = If[Length[allowed] > 0, Min[allowed[[All, 1]]], Missing["NoAllowedGridPoint"]];
  csv = FileNameJoin[{outdir, "highpt_" <> label <> "_L" <> ToString[Round[lam]] <> "GeV.csv"}];
  Export[csv, Prepend[grid, {"C", "chi2"}], "CSV"];
  {label, lam/1000., best[[1]], chiMin, bound, 3.84, csv}
];

rows = Flatten[Table[
   {scanOne["l2u2D22", lam], scanOne["l2d2D22", lam]},
   {lam, {3000., 4000., 5000.}}], 1];

summary = Prepend[rows,
  {"label", "Lambda_TeV", "C_best", "chi2_min", "C95_lower_positive_branch_Cle0", "DeltaChi2", "grid_csv"}];
Export[FileNameJoin[{outdir, "highpt_tautau_d8_v0_25.csv"}], summary, "CSV"];
Print[Grid[summary, Frame -> All]];
Print["WROTE ", FileNameJoin[{outdir, "highpt_tautau_d8_v0_25.csv"}]];
