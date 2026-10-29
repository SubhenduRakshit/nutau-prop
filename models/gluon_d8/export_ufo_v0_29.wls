(* Export the minimal Q_l2G2D[33] collider model to UFO.
   Usage:
     FEYNRULES_PATH=/path/to/FeynRules wolframscript -file models/gluon_d8/export_ufo_v0_29.wls
*)

fr = Environment["FEYNRULES_PATH"];
If[fr === $Failed || fr === "",
  Print["ERROR: set FEYNRULES_PATH to the FeynRules installation directory."]; Exit[2]
];

repo = DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
model = FileNameJoin[{repo,"models","gluon_d8","GLUON_D8.fr"}];
out = FileNameJoin[{repo,"models","gluon_d8","SM_GLUON_D8_UFO"}];
sm = FileNameJoin[{fr,"Models","SM","SM.fr"}];
If[!FileExistsQ[sm], sm = FileNameJoin[{fr,"Models","SM.fr"}]];

SetDirectory[fr];
<< FeynRules`;
FeynmanGauge = False;
LoadModel[sm, model];

Print["=== Hermiticity check: BSM term ==="];
CheckHermiticity[LGLUOND8];

Print["=== Requested gg tau tau vertex ==="];
verts = FeynmanRules[LGLUOND8,
  SelectParticles -> {{G,G,ta,tabar}},
  ScreenOutput -> False
];
Print[verts];
If[Length[verts] < 1, Print["ERROR: gg tau tau vertex not generated"]; Exit[3]];

If[DirectoryQ[out], DeleteDirectory[out, DeleteContents -> True]];
WriteUFO[LSM + LGLUOND8, Output -> out];
Print["WROTE ", out];
