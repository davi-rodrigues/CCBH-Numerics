(* ::Package:: *)

(* ::Title:: *)
(*CCBH and Power Law + Peak: Exclusion plots*)


(* ::Author:: *)
(*Davi  C. Rodrigues (this code). July, 2023.*)


(* ::Input:: *)
(*ResourceFunction["DarkMode"][]*)


(* ::Text:: *)
(*CCBH = Cosmologically coupled black hole (general case)*)
(*DEBH = Dark energy coupled Black Hole (implies k = - 3w)*)
(*GBH = Growing BH, k is arbitrary and without direct dark energy impact.*)
(*m1 = primary (larger) mass of the BBH or NSBH pair.*)
(**)
(*Naming conventions:*)
(*All defined variables/functions start with lower case letter.*)
(*"I" at the end of a function name: interpolated version. *)
(*"Raw": a function without options. It is useful to store results. The version with options is based on the raw version.*)


(* ::Section:: *)
(*Starting...*)


SetDirectory[NotebookDirectory[]];
Get[FileNameJoin[{"codes", "Directories.wl"}]];

(*
  Calling Code files
  ****************** 
*)
getCode["Constants.wl"];
baseSimPoints = 3 10^5; (*The base number of points to be simmulated for each dimension, commonly 10^5-10^7, depends on the purpose.*)

getCode["ObsDataPreparationGWTC-3.wl"]; (*Not essential for this notebook*)
getCode["PowerLawPlusPeakDefinition.wl"];
getCode["Cosmology.wl"];
getCode["MassFactorCorrection.wl"];

Print["End of running wl files."];

(*
  GENERAL PURPOSE FUNCTIONS
  **************************
*)

findSigma::usage = "findSigma[probability] yields the number of \[Sigma]'s for a given probability.";
findSigma[prob_] := Block[{sigma}, 
  sigma /. FindRoot[NProbability[Abs[x] > sigma , x\[Distributed]NormalDistribution[]] == prob, {sigma, 2}]
];

Clear[outSigma];
outSigma[n_] := outSigma[n] = NProbability[Abs[x] > n, x\[Distributed]NormalDistribution[]];



(* ::Section:: *)
(*Corner plots GBH and DEBH approaches (baseSimPoints = 3 x (10^5))*)


(* ::Subsection:: *)
(*Probability definition and testing*)


Clear[\[ScriptCapitalD]M1FormationEmpRaw,
  \[ScriptCapitalD]M1FormationEmp,
  probabilityMminApprox,
  probabilityMminDoubleApprox
];

position[list_, value_] := Position[list, value] /. {{x_}} -> x; (*useful for loading the data.*)

\[ScriptCapitalD]M1FormationEmpRaw[zObs_, tdmin_, tdmax_, k_, zMax_, w_] := \[ScriptCapitalD]M1FormationEmpRaw[zObs, tdmin, tdmax, k, zMax, w] = 
 EmpiricalDistribution[
   dataM1FormationRaw[zObs, tdmin, tdmax, k, zMax, w]
 ];

\[ScriptCapitalD]M1FormationEmp[zObs_, opts:OptionsPattern[optionstkw]] := \[ScriptCapitalD]M1FormationEmpRaw[
  zObs, 
  OptionValue @ tdmin, 
  OptionValue @ tdmax, 
  OptionValue @ k, 
  OptionValue @ zMax, 
  OptionValue @ w
];

probabilityMminApprox[mMin_, opts:OptionsPattern[optionstkw]] := SurvivalFunction[\[ScriptCapitalD]M1FormationEmp[0.4, opts], mMin]^72;
probabilityMminDoubleApprox[mMin_, opts:OptionsPattern[optionstkw]] := SurvivalFunction[\[ScriptCapitalD]M1FormationEmp[0.4, opts], mMin]^(2*72);


Echo[{probabilityMminApprox[2], findSigma[probabilityMminApprox[2]]}, "Probability and tension level in sigmas, current data: "];
Echo[{probabilityMminDoubleApprox[2], findSigma[probabilityMminDoubleApprox[2]]}, "Probability and tension level in sigmas, future data: "];
Echo[baseSimPoints, "baseSimPoints: "];


kValues = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.0}; 
wValues = Select[listw, -0.6 >= # >= -1&];
tdmaxValues = Range[3, 13, 1];
tdminValues = {0.005, 0.01, 0.02, 0.03, 0.04, 0.05}; 
zMaxValues = Range[2, 10, 1];


(* ::Subsection:: *)
(*Computing... And defining the probabilities*)


(* ::Text:: *)
(*Organized in 3 groups: *)
(*1. KT, KZ, KTmin*)
(*2. WT, WZ, WTmin*)
(*3. TZ, TTmin, ZTmin*)


(* ::Text:: *)
(*IMPORTANT: The reading of the files need to be updated. *)


(* ::Subsubsection:: *)
(*Computing KT, KZ, KTmin*)


isComputeTable\[ScriptCapitalD]M1KT = False; (*Takes about 4 minutes with baseSimPoints = 3 10^5*)

If[isComputeTable\[ScriptCapitalD]M1KT, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1KT = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> k1, tdmax-> tdmax1], 
    {k1, kValues}, 
    {tdmax1, tdmaxValues}
  ]; // EchoTiming;
  dumpsave["tableDM1KT.mx", table\[ScriptCapitalD]M1KT],
  (*else*)
  getAux["tableDM1KT.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> k1, tdmax-> tdmax1] = table\[ScriptCapitalD]M1KT[[position[kValues, k1], position[tdmaxValues, tdmax1]]], 
    {k1, kValues}, 
    {tdmax1, tdmaxValues}
  ];
  Echo["table\[ScriptCapitalD]M1KT loaded."];
];

tableProbKT = Table[
  {{k1, tdmax1}, probabilityMminApprox[2, k-> k1, tdmax -> tdmax1]}, 
  {k1, kValues}, 
  {tdmax1, tdmaxValues}
]; // EchoTiming

probabilityKT[k_, tdmax_] = Interpolation[Flatten[tableProbKT,1], InterpolationOrder -> 1][k, tdmax]; 


(* ::Text:: *)
(*The codes below do not fu*)


isComputeTable\[ScriptCapitalD]M1KZ = False; (*About 7 minutes with baseSimPoints = 3 10^5*)

If[isComputeTable\[ScriptCapitalD]M1KZ, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1KZ = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> k1, zMax-> zMax1], 
    {k1, kValues}, 
    {zMax1, zMaxValues}
  ]; // EchoTiming;
  dumpsave["tableDM1KZ.mx", table\[ScriptCapitalD]M1KZ],
  (*else*)
  getAux["tableDM1KZ.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> k1, zMax-> zMax1] = table\[ScriptCapitalD]M1KT[[position[kValues, k1], position[zMaxValues, zMax1]]], 
    {k1, kValues}, 
    {zMax1, zMaxValues}
  ];
  Echo["table\[ScriptCapitalD]M1KZ loaded."];
];

tableProbKZ = Table[
  {{k1, zMax1}, probabilityMminApprox[2, k-> k1, zMax-> zMax1]}, 
  {k1, kValues}, 
  {zMax1, zMaxValues}
]; // EchoTiming

probabilityKZ[k_, zMax_] = Interpolation[Flatten[tableProbKZ,1], InterpolationOrder -> 1][k, zMax]; 


isComputeTable\[ScriptCapitalD]M1KTmin = False;  (*About 4 min*)

If[isComputeTable\[ScriptCapitalD]M1KTmin, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1KTmin = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> k1, tdmin-> tdmin1], 
    {k1, kValues}, 
    {tdmin1, tdminValues}
  ]; // EchoTiming;
  dumpsave["tableDM1KTmin.mx", table\[ScriptCapitalD]M1KTmin],
  (*else*)
  getAux["tableDM1KTmin.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> k1, tdmin-> tdmin1] = table\[ScriptCapitalD]M1KT[[position[kValues, k1], position[tdminValues, tdmin1]]], 
    {k1, kValues}, 
    {tdmin1, tdminValues}
  ];
  Echo["table\[ScriptCapitalD]M1KTmin loaded."];
];

tableProbKTmin = Table[
  {{k1, tdmin1}, probabilityMminApprox[2, k-> k1, tdmin-> tdmin1]}, 
  {k1, kValues}, 
  {tdmin1, tdminValues}
]; // EchoTiming

probabilityKTmin[k_, tdmin_] = Interpolation[Flatten[tableProbKTmin,1], InterpolationOrder -> 1][k, tdmin]; 


(* ::Subsubsection:: *)
(*Computing: wT, wZ, wTmin, *)


isComputeTable\[ScriptCapitalD]M1WT = False; (*Takes about 3 minutes with baseSimPoints = 2 10^5*)

If[isComputeTable\[ScriptCapitalD]M1WT, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1WT = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> -3 w1, w -> w1, tdmax-> tdmax1], 
    {w1, wValues}, 
    {tdmax1, tdmaxValues}
  ]; // EchoTiming;
  dumpsave["tableDM1WT.mx", table\[ScriptCapitalD]M1WT],
  (*else*)
  getAux["tableDM1WT.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, w-> w1, tdmax-> tdmax1] = table\[ScriptCapitalD]M1KT[[position[wValues, w1], position[tdmaxValues, tdmax1]]], 
    {w1, wValues}, 
    {tdmax1, tdmaxValues}
  ];
  Echo["table\[ScriptCapitalD]M1WT loaded."];
];

tableProbWT = Table[
  {{w1, tdmax1}, probabilityMminApprox[2, k-> -3 w1, w -> w1, tdmax-> tdmax1]}, 
  {w1, wValues}, 
  {tdmax1, tdmaxValues}
]; // EchoTiming

probabilityWT[w_, tdmax_] = Interpolation[Flatten[tableProbWT,1], InterpolationOrder -> 1][w, tdmax]; 


isComputeTable\[ScriptCapitalD]M1WZ = False; (*About 5 minutes *)

If[isComputeTable\[ScriptCapitalD]M1WZ, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1WZ = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> - 3 w1, w -> w1, zMax-> zMax1], 
    {w1, wValues}, 
    {zMax1, zMaxValues}
  ]; // EchoTiming;
  dumpsave["tableDM1WZ.mx", table\[ScriptCapitalD]M1WZ],
  (*else*)
  getAux["tableDM1WZ.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, w-> w1, zMax-> zMax1] = table\[ScriptCapitalD]M1KT[[position[wValues, w1], position[zMaxValues, zMax1]]], 
    {w1, wValues}, 
    {zMax1, zMaxValues}
  ];
  Echo["table\[ScriptCapitalD]M1WZ loaded."];
];

tableProbWZ = Table[
  {{w1, zMax1}, probabilityMminApprox[2, k-> - 3 w1, w -> w1, zMax-> zMax1]}, 
  {w1, wValues}, 
  {zMax1, zMaxValues}
]; // EchoTiming

probabilityWZ[w_, zMax_] = Interpolation[Flatten[tableProbWZ,1], InterpolationOrder -> 1][w, zMax]; 


isComputeTable\[ScriptCapitalD]M1WTmin = False; (*2-3 min*)

If[isComputeTable\[ScriptCapitalD]M1WTmin, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1WTmin = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, k-> - 3 w1, w-> w1, tdmin-> tdmin1], 
    {w1, wValues}, 
    {tdmin1, tdminValues}
  ]; // EchoTiming;
  dumpsave["tableDM1WTmin.mx", table\[ScriptCapitalD]M1WTmin],
  (*else*)
  getAux["tableDM1WTmin.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, w-> w1, tdmin-> tdmin1] = table\[ScriptCapitalD]M1KT[[position[wValues, w1], position[tdminValues, tdmin1]]], 
    {w1, wValues}, 
    {tdmin1, tdminValues}
  ];
  Echo["table\[ScriptCapitalD]M1WTmin loaded."];
];

tableProbWTmin = Table[
  {{w1, tdmin1}, probabilityMminApprox[2, k-> - 3 w1, w-> w1, tdmin-> tdmin1]}, 
  {w1, wValues}, 
  {tdmin1, tdminValues}
]; // EchoTiming

probabilityWTmin[w_, tdmin_] = Interpolation[Flatten[tableProbWTmin,1], InterpolationOrder -> 1][w, tdmin]; 


(* ::Subsubsection:: *)
(*Computing: TZ, TTmin, ZTmin*)


isComputeTable\[ScriptCapitalD]M1TZ = False; (*6 min*)

If[isComputeTable\[ScriptCapitalD]M1TZ, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1TZ = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, tdmax-> tdmax1, zMax-> zMax1], 
    {tdmax1, tdmaxValues}, 
    {zMax1, zMaxValues}
  ]; // EchoTiming;
  dumpsave["tableDM1TZ.mx", table\[ScriptCapitalD]M1TZ],
  (*else*)
  getAux["tableDM1TZ.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, tdmax-> tdmax1, zMax-> zMax1] = table\[ScriptCapitalD]M1KT[[position[tdmaxValues, tdmax1], position[zMaxValues, zMax1]]], 
    {tdmax1, tdmaxValues}, 
    {zMax1, zMaxValues}
  ];
  Echo["table\[ScriptCapitalD]M1TZ loaded."];
];

tableProbTZ = Table[
  {{tdmax1, zMax1}, probabilityMminApprox[2, tdmax-> tdmax1, zMax-> zMax1]}, 
  {tdmax1, tdmaxValues}, 
  {zMax1, zMaxValues}
]; // EchoTiming

probabilityTZ[tdmax_, zMax_] = Interpolation[Flatten[tableProbTZ,1], InterpolationOrder -> 1][tdmax, zMax]; 


isComputeTable\[ScriptCapitalD]M1TTmin = False;  (*3-4 min*)

If[isComputeTable\[ScriptCapitalD]M1TTmin, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1TTmin = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, tdmax-> tdmax1, tdmin-> tdmin1], 
    {tdmax1, tdmaxValues}, 
    {tdmin1, tdminValues}
  ]; // EchoTiming;
  dumpsave["tableDM1TTmin.mx", table\[ScriptCapitalD]M1TTmin],
  (*else*)
  getAux["tableDM1TTmin.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, tdmax-> tdmax1, tdmin-> tdmin1] = table\[ScriptCapitalD]M1KT[[position[tdmaxValues, tdmax1], position[tdminValues, tdmin1]]], 
    {tdmax1, tdmaxValues}, 
    {tdmin1, tdminValues}
  ];
  Echo["table\[ScriptCapitalD]M1TTmin loaded."];
];

tableProbTTmin = Table[
  {{tdmax1, tdmin1}, probabilityMminApprox[2, tdmax-> tdmax1, tdmin-> tdmin1]}, 
  {tdmax1, tdmaxValues}, 
  {tdmin1, tdminValues}
]; // EchoTiming

probabilityTTmin[tdmax_, tdmin_] = Interpolation[Flatten[tableProbTTmin,1], InterpolationOrder -> 1][tdmax, tdmin]; 


isComputeTable\[ScriptCapitalD]M1ZTmin = False; (*5 min*)

If[isComputeTable\[ScriptCapitalD]M1ZTmin, 
  Echo["Computing..."];
  table\[ScriptCapitalD]M1ZTmin = Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, zMax-> zMax1, tdmin-> tdmin1], 
    {zMax1, zMaxValues}, 
    {tdmin1, tdminValues}
  ]; // EchoTiming;
  dumpsave["tableDM1ZTmin.mx", table\[ScriptCapitalD]M1ZTmin],
  (*else*)
  getAux["tableDM1ZTmin.mx"];
  Table[
    \[ScriptCapitalD]M1FormationEmp[0.4, zMax-> zMax1, tdmin-> tdmin1] = table\[ScriptCapitalD]M1KT[[position[zMaxValues, zMax1], position[tdminValues, tdmin1]]], 
    {zMax1, zMaxValues}, 
    {tdmin1, tdminValues}
  ];
  Echo["table\[ScriptCapitalD]M1ZTmin loaded."];
];

tableProbZTmin = Table[
  {{zMax1, tdmin1}, probabilityMminApprox[2, zMax-> zMax1, tdmin-> tdmin1]}, 
  {zMax1, zMaxValues}, 
  {tdmin1, tdminValues}
]; // EchoTiming

probabilityZTmin[zMax_, tdmin_] = Interpolation[Flatten[tableProbZTmin,1], InterpolationOrder -> 1][zMax, tdmin]; 


(* ::Subsection:: *)
(*Plots*)


(* ::Subsubsection:: *)
(*Plots of k vs probability *)


frameTicksSmall = Flatten[Table[{n 10^i, "", tickSmall}, {n,1,10}, {i, -10, 0, 1}],1];
frameTicksLarge = Table[{10^i, Superscript[10, i], tickLarge}, {i, -10, 0, 1}];
frameTicks = Join[frameTicksSmall, frameTicksLarge];
tickLarge = 0.015;
tickSmall = 0.006;


plotkProbability = Block[{},
  LogPlot[
    {
      probabilityKT[k, tdmax0],
      probabilityKT[k, tdmax0]^2
    },
    {k, 0.5, 3},
    PlotPoints -> 20,
    MaxRecursion->0,
    PlotStyle->{{Thick, Black}, {Thick, Black, Dashed}},
    Background->White,
    Frame -> True,
    Axes -> False,
    FrameStyle-> Directive[Black, FontSize -> 13],
    FrameTicksStyle-> Directive[Black, FontFamily-> "Times"],
    FrameTicks->{{frameTicks, {{outSigma[1], "1\[Sigma]"}, {outSigma[2], "2\[Sigma]"}, {outSigma[3], "3\[Sigma]"}, {outSigma[4], "4\[Sigma]"}, {outSigma[5], "5\[Sigma]"}}}, {Automatic, Automatic}},
    GridLines->{None, {outSigma[1], outSigma[2], outSigma[3], outSigma[4], outSigma[5]}},
    GridLinesStyle->Directive[Gray, Dotted],
    PlotLegends-> Placed[
      {Style["current data", FontFamily->"Times", FontSize -> 14], Style["forecast", FontFamily->"Times", FontSize -> 14]}, 
      {0.3,0.3}
    ],
    PlotRange->{All, {2 10^-7, 2}}
  ]
]

(*exportOut["plotkProbability.pdf", plotkProbability]*)


FindRoot[(probabilityKT[k, tdmax0] //findSigma) == 2, {k, 1, 2}] // Quiet


(* ::Text:: *)
(*Now with Luca's function*)


pointsLuca = {{0.2`,0.9999999999979087`},{0.30000000000000004`,0.9999999999979087`},{0.4`,0.9999999999979087`},{0.5`,0.9979114858401865`},{0.6000000000000001`,0.9913076700861714`},{0.7`,0.9813534788396552`},{0.8`,0.9536970773327864`},{0.9000000000000001`,0.8958486726075061`},{1.`,0.8281748880287615`},{1.1`,0.7557367813610051`},{1.2`,0.6780132988000289`},{1.3`,0.5959458338676502`},{1.4000000000000001`,0.511524384871613`},{1.5`,0.42684725162977033`},{1.6`,0.34658463728522426`},{1.7`,0.276431967376146`},{1.8`,0.21638817172690022`},{1.9000000000000001`,0.1672760102091132`},{2.`,0.12851258171569166`},{2.1`,0.09820129368041103`},{2.2`,0.0748284009264067`},{2.3000000000000003`,0.0569030518896383`},{2.4000000000000004`,0.043213343356207384`},{2.5000000000000004`,0.03279174072861294`},{2.6000000000000005`,0.02487653404872536`},{2.7`,0.018874566684654762`},{2.8000000000000003`,0.014327865885673271`},{2.9000000000000004`,0.010885180158910268`},{3.`,0.008278492986420966`}};
fullptabLuca[k_] = Interpolation[pointsLuca][k];

frameTicksSmall = Flatten[Table[{n 10^i, "", tickSmall}, {n,1,10}, {i, -10, 0, 1}],1];
frameTicksLarge = Table[{10^i, Superscript[10, i], tickLarge}, {i, -10, 0, 1}];
frameTicks = Join[frameTicksSmall, frameTicksLarge];
tickLarge = 0.015;
tickSmall = 0.006;


plotkProbabilityDirectPLPP = Block[{},
  LogPlot[
    {
      fullptab[k, 13],
      fullptab[k,13]^2,
      probabilityKT[k, tdmax0],
      probabilityKT[k, tdmax0]^2
    },
    {k, 0.5, 3},
    PlotPoints -> 20,
    MaxRecursion->0,
    PlotStyle->{{Thick, ColorData["AvocadoColors"][0.5]}, {Thick, ColorData["AvocadoColors"][0.5], Dashed}, {Thick, Lighter@Blue}, {Thick, Lighter@Blue, Dashed}},
    Background->White,
    Frame -> True,
    Axes -> False,
    FrameStyle-> Directive[Black, FontSize -> 13],
    FrameTicksStyle-> Directive[Black, FontFamily-> "Times"],
    FrameTicks->{{frameTicks, {{outSigma[1], "1\[Sigma]"}, {outSigma[2], "2\[Sigma]"}, {outSigma[3], "3\[Sigma]"}, {outSigma[4], "4\[Sigma]"}, {outSigma[5], "5\[Sigma]"}}}, {Automatic, Automatic}},
    GridLines->{None, {outSigma[1], outSigma[2], outSigma[3], outSigma[4], outSigma[5]}},
    GridLinesStyle->Directive[Gray, Dotted],
    PlotLegends-> Placed[
      {Style["Direct, current data", FontFamily->"Times", FontSize -> 13], Style["Direct, forecast", FontFamily->"Times", FontSize -> 13], Style["PLPP, current data", FontFamily->"Times", FontSize -> 13], Style["PLPP, forecast", FontFamily->"Times", FontSize -> 13]}, 
      {0.27,0.25}
    ],
    PlotRange->{All, {2 10^-7, 2}}
  ]
]

exportOut["plotkProbabilityDirectPLPP.pdf", plotkProbabilityDirectPLPP]


(* ::Subsubsection:: *)
(*Exclusion plots*)


largeTick = 0.02;
smallTick = 0.01;

(*The ticks definitions:*)
frameTicksK = Table[If[IntegerPart[i]==i, {i/10., i/10., largeTick}, {i/10., "", smallTick}], {i, 0, 30, 2.5}];
frameTicksW = Table[If[IntegerPart[i]==i, {i/10., i/10., largeTick}, {i/10., "", smallTick}], {i, -10, -6.0, 0.5}];
frameTicksZ = Table[If[IntegerPart[i]==i, {i, i, largeTick}, {i, "", smallTick}], {i, 2., 10, 0.5}];
frameTicksTmin = Table[If[IntegerPart[i]==i, {i/100, i 10 (*Myr unit*), largeTick}, {i/100, "", smallTick}], {i, 0.0, 5, 0.5}];
frameTicksTmax = Table[If[IntegerPart[i]==i, {i, i, largeTick}, {i, "", smallTick}], {i, 3, 13, 0.5}];

(*Here the "dumb" ticks*)
frameTicksKD = frameTicksK /. {x_, y_, z_} -> {x, "", z}; 
frameTicksWD = frameTicksW /. {x_, y_, z_} -> {x, "", z}; 
frameTicksZD = frameTicksZ /.{x_, y_, z_} -> {x, "", z};
frameTicksTminD = frameTicksTmin /. {x_, y_, z_} -> {x, "", z};
frameTicksTmaxD = frameTicksTmax /.{x_, y_, z_} -> {x, "", z};

(*Sizes and spaces definitions*)
plotSize = 200;
spaceLeft = 30;
spaceBottom = 20;
spaceBetween = 10;

(*Setting general options for the RegionPlots*)
frameStyle = Directive[Black, FontSize -> 16, Thickness[0.005]];
frameTicksStyle = Directive[Black, FontFamily-> "Times"];
boundaryStyle = {1 -> Directive[White, Thin], 2 -> Directive[White, Thin], 3-> Directive[White, Thin]};
plotStyle = {Directive[Lighter@Blue, Opacity[0.3]], Directive[Lighter@Blue, Opacity[0.5]], Directive[Blue, Opacity[0.5]]};

SetOptions[RegionPlot, 
  PlotStyle -> plotStyle, 
  BoundaryStyle -> boundaryStyle, 
  FrameTicksStyle -> frameTicksStyle,
  FrameStyle -> frameStyle,
  PlotRangePadding -> None,
  Background -> White
];

(*The plots...*)

regionPlotKT2 = RegionPlot[
  {
    probabilityKT[k,tdmax] < outSigma[1],
    probabilityKT[k,tdmax] < outSigma[2],
    probabilityKT[k,tdmax] < outSigma[3]
  },
  {k, .5, 3}, {tdmax, First@tdmaxValues, tdmax0},
  FrameTicks -> {{frameTicksTmax, frameTicksTmaxD}, {frameTicksKD, frameTicksKD}},
  PlotRange -> {{0,3}, {First@tdmaxValues, tdmax0}},
  ImagePadding -> {{spaceLeft, spaceBetween}, {spaceBetween, spaceBetween}},
  ImageSize -> {plotSize +  spaceBetween + spaceLeft, plotSize + 2 spaceBetween}
]

regionPlotKZ2 = RegionPlot[
  {
    probabilityKZ[k,zMax] < outSigma[1],
    probabilityKZ[k,zMax] < outSigma[2],
    probabilityKZ[k,zMax] < outSigma[3]
  },
  {k, .5, 3}, {zMax, First@zMaxValues, zMax0},
  FrameTicks->{{frameTicksZ, frameTicksZD}, {frameTicksKD, frameTicksKD}},
  FrameTicksStyle-> Directive[Black, FontFamily-> "Times"],
  ImagePadding->{{spaceLeft, spaceBetween}, {spaceBetween, spaceBetween}},
  ImageSize->{plotSize + spaceBetween + spaceLeft, plotSize + 2 spaceBetween},
  PlotRange->{{0,3}, {First@zMaxValues, zMax0}}
]

regionPlotKTmin2 = RegionPlot[
  {
    probabilityKTmin[k,tdmin] < outSigma[1],
    probabilityKTmin[k,tdmin] < outSigma[2],
    probabilityKTmin[k,tdmin] < outSigma[3]
  },
  {k, .5, 3}, {tdmin, First@tdminValues, tdmin0},
  FrameTicks -> {{frameTicksTmin, frameTicksTminD},{frameTicksK, frameTicksKD}},
  ImagePadding->{{spaceLeft, spaceBetween}, {spaceBottom, spaceBetween}},
  ImageSize->{plotSize + spaceBetween + spaceLeft, plotSize + spaceBottom + spaceBetween},
  PlotLegends -> Placed[{
    Style["1\[Sigma]", FontFamily->"Times", FontSize -> 18], 
    Style["2\[Sigma]", FontFamily->"Times", FontSize -> 18], 
    Style["3\[Sigma]", FontFamily->"Times", FontSize -> 18]}, 
    {0.15, 0.18}
  ],
  PlotRange->{{0,3}, {First@tdminValues, tdmin0}}
]




regionPlotWT2 = RegionPlot[
  {
    probabilityWT[w,tdmax] < outSigma[1],
    probabilityWT[w,tdmax] < outSigma[2],
    probabilityWT[w,tdmax] < outSigma[3]
  },
  {w,-1, -0.6}, {tdmax, First@tdmaxValues, tdmax0},
  FrameTicks -> {{frameTicksTmaxD, frameTicksTmaxD}, {frameTicksWD, frameTicksWD}},
  PlotRange -> {{-0.6,-1}, {First@tdmaxValues, tdmax0}},
  ImagePadding -> {{spaceBetween, spaceBetween}, {spaceBetween, spaceBetween}},
  ImageSize -> {plotSize + 2 spaceBetween, plotSize + 2 spaceBetween},
  ScalingFunctions -> {"Reverse", Identity}
]

regionPlotWZ2 = RegionPlot[
  {
    probabilityWZ[w,zMax] < outSigma[1],
    probabilityWZ[w,zMax] < outSigma[2],
    probabilityWZ[w,zMax] < outSigma[3]
  },
  {w, -0.6, -1}, {zMax, First@zMaxValues, zMax0},
  FrameTicks->{{frameTicksZD, frameTicksZD}, {frameTicksWD, frameTicksWD}},
  FrameTicksStyle-> Directive[Black, FontFamily-> "Times"],
  ImagePadding->{{spaceBetween, spaceBetween}, {spaceBetween, spaceBetween}},
  ImageSize->{plotSize + 2 spaceBetween, plotSize + 2 spaceBetween},
  PlotRange->{{-1, -0.6}, {First@zMaxValues, zMax0}},
  ScalingFunctions -> {"Reverse", Identity}
]

regionPlotWTmin2 = RegionPlot[
  {
    probabilityWTmin[w,tdmin] < outSigma[1],
    probabilityWTmin[w,tdmin] < outSigma[2],
    probabilityWTmin[w,tdmin] < outSigma[3]
  },
  {w, -0.6, -1}, {tdmin, First@tdminValues, tdmin0},
  FrameTicks -> {{frameTicksTminD, frameTicksTminD},{frameTicksW, frameTicksWD}},
  ImagePadding->{{spaceBetween,  spaceBetween}, {spaceBottom, spaceBetween}},
  ImageSize->{plotSize + 2 spaceBetween, plotSize + spaceBottom + spaceBetween},
  PlotRange->{{-1,-0.6}, {First@tdminValues, tdmin0}},
  ScalingFunctions -> {"Reverse", Identity}
]



regionPlotTZ2 = RegionPlot[
  {
    probabilityTZ[tdmax,zMax] < outSigma[1],
    probabilityTZ[tdmax,zMax] < outSigma[2],
    probabilityTZ[tdmax,zMax] < outSigma[3]
  },
  {tdmax, First@tdmaxValues, tdmax0}, {zMax, First@zMaxValues, zMax0},
  FrameTicks->{{frameTicksZD, frameTicksZD},{frameTicksTmaxD, frameTicksTmaxD}},
  ImagePadding->{{spaceBetween, spaceBetween}, {spaceBetween, spaceBetween}},
  ImageSize->{plotSize + 2 spaceBetween, plotSize + 2 spaceBetween},
  PlotRange->{{First@tdmaxValues, tdmax0}, {First@zMaxValues, zMax0}}
]

regionPlotTTmin2 = RegionPlot[
  {
    probabilityTTmin[tdmax,tdmin] < outSigma[1],
    probabilityTTmin[tdmax,tdmin] < outSigma[2],
    probabilityTTmin[tdmax,tdmin] < outSigma[3]
  },
  {tdmax, First@tdmaxValues, tdmax0}, {tdmin, First@tdminValues, tdmin0},
  FrameTicks->{{frameTicksTminD, frameTicksTminD},{frameTicksTmax, frameTicksTmaxD}},
  ImagePadding->{{spaceBetween, spaceBetween}, {spaceBottom, spaceBetween}},
  ImageSize->{plotSize + 2 spaceBetween, plotSize + spaceBottom + spaceBetween},
  PlotRange->{{First@tdmaxValues, tdmax0}, {First@tdminValues, tdmin0}}
]

regionPlotZTmin2 = RegionPlot[
  {
    probabilityZTmin[zMax,tdmin] < outSigma[1],
    probabilityZTmin[zMax,tdmin] < outSigma[2],
    probabilityZTmin[zMax,tdmin] < outSigma[3]
  },
  {zMax, First@zMaxValues, zMax0}, {tdmin, First@tdminValues, tdmin0},
  FrameTicks->{{frameTicksTminD, frameTicksTminD},{frameTicksZ, frameTicksZD}},
  ImagePadding->{{spaceBetween, spaceBetween}, {spaceBottom, spaceBetween}},
  ImageSize->{plotSize + 2 spaceBetween, plotSize + spaceBottom + spaceBetween},
  PlotRange->{{First@zMaxValues, zMax0}, {First@tdminValues, tdmin0}}
]


cornerPlotCCBH = Grid[
  {
    {regionPlotKT2, regionPlotWT2},
    {regionPlotKZ2, regionPlotWZ2, regionPlotTZ2}, 
    {regionPlotKTmin2, regionPlotWTmin2, regionPlotTTmin2, regionPlotZTmin2}
  }, 
  Spacings -> {0,0}, 
  Background->White
]

exportOut["cornerPlotCCBH2.jpg", cornerPlotCCBH]





