(* ::Package:: *)

(* ::Title:: *)
(*CCBH and Power Law + Peak: Formation PDF for m2*)


(* ::Author:: *)
(*Davi  C. Rodrigues (this code).  2023.*)


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
(**)
(**)


(* ::Section:: *)
(*Starting...*)


SetDirectory[NotebookDirectory[]];
<<"Directories.wl";

(*
  Calling wl files
  ***************** 
*)
getCode["Constants.wl"];
baseSimPoints = 10^7; (*The base number of points to be simmulated for each dimension, commonly 10^5-10^7, depends on the purpose.*)

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



Off[NIntegrate::slwcon];
(*plpp`listPiM2; *)(*Takes about 15 minutes. For sure this can be faster, but it works. Parallelizing is not very helpful*)
(*dumpsave["plpplistPiM2.mx", plpp`listPiM2]*)
getAux["plpplistPiM2.mx"];
plpp`listPiM2


Plot[{plpp`pi[m1], plpp`piM2Inorm[m1]}, {m1,1,100}, PlotPoints->40, PlotRange->All]


#[RandomVariate[plpp`dist2[], 10^4]]  & /@  {Median, Mean} 
#[RandomVariate[plpp`dist[], 10^4]]  & /@  {Median, Mean} 


(* ::Section:: *)
(*Formation mass distribution*)


ClearAll[
  \[ScriptCapitalD]M1Formation, 
  \[ScriptCapitalD]M1FormationRaw
];

Clear[dataM2FormationRaw, dataM2Formation];
\[ScriptCapitalD]plp2 = plpp`\[ScriptCapitalD]2[]; 
dataM2Merger = RandomVariate[\[ScriptCapitalD]plp2, 2 baseSimPoints] ~ EchoTiming ~ "dataM2Merger"; 

dataM2FormationRaw[zObs_, tdmin_, tdmax_, k_, zMax_, w_]:= dataM2FormationRaw[zObs, tdmin, tdmax, k, zMax, w]= Block[
  {length, dataM2MergerSameLength},
  length = Length @ dataMassFactor[zObs, tdmin, tdmax, k, zMax, w];
  dataM2MergerSameLength = Take[dataM2Merger, length];
  1/dataMassFactor[zObs, tdmin, tdmax, 1, zMax, w]^k dataM2MergerSameLength 
];

dataM2Formation[zObs_, opts:OptionsPattern[optionstkw]] := dataM2FormationRaw[
  zObs, 
  OptionValue @ tdmin, 
  OptionValue @ tdmax, 
  OptionValue @ k, 
  OptionValue @ zMax, 
  OptionValue @ w
]


\[ScriptCapitalD]M1FormationRaw[zObs_, tdmin_, tdmax_, k_, zMax_, w_] := \[ScriptCapitalD]M1FormationRaw[zObs, tdmin, tdmax, k, zMax, w] = SmoothKernelDistribution[
  dataM1FormationRaw[zObs, tdmin, tdmax, k, zMax, w],
  0.05, 
  InterpolationPoints -> 10^5, 
  MaxExtraBandwidths -> {2,2}, (*This slightly extends the profile beyond data (about 0.1 solar mass). 
  It is only used for the plot. Use {0,0} for no extension.*)
  PerformanceGoal->"Quality"
];

\[ScriptCapitalD]M1Formation[zObs_, opts:OptionsPattern[optionstkw]] := \[ScriptCapitalD]M1FormationRaw[
  zObs, 
  OptionValue @ tdmin, 
  OptionValue @ tdmax, 
  OptionValue @ k, 
  OptionValue @ zMax, 
  OptionValue @ w
];

\[ScriptCapitalD]M2FormationRaw[zObs_, tdmin_, tdmax_, k_, zMax_, w_] := \[ScriptCapitalD]M2FormationRaw[zObs, tdmin, tdmax, k, zMax, w] = SmoothKernelDistribution[
  dataM2FormationRaw[zObs, tdmin, tdmax, k, zMax, w],
  0.05, 
  InterpolationPoints -> 10^5, 
  MaxExtraBandwidths -> {2,2}, (*This slightly extends the profile beyond data (about 0.1 solar mass). 
  It is only used for the plot. Use {0,0} for no extension.*)
  PerformanceGoal->"Quality"
];

\[ScriptCapitalD]M2Formation[zObs_, opts:OptionsPattern[optionstkw]] := \[ScriptCapitalD]M2FormationRaw[
  zObs, 
  OptionValue @ tdmin, 
  OptionValue @ tdmax, 
  OptionValue @ k, 
  OptionValue @ zMax, 
  OptionValue @ w
];

Clear[\[ScriptCapitalD]M2FormationEmpRaw, \[ScriptCapitalD]M2FormationEmp]
\[ScriptCapitalD]M2FormationEmpRaw[zObs_, tdmin_, tdmax_, k_, zMax_, w_] := \[ScriptCapitalD]M2FormationEmpRaw[zObs, tdmin, tdmax, k, zMax, w] = 
 EmpiricalDistribution[
   dataM2FormationRaw[zObs, tdmin, tdmax, k, zMax, w]
 ];

\[ScriptCapitalD]M2FormationEmp[zObs_, opts:OptionsPattern[optionstkw]] := \[ScriptCapitalD]M2FormationEmpRaw[
  zObs, 
  OptionValue @ tdmin, 
  OptionValue @ tdmax, 
  OptionValue @ k, 
  OptionValue @ zMax, 
  OptionValue @ w
];


Histogram[{dataM1Merger, dataM2Merger}, {1}, "PDF", PlotRange->All]


(* ::Section:: *)
(*Main results:*)


probabilityMminApprox2[mMin_, opts:OptionsPattern[optionstkw]] := SurvivalFunction[\[ScriptCapitalD]M2FormationEmp[0.4, opts], mMin]^69; (*Assumes 69 m2 BH masses.*)
probabilityMminDoubleApprox2[mMin_, opts:OptionsPattern[optionstkw]] := SurvivalFunction[\[ScriptCapitalD]M2FormationEmp[0.4, opts], mMin]^2*69;


{probabilityMminApprox2[2], findSigma@probabilityMminApprox2[2]}


0.00007159941350416846`^(72/69) // findSigma


{probabilityMminApprox2[2], findSigma@probabilityMminApprox2[2]}


{probabilityMminDoubleApprox2[2], findSigma@probabilityMminDoubleApprox2[2]}


tdmaxE = 5;
tdminE = 0.005;
zMaxE = 2;
\[ScriptCapitalD]M2Formation[0.4]; //EchoTiming
\[ScriptCapitalD]M2Formation[0.4, k-> 0.5]; //EchoTiming
\[ScriptCapitalD]M2Formation[0.4, tdmax -> tdmaxE]; //EchoTiming
\[ScriptCapitalD]M2Formation[0.4, tdmin -> tdminE]; //EchoTiming
\[ScriptCapitalD]M2Formation[0.4, zMax -> zMaxE]; //EchoTiming

legend1 = LineLegend[
  {Lighter@Orange, ColorData["AvocadoColors"][0.5], Lighter@Blue}, 
  {
    Style["\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)=3.0", FontFamily->"Times", FontSize -> 11], 
    Style["\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)=0.5", FontFamily->"Times", FontSize -> 11], 
    Style["\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)=0.0", FontFamily->"Times", FontSize -> 11]
  }, 
  LegendMarkerSize -> 20, 
  Spacings -> 0.15 (*It works, although in red*)
];

legend2 = LineLegend[
  {Directive[Lighter@Orange, Dashed], Directive[Lighter@Orange, Dotted], Directive[Lighter@Orange, DotDashed]}, 
  {
    Style["\!\(\*SubscriptBox[\(t\), \(max\)]\)=" <> ToString[Round@tdmaxE]<>" Gyr", FontFamily->"Times", FontSize -> 11], 
    Style["\!\(\*SubscriptBox[\(t\), \(min\)]\)=" <> ToString[Round[10^3 tdminE]]<>" Myr", FontFamily->"Times", FontSize -> 11],
    Style["\!\(\*SubscriptBox[\(z\), \(max\)]\)=" <> ToString[Round[zMaxE]], FontFamily->"Times", FontSize -> 11]
  }, 
  LegendMarkerSize -> 20, 
  Spacings -> 0.15 (*It works, although in red*)
];

tickLarge = 0.015;
tickSmall = 0.006;
frameTicksSmall = Flatten[Table[{n 10^i, "", tickSmall}, {n,1,10}, {i, -4, 0, 1}],1];
frameTicksLarge = Table[{10^i, Superscript[10, i], tickLarge}, {i, -4, 0, 1}];
frameTicksLargeDumb = Table[{10^i, "", tickLarge}, {i, -4, 0, 1}];
frameTicks = Flatten[Join[{frameTicksSmall, frameTicksLarge}],1];
frameTicksDumb = Flatten[Join[{frameTicksSmall, frameTicksLargeDumb}],1];


Off[General::munfl]; 
plotPDFformationDist2 = LogLogPlot[
  {
    PDF[\[ScriptCapitalD]M2Formation[0.4], x], 
    PDF[\[ScriptCapitalD]M2Formation[0.4, tdmax -> tdmaxE], x], 
    PDF[\[ScriptCapitalD]M2Formation[0.4, tdmin -> tdminE], x], 
    PDF[\[ScriptCapitalD]M2Formation[0.4, zMax -> zMaxE], x],
    PDF[\[ScriptCapitalD]M2Formation[0.4, k-> 0.5], x], 
    PDF[\[ScriptCapitalD]plp2, x]
  }, 
  {x, 0.5, plpp`mMax /. optionsPlpp}, 
  PlotRange->{{0.5, 500}, {8 10^-5, .3}},
  Background->White,
  Frame -> True,
  Axes -> False,
  FrameStyle-> Directive[Black, FontSize -> 13],
  FrameTicksStyle-> Directive[Black, FontFamily-> "Times"],
  FrameTicks->{{frameTicks, frameTicksDumb}, {Automatic, Automatic}},
  PlotPoints->5000,
  MaxRecursion->3,
  PlotStyle->{
    {Lighter@Orange}, 
    {Lighter@Orange, Dashed}, 
    {Lighter@Orange, Dotted}, 
    {Lighter@Orange, DotDashed}, 
    {ColorData["AvocadoColors"][0.5]}, 
    {Lighter@Blue}},
  PlotLegends-> {Placed[legend1, {0.60,0.81}], Placed[legend2, {0.85,0.81}]}
]
On[General::munfl];

exportOut["plotPDFformationDist2.pdf", plotPDFformationDist2]




