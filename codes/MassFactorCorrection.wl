
(*
  The Mass Factor Correction and Formation Mass
  This is part of PowerLawPlusPeakCCBH code
  Davi C. Rodrigues (2023)
  
  DEPENDENCIES:
  - Cosmology.wl 
  - PowerLawPlusPeakDefinition.wl
*)

(*
  WARNING ABOUT PARALLELIZING: 
  Computing anything tha depends directly on zFormationI should not be done in parallel.
  If necessary, one needs to define zFormationI directly from `Interpolation' in the subkernels. 
  Even doing the above, the speed was not good (tested in one computer only, may be due to a memory issue).
  Compmuting zFormation in parallel is not worth, since zFormationI is ~100 times faster.
*)

(*
  MASS FACTOR
  ***********
*)

(* ::Section:: *)
(* Mass Factor *)

ClearAll[
  massFactorDelay, 
  dataMassFactor,
  dataDelayTime,
  minDelayTime, 
  maxDelayTime
];

Echo[baseSimPoints, "Base number of simmulated points per dimension (baseSimPoints): "];

listDelayTime = Range[tdmin0/10, tdmax0, 0.015]; (* delayTime values to be computed. tdmin lower than tdmin0 will be considered*)
listzObs = Range[0.0, 1.0, 0.01];
listw = Range[-1.5, -0.5, 0.1]; 

minDelayTime = First[listDelayTime]; (*min delay time to be computed, it can be lower than tdmin0*)
maxDelayTime[zObs_, tdmax_?NumberQ, zMax_?NumberQ, w_] := maxDelayTime[zObs, tdmax, zMax, w] = Min[
  tdmax,
  timeZ[zObs, zMax, w]
]; 

massFactor[zObs_,zFormation_, k_] = ((1 + zFormation)/(1 + zObs))^k; (*(a/a_i)^k. This is the factor that corrects the mass*)

SetAttributes[massFactorDelay, Listable];
massFactorDelay[zObs_, delayTime_, k_, zMax_, w_] =  massFactor[zObs, zFormationI[zObs, delayTime, zMax, w], k];

randomLogDelayTime[zObs_, tdmin_, tdmax_, zMax_, w_, realizations_] := randomLogDelayTime[zObs, tdmin, tdmax,  zMax, w, realizations] = 
 RandomReal[{Log10[tdmin], Log10 @ maxDelayTime[zObs, tdmax, zMax, w]}, realizations]; (*log-flat distribution.*)

dataDelayTime[zObs_, tdmin_, tdmax_, zMax_, w_] := dataDelayTime[zObs, tdmin, tdmax, zMax, w] = 
 10^randomLogDelayTime[zObs, tdmin, tdmax, zMax, w,  proportionalBaseSimPoints[tdmin, tdmax]]; 

dataMassFactor[zObs_, tdmin_, tdmax_, k_, zMax_, w_] := dataMassFactor[zObs, tdmin, tdmax, k, zMax, w] = 
 massFactorDelay[zObs, dataDelayTime[zObs, tdmin, tdmax, zMax, w], k, zMax, w];
 
(*
  FORMATION MASS
  ***********
*)

(* ::Section:: *)
(* Formation Mass *)

Clear[dataM1FormationRaw, dataM1Formation];
\[ScriptCapitalD]plp = plpp`\[ScriptCapitalD][]; (*plpp`\[ScriptCapitalD] is defined in PowerLawPlusPeakDefinition.wl.*)
dataM1Merger = RandomVariate[\[ScriptCapitalD]plp, 2 baseSimPoints] ~ EchoTiming ~ "dataM1Merger"; 
(*
  "2 baseSimPoints" is important to guarantee that the number of dataM1Merger points will be larger than the number of points 
  of dataMassFactor, which depends on tdmin, tdmax.
*)

dataM1FormationRaw[zObs_, tdmin_, tdmax_, k_, zMax_, w_]:= dataM1FormationRaw[zObs, tdmin, tdmax, k, zMax, w]= Block[
  {length, dataM1MergerSameLength},
  length = Length @ dataMassFactor[zObs, tdmin, tdmax, k, zMax, w];
  dataM1MergerSameLength = Take[dataM1Merger, length];
  1/dataMassFactor[zObs, tdmin, tdmax, 1, zMax, w]^k dataM1MergerSameLength 
];(*Note that the k=1 inside dataMassFactor, but there is k in the exponent. It's relevant for function memmory, 
such that dataMassFactor needs not to be comptuted again if k changes.*)

dataM1Formation[zObs_, opts:OptionsPattern[optionstkw]] := dataM1FormationRaw[
  zObs, 
  OptionValue @ tdmin, 
  OptionValue @ tdmax, 
  OptionValue @ k, 
  OptionValue @ zMax, 
  OptionValue @ w
]
