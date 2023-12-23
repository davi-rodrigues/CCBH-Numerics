(* ::Package:: *)

(*ResourceFunction["DarkMode"][]*)


(*
  POWER-LAW-PLUS-PEAK DEFINITION
  ******************************* 
  Davi C. Rodrigues (2023)
  
  Based on 
  * Talbot & Thrane ApJ 2018 (Measuring the binary black hole mass spectrum...)
  * Abbott et al ApJ 2021 (Population properties of Compact Objects...)
  * Abbott et al PRX 2023 (Population of merging...)
*)

Begin["plpp`"];

b[m_, a_, mMin_, mMax_] =  bNorm[a, mMin, mMax] m^-a;

bNorm[a_, mMin_, mMax_] = Block[{k}, 
  k /. First @ Solve[
    Integrate[ k m^-a, {m, mMin, mMax}, Assumptions->{a>0, mMin >0, mMax > mMin}] == 1, 
    k
  ]
];

smoothing[m_, mMin_, mMax_, dm_] = Piecewise[{
  {0, m < mMin}, 
  {0, m > mMax}, (*It is convenient to implement the maximum mass in this function.*)
  {1/(ff[m-mMin, dm]+1), mMin <= m < mMin+dm}, 
  {1, m >= mMin + dm}
}];

ff[m_,dm_] = Exp[dm/m + dm/(m-dm)];

Clear[plpp];
(* Options[plpp] = {
  a -> 2.63,
  mMin -> 4.59,
  dm -> 4.82,
  mMax -> 86.22,
  mu -> 33.07,
  s -> 5.69,
  l -> 0.10
}; From Fig. 16 of Abbott et al ApJ 2021 *)

Options[plpp] = {
  a -> 3.549,
  mMin -> 4.816,
  dm -> 5.454,
  mMax -> 83.140,
  mu -> 34.467,
  s -> 1.867,
  l -> 0.019,
  betaQ -> 0.760
}; (*Best values computed by Sumit Kumar based on GWTC-3, Abbott et al PRX 2023.*)

(*
  m1 PDF and distribution
  ***********************
*)

Clear[piUnnorm, piNorm, pi];
pi::usage = "pi[m, options] from the PLPP context is the (normalized) PDF of the power-law-plus-peak distribution.";
piUnnorm::usage = "piUnnorm[m, options] from the PLPP context is the unnormalized PDF of the power-law-plus-peak distribution.";
piNorm::usage = "piNorm[options] from the PLPP context is the normalization factor.";

piUnnorm[m_, OptionsPattern[plpp]] = Block[
  {
    lV = OptionValue @ l,
    aV = OptionValue @ a,
    mMinV = OptionValue @ mMin,
    dmV = OptionValue @ dm,
    mMaxV = OptionValue @ mMax,
    muV = OptionValue @ mu,
    sV = OptionValue @ s
  },
  10^3 ((1 - lV) b[m, aV, mMinV, mMaxV] + lV PDF[NormalDistribution[muV, sV], m]) smoothing[m, mMinV, mMaxV, dmV]
];
piNorm[opts:OptionsPattern[plpp]] := piNorm[opts] = 1/NIntegrate[
  piUnnorm[m, opts], 
  {m, 0, OptionValue@mMax},
  AccuracyGoal -> Infinity,
  PrecisionGoal -> 5,
  MaxRecursion -> 15
];
pi[m_, opts:OptionsPattern[plpp]] := pi[m, opts] = piNorm[opts] piUnnorm[m, opts];

Clear[dist, \[ScriptCapitalD]];
dist::usage = "dist[options] or \[ScriptCapitalD][options] from the PLPP context stands for the PLPP distribution";
dist[opts:OptionsPattern[plpp]] := dist[opts] = ProbabilityDistribution[pi[m, opts], {m, 0, OptionValue@mMax}];
\[ScriptCapitalD][opts:OptionsPattern[plpp]] := dist[opts];

(*
  m2 PDF and distribution
  ***********************
*)

piM2M1[m2_, m1_, opts:OptionsPattern[plpp]] := (m2/m1)^OptionValue@betaQ smoothing[m2, OptionValue@mMin, OptionValue@mMax, OptionValue@dm] HeavisideTheta[m1-m2]; (*PDF for M2 given M1*)

piM2M1normFactor[m1_?NumberQ]:= piM2M1normFactor[m1] = 1/NIntegrate[piM2M1[m2, m1], {m2, 1, 100},
  PrecisionGoal-> 3, 
  MaxRecursion -> 12, 
  AccuracyGoal -> 10,
  Method -> {Automatic, "SymbolicProcessing" -> 5}
]

piM2[m2_] := NIntegrate[piM2M1[m2, m1] piM2M1normFactor[m1] pi[m1], {m1, m2, 100}, 
  PrecisionGoal -> 3, 
  MaxRecursion -> 12, 
  AccuracyGoal -> 10, 
  Method -> {Automatic, "SymbolicProcessing" -> 0}
];

listPiM2 := listPiM2 = Map[{#, piM2[#]} &, Range[1, 100, 0.5]];

piM2I := piM2I = Interpolation[listPiM2];

piM2InormFactor := piM2InormFactor = 1/ NIntegrate[piM2I[m2], {m2, 1, 100}];

piM2Inorm[m2_] := piM2InormFactor piM2I[m2];


(* 

piM2Unnorm[m2_?NumberQ, opts:OptionsPattern[plpp]] := NIntegrate[
  10^3 piM2M1[m2, m1, opts] pi[m1, opts], {m1, 1, 100}, (*Requires mMin > 1 and mMax < 100.*)
  MinRecursion -> 3, 
  PrecisionGoal -> 5, 
  AccuracyGoal->Infinity, 
  Method-> {Automatic, "SymbolicProcessing" -> 0}, 
  MaxRecursion->20
];

piM2Norm[opts:OptionsPattern[plpp]] := piM2Norm[opts] =  1/NIntegrate[
  piM2Unnorm[m2], {m2, 1, 100}, 
  Method -> {Automatic, "SymbolicProcessing"->0}, 
  MaxRecursion->12, 
  PrecisionGoal->5
];

piM2[m2_, opts:OptionsPattern[plpp]] := piM2Unnorm[m2, opts] piM2Norm[opts];

piM2I[m2_, opts:OptionsPattern[plpp]] := piM2I[m2, opts] = Interpolation[
  {#, piM2[#, opts]} & /@ Range[1, 100, 0.025],
  InterpolationOrder -> 1
][m2]; *)

dist2[opts:OptionsPattern[plpp]] := dist2[opts] = ProbabilityDistribution[piM2Inorm[m2], {m2, 1, 100}];
\[ScriptCapitalD]2[opts:OptionsPattern[plpp]] := dist2[opts];


Remove[lV, aV, mMinV, mMaxV, dmV, muV, sV];

End[];

DistributeDefinitions["plpp`"];

optionsPlpp = Options[plpp`plpp];
Echo["","Use plpp`pi[m, options] for the power-law-plus-peak PDF and plpp`\[ScriptCapitalD][options] for the distribution. Use optionsPlpp to see the options."];
