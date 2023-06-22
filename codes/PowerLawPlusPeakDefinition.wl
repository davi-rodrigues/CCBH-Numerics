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
piNorm[opts:OptionsPattern[plpp]] := piNorm[opts] = 1/NIntegrate[piUnnorm[m, opts], {m, 0,OptionValue@mMax}];
pi[m_, opts:OptionsPattern[plpp]] := pi[m, opts] = piNorm[opts] piUnnorm[m, opts];

Clear[dist, \[ScriptCapitalD]];
dist::usage = "dist[options] or \[ScriptCapitalD][options] from the PLPP context stands for the PLPP distribution";
dist[opts:OptionsPattern[plpp]] := dist[opts] = ProbabilityDistribution[pi[m, opts], {m, 0, OptionValue@mMax}];
\[ScriptCapitalD][opts:OptionsPattern[plpp]] := dist[opts];

Remove[lV, aV, mMinV, mMaxV, dmV, muV, sV];

End[];

optionsPlpp = Options[plpp`plpp];
Echo["","Use plpp`pi[m, options] for the power-law-plus-peak PDF and plpp`\[ScriptCapitalD][options] for the distribution. Use optionsPlpp to see the options."];
