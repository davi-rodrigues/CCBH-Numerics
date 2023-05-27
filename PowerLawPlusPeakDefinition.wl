(* ::Package:: *)

(*ResourceFunction["DarkMode"][]*)


(*
  POWER-LAW-PLUS-PEAK DEFINITION
  ******************************* 
*)

\[Beta][m1_, \[Alpha]_, mMin_, mMax_] =  Piecewise[
  {
    {\[Beta]normalization m1^-\[Alpha], m1 <= mMax}, 
    {0, m1 > mMax}
  }
];

\[Beta]normalization = k /. First @ Solve[
  Integrate[ k m1^-\[Alpha], {m1, mMin, mMax}, Assumptions->{\[Alpha]>0, mMin >0, mMax > mMin}] == 1, 
  k
];

smoothing[m_, mMin_, \[Delta]m_] = Piecewise[{
  {0, m < mMin}, 
  {1/(ff[m-mMin, \[Delta]m]+1), mMin <= m < mMin+\[Delta]m}, 
  {1, m >= mMin + \[Delta]m}
}];

ff[m_,\[Delta]m_] = Exp[\[Delta]m/m + \[Delta]m/(m-\[Delta]m)];

ClearAll[powerLawPeak];
Options[powerLawPeak] = {
  \[Lambda] -> 0.10,
  \[Alpha] -> 2.63,
  mMin -> 4.59,
  \[Delta]m -> 4.82,
  mMax -> 86.22,
  \[Mu]m -> 33.07,
  \[Sigma]m -> 5.69
}; (*From Fig. 16 of GWTC-2 populations paper, 2010.14533*)

(*Options[powerLawPeak] = {
  \[Lambda] -> 0.05,
  \[Alpha] -> 3.3,
  mMin -> 4.59,
  \[Delta]m -> 4.82,
  mMax -> 100,
  \[Mu]m -> 33.07,
  \[Sigma]m -> 5.69
};
*) (*Values with Sumit*)

powerLawPeak[m1_, OptionsPattern[]] := Block[
  {
    \[Lambda] = OptionValue @ \[Lambda],
    \[Alpha] = OptionValue @ \[Alpha],
    mMin = OptionValue @ mMin,
    \[Delta]m = OptionValue @ \[Delta]m,
    mMax = OptionValue @ mMax,
    \[Mu]m = OptionValue @ \[Mu]m,
    \[Sigma]m = OptionValue @ \[Sigma]m,
    nBHs = 1,
    plp,
    normalization,
    m 
  },
  10^3 ((1 - \[Lambda]) \[Beta][m1, \[Alpha], mMin, mMax] + \[Lambda] PDF[NormalDistribution[\[Mu]m, \[Sigma]m], m1]) smoothing[m1, mMin, \[Delta]m]  
  (*Unnormalized version, faster, 10^3 is arbitrary. Normalization comes later.*)
];

plpNormalization = NIntegrate[powerLawPeak[m1], {m1, 0,200}]; (*The PDF becomes 0 after mMax, hence the normalization can be computed with high mMax.*)
\[ScriptCapitalD]plp = ProbabilityDistribution[powerLawPeak[m1]/ plpNormalization, {m1, 0, 100}];

Echo[NProbability[200>m>0, m \[Distributed] \[ScriptCapitalD]plp], "Probability of 0<m1<200: "];

(*Row[{Plot[PDF[\[ScriptCapitalD]plp, m1], {m1, 0, 100}, PlotRange->All, ImageSize-> Medium], LogPlot[PDF[\[ScriptCapitalD]plp, m1], {m1, 0, 100}, ImageSize-> Medium]}]*) 
