
(*
  Cosmology
  This is part of PowerLawPlusPeakCCBH code
  Davi C. Rodrigues (2023)
*)


(*
  COSMOLOGY
  *********
*)

(* ::Section:: *)
(* Cosmology *)

Clear[zFormation, hubble, timeZ, w, z];

hubble[z_, w_] = H0Gy Sqrt[\[CapitalOmega]m (1+z)^3+(1-\[CapitalOmega]m)(1+z)^(3(1+w))]; 

Clear[timeZ];
timeZ[z1_?NumberQ, z2_?NumberQ, w_] := NIntegrate[
  ((1+z) hubble[z, w])^-1, 
  {z, z1, z2}, 
  Method-> {Automatic, "SymbolicProcessing" -> 0}, (*Improves computational time by a factor 5-10, without changes in the answers.*)
  PrecisionGoal -> 5,
  AccuracyGoal -> Infinity,
  MaxRecursion -> 10  (*Removing this may improve speed, it is necessary to check*)
];

zFormation[zObs_?NumberQ, delayTime_?NumberQ, zMax_, w_] := If[
  timeZ[zObs, zMax, w] > delayTime, (*Time between zObs and zMax > td*)
  (*then*)
  FindRoot[
    timeZ[zObs, zformation, w] == delayTime, 
    {zformation, zObs, zMax}, 
    PrecisionGoal -> 5, (*Reminder: AccuracyGoal behaves differently in FindRoot. This is just a minimum precision requirement, probably it will be higher.*)
    Method -> "Brent" (*This is the automatic method for this case, but I am stating it explicitly.*)
  ][[1,2]],
  (*else*)
  zMax  
];

(*
  COSMOLOGY: Defining zFormationI 
*)

(* ::Subsection:: *)
(*Defining zFormationI *)

isComputeTableZformation = False; (*Change to True to compute and regenerate the mx file. About 10 min.*)

(*Even if the above is False, the table is computed in case the file tableZformation.mx does not exist.*)
If[
   FileExistsQ[FileNameJoin[{pathAux, "tableZformation.mx"}]],
   Null,
   (*else*)
   Echo["File tableZformation.mx not found. It will be computed."];
   isComputeTableZformation = True
];

If[isComputeTableZformation, 
  Echo["Computing tableZformation..."];
  tableZformation = ParallelTable[
    {{zObs, delayTime, w1}, zFormation[zObs, delayTime, zMax0, w1]}, (*this table uses zMax0, zFornationI can restrict that*)
    {zObs, listzObs},
    {delayTime, listDelayTime}, (*It is computed up to tdmax0*)
    {w1, listw}
  ]; // EchoTiming;
  tableZformationFlat = Flatten[tableZformation, 2];
  dumpsave["tableZformation.mx", tableZformationFlat],
  (*else*)
  getAux["tableZformation.mx"];
  Echo["tableZformation loaded."];
];

Clear[zFormationIraw, zFormationIfast, zFormationI]
zFormationInoZmax[zObs_, delayTime_, w_] = Interpolation[
  tableZformationFlat,
  InterpolationOrder -> 1 
][zObs, delayTime, w]; (*the "I" in the name stands for interpolated. "raw" for no options.*)

SetAttributes[zFormationI, Listable]; (*zFormationInoZmax is automatically listable, but the Min function is not Listable*)
zFormationI[zObs_, delayTime_, zMax_, w_] = Min[
  zFormationInoZmax[zObs, delayTime, w],
  zMax
];
