(* ::Package:: *)

(*
  Constants.wl
  This code is part of CCBHminMass-PLPP
  Davi C. Rodrigues (2023)
*)

(*
  PHYSICAL CONSTANTS
  ******************
*)

(*\[CapitalOmega]\[CapitalLambda] = 0.6889;*) (*2018 data from Planck*)
(*\[CapitalOmega]m = 0.315;*) (*2018 data from Planck*)
(*H0Gy = 67.66 / (kpc 10^3) gyr;*) (*1/Gyr, 2018 data from Planck, 67.66 km/s Mpc^-1*)

\[CapitalOmega]m = 0.32;
H0Gy = 70.0 / (kpc 10^3) gyr; (* With these rounded values, the universe age is 13.2 GY*)
gyr = 10^9 60 60 24 365; (*Gyr in seconds*)
kpc = 3.08568 10^16; (*kpc in km*)

(* ------ *)
w0 = -1; (*w0 is the standard value for w, the DE equation of state parameter*)
k0 = 3; (*k0 is the standard growth exponent for CBH.*)
tdmax0 = 13.0; (*Gyrs. The standard maximum delay time. 
13.5 Gyrs from GWTC-3 2111.03634v4 (population paper), but we adopt here a slightly smaller value*)
tdmin0 = 0.05; (*Gyrs. The standard minimum delay time. From GWTC-3 population paper*)
zMax0 = 10; (* Highest redshift where BBH start to form. In accordance with GWTC-3, arXiv:2111.03634v4 (population paper) *)


(*
  META PARAMETERS AND OPTIONS
  ****************************
*)

optionstkw = {
  tdmin -> tdmin0, 
  tdmax -> tdmax0, 
  k -> k0, 
  w -> w0, 
  zMax -> zMax0
};

Clear[baseSimPoints, proportionalBaseSimPoints];
proportionalBaseSimPoints[tdmin_, tdmax_] = Round[(tdmax - tdmin)/(tdmax0 - tdmin0/10) baseSimPoints]; 
(*The maximum number of realizations is baseSimPoints, 
but this number can get smaller proportionally to tdmax - tdmin. *)

$HistoryLength = 5; (*Reduces memory usage*)

(*
  WARNING:
  baseSimPoints need to be defined AFTER calling this wl file.
  Different notebooks commonly use different baseSimPoints.
*)
