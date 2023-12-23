(* 
  GWTC-3 Data Preparation
  This is part of the Cosmologically Coupled Black Holes code
  
  Davi C. Rodrigues (2023)
*)


(*
  IMPORTING DATA
  **************
*)

(*Events and classification from 2111.03634 with FAR < 1/year*)
datasetGWpop = Import[FileNameJoin[{pathIn, "GWlistPopulationExtended.csv"}],"Dataset"];

(*GWTC-3, data from https://www.gw-openscience.org/eventapi/html/GWTC *)
datasetGWTC = Import[FileNameJoin[{pathIn, "GWTC.csv"}], "Dataset", "HeaderLines"->{1,1}];


(*
  PREPARING THE DATA
  ******************
*)

(*
  PREPARING THE DATA: Removing unnecessary data and tagging the BHNS events
*)

(*Echo[
  Complement[datasetGWpop[All,1] // Normal, gwtcEvents], 
  "Event that is in gwPopEvents, but not in gwtcEvents: "
];*)

pos200105 = Position[datasetGWpop,"GW200105_162426"] // First // First;

(* Echo[
  pos200105 = Position[datasetGWpop,"GW200105_162426"] // First // First, 
  "Event that is in gwPopEvents, but not in gwtcEvents (pos200105): "
]; *)

pos190426 = Position[datasetGWpop,"GW190426_152155"] // First // First;

(* Echo[
  pos190426 = Position[datasetGWpop,"GW190426_152155"] // First // First, 
  "Event that is in gwPopEvents, but not in gwtcEvents (pos190426): "
]; *)

(*
  There are two events that are in 2111.03634, but not in the gtw-3 provided data. 
  According to 2111.03634 (Population of mergin...) these events have a low p_astro, of only 0.36 and 0.14. 
  There are no other events with p_astro in this list with p_astro < 0.49. One is a BHNS event, the other BBH.
*)

bnsEvents = Cases[datasetGWpop, {_, "BNS"} ][All,1] // Normal;
posBNSlines = Block[{pos}, 
  pos = Position[datasetGWpop, "BNS"];
  Partition[pos[[All,1]], 1]
];
(* Echo[bnsEvents, "BNS events: "];
Echo[posBNSlines, "BNS line positions in datasetGWpop: "]; *)

bhnsEvents = Cases[datasetGWpop, {_, "BHNS"} ][All,1] // Normal;
posBHNSlines = Block[{pos}, 
  pos = Position[datasetGWpop, "BHNS"];
  Partition[pos[[All,1]], 1]
];
(* Echo[bhnsEvents, "BHNS events: "];
Echo[posBHNSlines, "BHNS line positions in datasetGWpop: "]; *)

posNeglected = Join[posBNSlines, {{pos200105}, {pos190426}}];  (*These two events do not appear in GWTC-3, since their p_astro is too small*)

(* Echo[
  posNeglected = Join[posBNSlines, {{pos200105}, {pos190426}}],  (*These two events do not appear in GWTC-3, since their p_astro is too small*)
  "datasetGWpop lines to be neglected: "
]; *)

(*
  PREPARING THE DATA: Defining the relevant dataset by removing the unnecessary events
*)

datasetGWpopSelec = Delete[datasetGWpop, posNeglected];
bhnsEventsSelec = Cases[datasetGWpopSelec, {_, "BHNS"} ][All,1] // Normal;
posBHNSselectLines = Block[{pos}, 
  pos = Position[datasetGWpopSelec, "BHNS"];
  Partition[pos[[All,1]], 1]
];

(*
  We will be interested only on the GWTC-3 events that appear in datasetGWpopSelec.
  To select such, we need to better prepare the datasetGWT key names, since datasetGWTC presents the events names with version labels.
  In the following, we remove such version labels and select the data from datasetGWTC that corresponds to the events in datasetGWpopSelec.
*)

datasetGWTCnoV = datasetGWTC[KeyMap[#/. x_ :>  StringReplace[x, "-v"~~_ -> ""] & , #]&];
datasetGWTrelev = datasetGWTCnoV[datasetGWpopSelec[All,1] // Normal] 

dataGWexport = datasetGWTrelev[All, 
  {
    "redshift", 
    "redshift_lower", 
    "redshift_upper", 
    "mass_1_source", 
    "mass_1_source_lower", 
    "mass_1_source_upper", 
    "mass_2_source", 
    "mass_2_source_lower", 
    "mass_2_source_upper"
  }
];

(* exportOut["dataGWexport.csv", dataGWexport]; *) (*Uncomment this to export the table.*)

(*
  OBSERVATIONAL DATA: PREPARING AND PLOTTING
  ******************************************
*)

mXzData1 = datasetGWTrelev[All, {"redshift", "redshift_lower","redshift_upper", "mass_1_source","mass_1_source_lower", "mass_1_source_upper"}][Values, Values] // Normal;
mXzDataForPlotting1 = {Around[#1, {#2, #3}], Around[#4, {#5, #6}]} & @@@ mXzData1;

(*For mXzData2, we Do NOT at this point remove the NS data from the BHNS event, use Delete[...., posBHNSselectLines] to do that.*)

mXzData2= datasetGWTrelev[All, {"redshift", "redshift_lower","redshift_upper" , "mass_2_source","mass_2_source_lower", "mass_2_source_upper"}][Values, Values] // Normal;
mXzDataForPlotting2= {Around[#1,{#2, #3}],Around[#4,{#5,#6}]}& @@@ mXzData2;

Echo[Length @ mXzDataForPlotting1, "Number of m1 black holes: "];



