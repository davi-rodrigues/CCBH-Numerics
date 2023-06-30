(* 
  GWTC-3 Data Preparation
  This is part of the Cosmologically Coupled Black Holes code
  
  Davi C. Rodrigues (2023)
*)


(*
  IMPORTING DATA
  **************
*)


(*Events and classification from 2112.06861, only events with FAR < 0.25 / year:*)
datasetGWpop = Import[FileNameJoin[{pathIn, "GWlistPopulation.csv"}], "Dataset"];

(*The extended version considers events with FAR < 1/year*)
datasetGWpopExtended = Import[FileNameJoin[{pathIn, "GWlistPopulationExtended.csv"}],"Dataset"];

(*GWTC-3, data from https://www.gw-openscience.org/eventapi/html/GWTC/:*)
datasetGWTC = Import[FileNameJoin[{pathIn, "GWTC.csv"}], "Dataset", "HeaderLines"->{1,1}];

datasetGWpop = datasetGWpopExtended; (*Only uses the extended version*)


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

(* mzListPlot[x__] := ListPlot[x, 
  PlotRange->All, 
  Axes -> False, 
  Frame -> True, 
  PlotMarkers->Graphics[{Black, Opacity[0.6],Disk[]}, ImageSize -> 11, ImagePadding->1], 
  IntervalMarkersStyle-> {Opacity[0.5],Darker@Blue},
  ImageSize->500, 
  GridLines->Automatic, 
  GridLinesStyle->Dotted,
  FrameTicksStyle->Directive[FontFamily->"Times New Roman",FontSize->17, FontColor -> "Black"],
  Background -> White
]; *)

(* plot1=datasetGWTrelev[All, {"redshift", "mass_1_source"}]//Values // mzListPlot;
plot2 = mzListPlot[datasetGWTrelev[All, {"redshift", "mass_2_source"}]//Values , 
  PlotMarkers->Graphics[{Blue, Opacity[0.6],Disk[]}, ImageSize -> 11, ImagePadding->1]
]; *)
(* Show[plot1,plot2] *)


(*
  From this point, the m2 masses will not be used.
  To include m2 masses, uncomment the mXzData2 code part below.
*)

mXzData1 = datasetGWTrelev[All, {"redshift", "redshift_lower","redshift_upper", "mass_1_source","mass_1_source_lower", "mass_1_source_upper"}][Values, Values] // Normal;
mXzDataForPlotting1 = {Around[#1, {#2, #3}], Around[#4, {#5, #6}]} & @@@ mXzData1;

(*For mXzData2, we Do NOT remove the NS data from the BHNS event, use Delete[...., posBHNSselectLines] to do that.*)

mXzData2= datasetGWTrelev[All, {"redshift", "redshift_lower","redshift_upper" , "mass_2_source","mass_2_source_lower", "mass_2_source_upper"}][Values, Values] // Normal;
mXzDataForPlotting2= {Around[#1,{#2, #3}],Around[#4,{#5,#6}]}& @@@ mXzData2;
(* mXzDataForPlotting = Join[mXzDataForPlotting1,mXzDataForPlotting2]; *)


(* mXzDataForPlotting = mXzDataForPlotting1; (*Selects only m1 data.*)
plotZxM = mzListPlot[
  mXzDataForPlotting, 
  (* FrameLabel-> {Style["z", FontFamily->"Times", 20, Italic], Style["M", FontFamily->"Times", 20, Italic]}, *) (*Labels will be generated in the tex file*)
  IntervalMarkersStyle -> {Opacity[0.3], Gray},
  PlotRange-> {{-0.01,1.0},{-1, 110}}
]; *)
Echo[Length @ mXzDataForPlotting, "Number of m1 black holes: "];

(* Print[plotZxM]; *)


(*
  DATA BINNING AND PLOTTING (probably unnecessary)
  *************************
*)

(* I will consider bins of size 0.12, which is about the median value of the uncertainties.*)

binWidth = 0.125;

bin1 = Select[mXzDataForPlotting , #[[1,1]] < binWidth & ] // Sort;
bin[1] = bin1;
bin[i_] := Select[mXzDataForPlotting, binWidth (i-1)< #[[1,1]] < binWidth i & ] // Sort ;

listBins = {1,2,3,4,5,6, 8}; (*The bin numbers that have data *)

(* OLD approach: does not use the observational errors and considers StandardDeviation for the uncertainties:
  binCentral[i_]:=bin[i][[All,1;;2,1]]; (*It is bin[i], but without the uncertanties*)
  listBinnedMass = Mean[binCentral[#]] & /@ listBins;

  listStd=(StandardDeviation[binCentral[#]] & /@ Drop[listBins, -1] )[[All,2]]; 
  (*Since bin[7] only has one data point, StandardDeviation cannot be applied to it, that's why a Drop is used.*)
  listStd = Append[listStd, Mean[bin[7][[1, 2, 2]]]];
  (*In the above, the corresponding StandardDeviation for bin[7] is taken to be the Mean of the upper and lower uncertainties.*)

  binnedDataWerrors = MapThread[{#1[[1]], Around[#1[[2]], #2]} & , {listBinnedMass, listStd}];
*)

binnedDataWerrors = {First @ Mean[bin[#][[All, 1]]],  Mean[bin[#][[All, 2]]]} & /@ listBins; (* MeanAround provides a weighted evaliation of the mean. 
For Mean and MeanAround, error bars show the error on the mean (assuming that there is a common value for all these data), not the data dispersion*)

plotZxMbinned = mzListPlot[
  binnedDataWerrors, 
  PlotMarkers -> Graphics[{Red,Opacity[0.6],Thickness[0.15],Circle[]}, ImageSize -> 17, ImagePadding->2], 
  IntervalMarkersStyle -> {Opacity[1],Thickness[0.002], Red}, 
  Joined -> True, 
  PlotStyle -> {Thickness[0.002], Red}
];

plotZxMwithBinned = Show[plotZxM, plotZxMbinned]; (*Seems that this binned version will not be relevant*)


