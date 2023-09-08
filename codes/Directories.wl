
(*
  Directories.wl
  This is part of the CCBHminMass-PLPP code
  Davi C. Rodrigues (2023)
*)

(* 
  DIRECTORY STRUCTURE
  *******************
*)
pathBase = Directory[];
pathOut = FileNameJoin[{pathBase, "output"}];
pathIn =  FileNameJoin[{pathBase, "input"}];
pathAux = FileNameJoin[{pathBase, "auxiliary"}];
pathCodes = FileNameJoin[{pathBase, "codes"}];

exportOut[fileName_, variable_] := (
  Export[FileNameJoin[{pathOut, fileName}], variable];
  "pathOut/"<>fileName
);

exportAux[fileName_, variable_] := (
  Export[FileNameJoin[{pathAux, fileName}], variable];
  "pathAux/"<>fileName
);

SetAttributes[dumpsave, HoldRest];
dumpsave[fileName_, variable_] := (
  DumpSave[FileNameJoin[{pathAux, fileName}], variable];
  "pathAux/"<>fileName
);

getAux[fileName_] := Get[FileNameJoin[{pathAux, fileName}]];

getCode[fileName_] := (
  Print["Starting "<>ToString[fileName]<> ": "]; 
  Get[FileNameJoin[{pathCodes, fileName}]]
);

SetDirectory[pathBase];