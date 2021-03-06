(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



(* ::Input::Initialization:: *)
GenerateSPhenoCTcouplings:=Block[{i,currentRegime,readRegime,tt2},
Print["   Writing Counter-Terms for Couplings"];

sphenoCoup=OpenWrite[ToFileName[$sarahCurrentSPhenoDir,"CouplingsCT_"<>ModelName<>".f90"]];

WriteCopyRight[sphenoCoup];

WriteString[sphenoCoup,"Module CouplingsCT_"<>ModelName<>"\n \n"];
WriteString[sphenoCoup, "Use Control \n"];
WriteString[sphenoCoup, "Use Settings \n"];
WriteString[sphenoCoup, "Use Model_Data_"<>ModelName<>" \n"];
WriteString[sphenoCoup, "Use Mathematics, Only: CompareMatrices, Adjungate \n \n"];

WriteString[sphenoCoup, "Contains \n \n "];

InitCouplingCT;

WriteSPhenoAllCouplings[SPhenoCoupCT,parametersAllCT,namesAllCT,"CalculateCouplingCT",""];
WriteSPhenoCouplings[SPhenoCoupCT,False,""];

WriteString[sphenoCoup,"End Module CouplingsCT_"<>ModelName<>" \n"];
Close[sphenoCoup];

];

InitCouplingCT:=Block[{i,j,temp,cnames,pos,ctname,ctsub},
namesAllCT={};
parametersAllCT=Flatten[Join[parametersAll,parametersAll/.subListCounterNamesSingle]];
SPhenoCoupCT={};

For[i=1,i<=Length[SPhenoCouplingsAll],
cnames=Drop[SPhenoCouplingsAll[[i,2]],1];
ctsub={SPhenoCouplingsAll[[i,2,1]]->"CT_"<>SPhenoCouplingsAll[[i,2,1]]};
For[j=1,j<=Length[cnames],
pos=Position[SPhenoParameters,cnames[[j]]][[1,1]];
ctname=ToExpression["ct"<>ToString[cnames[[j]]]];
SPhenoParameters=Join[SPhenoParameters,{SPhenoParameters[[pos]] /. cnames[[j]]->ctname}];
namesAllCT=Join[namesAllCT,{ctname}];
ctsub=Join[ctsub,{cnames[[j]]->ctname}];
j++;];
Switch[Length[SPhenoCouplingsAll[[i]]],
6,
SPhenoCoupCT=Join[SPhenoCoupCT,{{SPhenoCouplingsAll[[i,1]],SPhenoCouplingsAll[[i,2]]/.ctsub,SPhenoCouplingsAll[[i,3]],Flatten[Join[SPhenoCouplingsAll[[i,4]],SPhenoCouplingsAll[[i,4]] /.subListCounterNamesSingle]],ExpandCT[SPhenoCouplingsAll[[i,5]]],SPhenoCouplingsAll[[i,6]]}}];,
7,
SPhenoCoupCT=Join[SPhenoCoupCT,{{SPhenoCouplingsAll[[i,1]],SPhenoCouplingsAll[[i,2]]/.ctsub,SPhenoCouplingsAll[[i,3]],Flatten[Join[SPhenoCouplingsAll[[i,4]],SPhenoCouplingsAll[[i,4]] /.subListCounterNamesSingle]],ExpandCT[SPhenoCouplingsAll[[i,5]]],ExpandCT[SPhenoCouplingsAll[[i,6]]],SPhenoCouplingsAll[[i,7]]}}];
];
i++;];

];

ExpandCT[vert_]:=Block[{temp},
temp=Expand[Expand[TrigExpand[vert] /. subListCounter]-TrigExpand[vert]]/. conj[dd[a_]]->dd[conj[a]]/. sum->SUM //. SUM[a__, d_ SUM[f__,g_]]->SUM[a,SUM[f,d g]]//. SUM[a_,b_,c_,d_Times]:>SUM[a,b,c,Expand[d]] //.SUM[a_,b_,c_,d_ + e_]:>SUM[a,b,c,d]+SUM[a,b,c,e]//. sum[a_,b_,c_,d_Times]:>sum[a,b,c,Expand[d]] //.sum[a_,b_,c_,d_ + e_]:>sum[a,b,c,d]+sum[a,b,c,e];
temp=Expand[temp]//. dd[a__] dd[b__]->0 //.dd[a__]^b__->0//. dd[a_] SUM[___,y___ dd[b_]]->0 //. dd[a_] SUM[__,SUM[___,y___ dd[b_]]]->0//. SUM[__, y___ dd[a_]] SUM[___,w___ dd[b_]]->0 //. SUM[__, dd[a_]] SUM[___,  dd[b_]]->0 //. SUM[__,SUM[__, x_ dd[a_]]] SUM[__,SUM[___,y_  dd[b_]]]->0/. SUM->sum;
If[Head[temp]===Plus,temp=Plus@@(Select[(List@@temp),Count[#,dd[x_],99]<2 &]);];
Return[temp /.dd[conj[a_]]->conj[dd[a]]/. subListCounterNames];
];



GenerateQuarticsGaugeShift:=Block[{i,j,temp,cnames,pos,ctname,ctsub,neededc,allg},
allg=Transpose[BetaGauge][[1]];
subListCounter=Table[allg[[i]]->allg[[i]]+dd[allg[[i]]],{i,1,Length[allg]}];
subListCounterNames={dd[x_]:>ToExpression["d"<>ToString[x]]};
subListCounterNamesSingle=Table[allg[[i]]->(ToExpression["d"<>ToString[allg[[i]]]]),{i,1,Length[allg]}];

SPhenoParameters=Join[SPhenoParameters,Table[{(ToExpression["d"<>ToString[allg[[i]]]]),{},{}},{i,1,Length[allg]}]];
realVar=Join[realVar,Table[(ToExpression["d"<>ToString[allg[[i]]]]),{i,1,Length[allg]}]];

namesQGS={};
parametersQGS=Flatten[Join[parametersAll,Transpose[BetaGauge][[1]]/.subListCounterNamesSingle]];
SPhenoCoupQGS={};

neededc=Intersection[Select[Flatten[Transpose[HighscaleMatchingConditions][[2]]/.Times->List/.Plus->List] ,Head[#]==Dot&]/.Dot->C /.A_[b_]:>A /;(FreeQ[PART[S],A]==False  && A=!=List)];
neededc=Select[SPhenoCouplingsAllreallyAll,FreeQ[neededc,#[[1,1]]]==False&];


For[i=1,i<=Length[neededc],
cnames=Drop[neededc[[i,2]],1];
ctsub={neededc[[i,2,1]]->"MSDR_"<>neededc[[i,2,1]]};
For[j=1,j<=Length[cnames],
pos=Position[SPhenoParameters,cnames[[j]]][[1,1]];
ctname=ToExpression["d"<>ToString[cnames[[j]]]];
SPhenoParameters=Join[SPhenoParameters,{SPhenoParameters[[pos]] /. cnames[[j]]->ctname}];
namesQGS=Join[namesQGS,{ctname}];
ctsub=Join[ctsub,{cnames[[j]]->ctname}];
j++;];
Switch[Length[neededc[[i]]],
6,
SPhenoCoupQGS=Join[SPhenoCoupQGS,{{neededc[[i,1]],neededc[[i,2]]/.ctsub,neededc[[i,3]],Intersection[Flatten[Join[neededc[[i,4]],neededc[[i,4]] /.subListCounterNamesSingle]]],ExpandCT[neededc[[i,5]]],neededc[[i,6]]}}];,
7,
SPhenoCoupQGS=Join[SPhenoCoupQGS,{{neededc[[i,1]],neededc[[i,2]]/.ctsub,neededc[[i,3]],Intersection[Flatten[Join[neededc[[i,4]],neededc[[i,4]] /.subListCounterNamesSingle]]],ExpandCT[neededc[[i,5]]],ExpandCT[neededc[[i,6]]],neededc[[i,7]]}}];
];
i++;];


];

