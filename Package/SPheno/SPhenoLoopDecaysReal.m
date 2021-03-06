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
getbremsgeneral[e1_,p1_,p2_,gaugegroupnr_]:=Block[{i,j,k,ps,cps,pos,toordered,coup,colstr,colindices,colindices1,colindices2,colindices3,colindices4,treecf,cf12,cf13,cf23,otherindices,otherfactor,totherfactor,cgsub,Casvec,Dynvec,colpair,colordimp,Colourpresent,dimvec,tdimvec,conjcoup,gen1,gen2,Coloureasy,cin1,cin2,cin3,cin4,cin5,cin6,cin7,ctvlist,cttarlist,ctrepr,tpos,cinlist,ggname,ggtype,pindices,othergroups,sumline,Casarr,jmap,colstrord,treetosum},
tdimvec={1,0,0};
Casarr=Array[0&,{6}];
cps={AntiField[e1],p1,p2};
ps={e1,p1,p2};
pindices=getIndizesWI/@{e1,p1,p2};
ggname=Gauge[[gaugegroupnr,3]];
ggtype=Gauge[[gaugegroupnr,2]];
If[FreeQ[pindices,ggname],Return[{0,Casarr}];];
cgsub={CG[SU[3],{{0,1},{1,0},{1,1}}][a_,b_,c_]->Lam[c,b,a]/2,CG[SU[3],{{0,1},{1,1},{1,0}}][a_,b_,c_]->Lam[b,c,a]/2,CG[SU[3],{{1,1},{0,1},{1,0}}][a_,b_,c_]->Lam[a,c,b]/2,CG[SU[3],{{1,0},{0,1},{1,1}}][a_,b_,c_]->Lam[c,a,b]/2,CG[SU[3],{{1,0},{1,1},{0,1}}][a_,b_,c_]->Lam[b,a,c]/2,CG[SU[3],{{1,1},{1,0},{0,1}}][a_,b_,c_]->Lam[a,b,c]/2};
(*First get the vertex*)
pos=Position[VerticesInv[All],C[AntiField[e1],p1,p2]];
If[pos==={},Vertex::DoesNotExist="The Vertex `` does not exist!";
Message[Vertex::DoesNotExist,C[AntiField[e1],p1,p2]];
Abort[];,pos=pos[[1,1]];];
toordered=VerticesOrg[All][[pos]]/.A_[{b__}]->A;
ctvlist={ct1,ct2,ct3};
cinlist=Array[Symbol["cin"<>ToString[#]]&,{3}];
cttarlist=Table[{cps[[i]],cinlist[[i]]},{i,1,3}];

ctrepr={};
For[i=1,i<=3,i++,
tpos=Position[cttarlist,toordered[[i]],2][[1,1]];(* the 2 here tells it to look only at level 2 \[Rule] otherwise we match conj[Su] to Su etc *)

(*Print[toordered[[i]],tpos,cttarlist];*)
AppendTo[ctrepr,ctvlist[[i]]->cttarlist[[tpos,2]]];
cttarlist=Delete[cttarlist,tpos];
];
coup=VerticesVal[All][[pos]][[1]];

(*For all unbroken nonabelian symmetries*)
otherindices=DeleteCases[DeleteCases[pindices,{generation,_},3],{lorentz,_},3];
othergroups=Intersection[Transpose[Flatten[otherindices,1]][[1]]];
treecf=1;
For[i=1,i<=Length[othergroups],i++,
	totherfactor=ExtractStructure[coup,othergroups[[i]]][[1,1]]/.{Lam[a__]->Lam[a]/2}/.cgsub;
	sumline={};
	tdimvec[[1]]=1;
	For[j=1,j<=3,j++,
	tpos=Position[otherindices[[j]],othergroups[[i]]];
	If[tpos!={},tpos=tpos[[1,1]];
	tdimvec[[j]]=otherindices[[j,tpos,2]];
	If[j>1,AppendTo[sumline,{Symbol["cin"<>ToString[j]],1,tdimvec[[j]]}];
	];
	];
	];

	If[othergroups[[i]]==ggname,
	dimvec=tdimvec;
	colstr=totherfactor;
	];

	treetosum=(totherfactor*conj[totherfactor])/.ctrepr/.{cin1->1} ;
	treecf=treecf*Sum[treetosum,Evaluate[Sequence@@sumline]];
	]; (* end for[i=1,...*)
Casvec=Array[0&,{3}];
Casarr=Array[0&,{6}];
Coloureasy=False;
For[i=1,i<=3,i++,
	If[dimvec[[i]]>1,
	Colourpresent=True;
	Casvec[[i]]=SA`Casimir[RE[ps[[i]]],gaugegroupnr];,Coloureasy=True;
	];
	];
jmap={1,4,6};
For[i=1,i<=3,i++,
	Casarr[[jmap[[i]]]]=treecf*Casvec[[i]];
	];
	If[Coloureasy,
		(*
		For[i=1,i\[LessEqual]2,i++,
		For[j=i+1,j\[LessEqual]3,j++,
		If[Casvec[[i]]=Casvec[[j]],
		Casarr[[i]]=treecf*Casvec[[i]]
		]
		];

		];
*)
	If[Casarr[[4]]==Casarr[[6]],
	Casarr[[5]]=-Casarr[[4]],Casarr[[5]]=0];

Casarr[[2]]=Casarr[[5]]+Casarr[[4]];
Casarr[[3]]=Casarr[[5]]+Casarr[[6]];


	Return[{treecf,Casarr}];
	];(* if[Coloureasy*)

(*If here then all states coloured\[Rule]more complicated*)
conjcoup=conj[colstr]/.ctrepr/.{cin1->1,cin2->cin4,cin3->cin5};
colstrord=colstr/.ctrepr/.cin1->1;
gen1=Generator[ggtype,SA`DynL[p1,gaugegroupnr]][cin6,cin2,cin4];
gen2=Generator[ggtype,SA`DynL[p2,gaugegroupnr]][cin6,cin3,cin5];
(*cf12=Sum[colstrord*conjcoup*gen1*gen2,{cin2,1,dimvec[[2]]},{cin3,1,dimvec[[3]]},{cin4,1,dimvec[[2]]},{cin5,1,dimvec[[3]]},{cin6,1,getNumberStatesAdjoint[ggtype]}];
cf13=-Casarr[[1]]-cf12;
cf23=Conjugate[-Casarr[[4]]-cf12];*)

(*NB,ordering of the diagrams is 1\[Rule]V (photon or gluon) on incoming line; 2\[Rule]V on particle 2,3\[Rule]V on particle 3*)cf23=Sum[colstrord*conjcoup*gen1*gen2,{cin2,1,dimvec[[2]]},{cin3,1,dimvec[[3]]},{cin4,1,dimvec[[2]]},{cin5,1,dimvec[[3]]},{cin6,1,getNumberStatesAdjoint[ggtype]}];
(*NB,-M1+M2+M3=0 due to incoming vs outgoing particles*)
       (*so M13=M23+M33,M12=M22+M32*)

cf13=cf23+Casarr[[6]];
       cf12=Conjugate[cf23+Casarr[[4]]];

Casarr[[2]]=cf12;
Casarr[[3]]=cf13;
Casarr[[5]]=cf23;
Return[{treecf,Casarr}]];
MakeRealCorrections[NModule_,DecayingParticle_,twobodydecays_,OneLtbds_,type_]:=Block[{i,j,k,replaceexternalparticles,NRoutine,treecolourfactor,Casarr,Q2arr,Casvec,Chargevec,ps,colourpresent,chargepresent,coloreasy,chargeeasy,colordimp},
replaceexternalparticles={f1->DecayingParticle,f2->OneLtbds[[1]],f3->OneLtbds[[2]]};
NRoutine=NModule<>"To"<>ToString[If[(ToString[Head[f2/.replaceexternalparticles]]==ToString[conj])||(ToString[Head[f2/.replaceexternalparticles]]==ToString[bar]),c,""]]<>ToString[RE[f2/.replaceexternalparticles]]<>ToString[If[(ToString[Head[f3/.replaceexternalparticles]]==ToString[conj])||(ToString[Head[f3/.replaceexternalparticles]]==ToString[bar]),c,""]]<>ToString[RE[f3/.replaceexternalparticles]];
WriteRealCorrections[NRoutine,DecayingParticle,twobodydecays,OneLtbds,type];
SA`SavedInformationOneLoopDecaysReal=Join[SA`SavedInformationOneLoopDecaysReal,{{DecayingParticle,{OneLtbds[[1]],OneLtbds[[2]]},NRoutine,OneLtbds}}];];

WriteRealCorrections[NameRoutine_,DecayingParticle_,twobodydecays_,OneLtbds_,vertextype_]:=Block[{i,j,k,couplings,masses,p1,p2,ps,ind,coup,dim,dimI,treecolfac,Casarr,Chargearr,Chargevec,chargepresent,Casvec,chargeeasy,colourpresent,pCas,pDyn,routine,type},

couplings=getSPhenoCoupling[OneLtbds[[3]]][[1]];
masses=DeleteCases[Intersection[SPhenoMass/@List@@OneLtbds[[3]]],0.];



MakeSubroutineTitle["Gamma_Real_"<>ModelName<>"_"<>NameRoutine,Flatten[{couplings,masses}],{"MLambda","em","gs"},{"GammarealPhoton","GammarealGluon"},outputfortran];
WriteString[outputfortran,"Implicit None\n\n"];
MakeVariableList[couplings,", Intent(in)",outputfortran];
MakeVariableList[masses,", Intent(in)",outputfortran];
p1=OneLtbds[[1]];
p2=OneLtbds[[2]];
ps={DecayingParticle,p1,p2};

(*Define colour properties*)
Casvec=Array[0&,{3}];
Chargevec=Casvec;
chargepresent=False;
chargeeasy=False;

Chargevec=getElectricCharge/@{DecayingParticle,p1,p2};
If[Chargevec=!={0,0,0},chargepresent=True,chargeeasy=True];

If[chargepresent,Chargearr=Flatten[Table[OneLtbds[[4]]*Chargevec[[i]]*Chargevec[[j]],{i,1,3},{j,i,3}]];];

dim="("<>ToString[getGen[DecayingParticle]]<>","<>ToString[getGen[p1]]<>","<>ToString[getGen[p2]]<>")";
dim=Nest[StringReplace[#,{"(1,"->"(",",1,"->",",",1)"->")","(1)"->""}]&,dim,3];
WriteString[outputfortran,"Real(dp), Intent(in) :: MLambda, em, gs \n\n"];
WriteString[outputfortran,"Real(dp), Intent(out) :: GammarealPhoton"<>dim<>", GammarealGluon"<>dim<>" \n\n"];
WriteString[outputfortran,"Integer :: i1, i2, i3, gt1, gt2, gt3, kont \n"];
WriteString[outputfortran,"Real(dp) :: Mex1, Mex2, Mex3, ExtRMsq \n\n"];
WriteString[outputfortran,"Real(dp) :: Mloop1, Mloop2, Mloop3 \n"];

If[vertextype=!=FFS||getType[DecayingParticle]===F,
type=vertextype,
(* In the case that we have a scalar decaying to two fermions need to change the type round*)
type=SFF];


If[type=!=SSS&&type=!=SVV&&type=!=SSV,(*fermionic processes*)
WriteString[outputfortran,"Complex(dp) :: CoupL, CoupR \n \n"];,
(*else*)
WriteString[outputfortran,"Complex(dp) :: Coup \n \n"];
];

ind=MakeIndicesCoupling[{AntiField[DecayingParticle],i1},{p1,i2},{p2,i3},getSPhenoCoupling[OneLtbds[[3]],SPhenoCouplingsAll][[2]]];

coup=getSPhenoCoupling[OneLtbds[[3]],SPhenoCouplingsAll];

dimI="(i1,i2,i3)";
If[getGen[DecayingParticle]>1,WriteString[outputfortran,"Do i1="<>ToString[getGenSPhenoStart[DecayingParticle]]<>","<>ToString[getGenSPheno[DecayingParticle]]<>"\n"];,dimI=StringReplace[dimI,"i1,"->""];];
If[getGen[p1]>1,WriteString[outputfortran,"  Do i2="<>ToString[getGenSPhenoStart[p1]]<>","<>ToString[getGenSPheno[p1]]<>"\n"];,dimI=StringReplace[dimI,"i2,"->""];];
If[getGen[p2]>1,WriteString[outputfortran,"    Do i3="<>ToString[getGenSPhenoStart[p2]]<>","<>ToString[getGenSPheno[p2]]<>"\n"];,dimI=StringReplace[dimI,"i3)"->")"];];
dimI=StringReplace[dimI,{"()"->"",",)"->")"}];

If[type=!=SSS&&type=!=SVV&&type=!=SSV,(*fermionic processes*)
WriteString[outputfortran,"CoupL = "<>ToString[coup[[1,1]]]<>ind[[1]]<>"\n"];
WriteString[outputfortran,"CoupR = "<>ToString[coup[[1,2]]]<>ind[[1]]<>"\n"];
,(* Scalar processes *)
WriteString[outputfortran,"Coup = "<>ToString[coup[[1,1]]]<>ind[[1]]<>"\n"];
];

WriteString[outputfortran,"Mex1 = "<>SPhenoMass[DecayingParticle,i1]<>"\n"];

WriteString[outputfortran,"Mex2 = "<>SPhenoMass[p1,i2]<>"\n"];
WriteString[outputfortran,"Mex3 = "<>SPhenoMass[p2,i3]<>"\n"];
numberingaugegroupofcolor=(Position[Gauge,color])[[1,1]];
colourpresent=False;
treecolfac=0;
If[type=!=FFV&&type=!=SSV&&type=!=SVV,
{treecolfac,Casarr}=getbremsgeneral[ps[[1]],ps[[2]],ps[[3]],numberingaugegroupofcolor];
If[treecolfac!=0,colourpresent=True];
,
pDyn=SA`DynL[DecayingParticle,numberingaugegroupofcolor];
If[pDyn=={0},
colourpresent=False;,
colourpresent=True;
pCas=SA`Casimir[DecayingParticle,numberingaugegroupofcolor];
];
];


WriteString[outputfortran,"If (Abs(Mex1).gt.(Abs(Mex2)+Abs(Mex3))) Then \n"];
Switch[type,
FFV,

	If[chargepresent,
	If[getElectricCharge[p2]===0,routine="FFZ";,routine="FFW"]; 

	
	WriteString[outputfortran,"  Call hardphoton"<>routine<>"(Mex1,Mex2,Mex3,MLambda,"<>If[Head[getElectricCharge[DecayingParticle]]===Integer,SPhenoForm[1. getElectricCharge[DecayingParticle]],SPhenoForm[getElectricCharge[DecayingParticle]]]<>","<>If[Head[getElectricCharge[p1]]===Integer,SPhenoForm[1. getElectricCharge[p1]],SPhenoForm[getElectricCharge[p1]]]<>",CoupL,CoupR,(0,1)*em,GammaRealPhoton"<>dimI<>",kont)\n"];
	,
	WriteString[outputfortran,"  GammarealPhoton"<>dimI<>" = 0._dp \n\n"];
	];
		

		(*pCas=SA`Casimir[RE[DecayingParticle],numberingaugegroupofcolor];*)
	(*pDyn=SA`DynL[DecayingParticle,3];*)

		If[colourpresent==False,
		WriteString[outputfortran,"  GammarealGluon"<>dimI<>" = 0._dp \n"];,
		WriteString[outputfortran,"  Call hardgluonFFZW(Mex1,Mex2,Mex3,MLambda,"<>ToString[FortranForm[If[Head[pCas]===Integer,1.0*pCas,pCas]]]<>",gs,CoupL,CoupR,Gammarealgluon"<>dimI<>",kont)\n"];
		];(* if pDyn\[Equal]0*)

	,(*after FFV*)
FFS,
(*Photon..*)

If[!chargepresent,
WriteString[outputfortran," Gammarealphoton"<>dimI<>" = 0._dp \n"];,WriteString[outputfortran," Call hardradiationFFS(Mex1,Mex2,Mex3,MLambda,em,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Chargearr,","]]<>",CoupL,CoupR,Gammarealphoton"<>dimI<>",kont)\n"];
];
(*gluon...*)

If[!colourpresent,
WriteString[outputfortran,"  GammarealGluon"<>dimI<>" = 0._dp \n"];
,
(*
WriteString[outputfortran," Call hardradiationFFS(Mex1,Mex2,Mex3,MLambda,gs,"<>StringRiffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Casarr,","]<>",CoupL,CoupR,Gammarealgluon"<>dimI<>",kont)\n"];
	];
*)
WriteString[outputfortran," Call hardradiationFFS(Mex1,Mex2,Mex3,MLambda,gs,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Casarr,","]]<>",CoupL,CoupR,Gammarealgluon"<>dimI<>",kont)\n"];
	];


,(* end FFS*)
SFF,

If[!chargepresent,
WriteString[outputfortran," Gammarealphoton"<>dimI<>" = 0._dp \n"];,WriteString[outputfortran," Call hardradiationSFF(Mex1,Mex2,Mex3,MLambda,em,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Chargearr,","]]<>",CoupL,CoupR,Gammarealphoton"<>dimI<>",kont)\n"];
];
If[!colourpresent,
WriteString[outputfortran,"  GammarealGluon"<>dimI<>" = 0._dp \n"];
,
WriteString[outputfortran," Call hardradiationSFF(Mex1,Mex2,Mex3,MLambda,gs,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Casarr,","]]<>",CoupL,CoupR,Gammarealgluon"<>dimI<>",kont)\n"];
	];
, (* end SFF *)
SSS,

If[!chargepresent,
WriteString[outputfortran," Gammarealphoton"<>dimI<>" = 0._dp \n"];,WriteString[outputfortran," Call hardradiationSSS(Mex1,Mex2,Mex3,MLambda,em,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Chargearr,","]]<>",Coup,Gammarealphoton"<>dimI<>",kont)\n"];
];
If[!colourpresent,
WriteString[outputfortran,"  GammarealGluon"<>dimI<>" = 0._dp \n"];
,
WriteString[outputfortran," Call hardradiationSSS(Mex1,Mex2,Mex3,MLambda,gs,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Casarr,","]]<>",Coup,Gammarealgluon"<>dimI<>",kont)\n"];
	];
,(* end SSS *)

SSV,
If[!chargepresent,
WriteString[outputfortran," Gammarealphoton"<>dimI<>" = 0._dp \n"];,
WriteString[outputfortran," Call hardradiationSSV(Mex1,Mex2,Mex3,MLambda,em,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Chargevec,","]]<>",Coup,Gammarealphoton"<>dimI<>",kont)\n"];
];

If[!colourpresent,
WriteString[outputfortran," GammarealGluon"<>dimI<>" = 0._dp \n"];,
WriteString[outputfortran," Call hardradiationSSV(Mex1,Mex2,Mex3,MLambda,gs,1._dp,1._dp,0._dp,Coup,Gammarealgluon"<>dimI<>",kont)\n"];
WriteString[outputfortran," GammarealGluon"<>dimI<>" = GammarealGluon"<>dimI<> "*"<>ToString[FortranForm[pCas]]<>" \n"];
];
, (* End SSV *)

SVV,
WriteString[outputfortran,"  GammarealGluon"<>dimI<>" = 0._dp \n"];
If[!chargepresent,
WriteString[outputfortran," Gammarealphoton"<>dimI<>" = 0._dp \n"];,
WriteString[outputfortran," Call hardphotonSVV(Mex1,Mex2,Mex3,MLambda,em,"<>StringJoin[Riffle[ToString[FortranForm[If[Head[#]===Integer,1.0*#,#]]]&/@Chargevec,","]]<>",Coup,Gammarealphoton"<>dimI<>",kont)\n"];
];


];(* end Switch *)

WriteString[outputfortran,"Else \n"];
WriteString[outputfortran,"  GammarealGluon"<>dimI<>" = 0._dp \n"];
WriteString[outputfortran,"  GammarealPhoton"<>dimI<>" = 0._dp \n\n"];
WriteString[outputfortran,"End if \n"];

If[getGen[p2]>1,WriteString[outputfortran,"    End Do\n"];];
If[getGen[p1]>1,WriteString[outputfortran,"  End Do\n"];];
If[getGen[DecayingParticle]>1,WriteString[outputfortran,"End Do\n"];];

WriteString[outputfortran,"End Subroutine Gamma_Real_"<>ModelName<>"_"<>NameRoutine<>"\n\n\n"];

];







