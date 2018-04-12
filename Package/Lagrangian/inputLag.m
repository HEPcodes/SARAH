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



CreateLagrangian[terms_,HC_,Overwrite_]:=Block[{i,j,pos,temp,fields,list={},coup,numerical,res,lVVV,lVVVV,rVVV,rVVVV,rPot,rKin},
DynamicStatusAddTerms[terms]="initializing";
If[SupersymmetricModel===False,InitChargeFactors;];
DynamicStatusAddTerms[terms]="expand terms";
If[Head[Expand[terms]]===Plus,
TermList=(CreateTermList/@(List@@Expand[terms]));,
TermList={CreateTermList[Expand[terms]]};
];
LagInputIncludeHC = HC;

lVVV = Select[TermList,(getPartCode[#[[3]]]==30)&];
lVVVV = Select[TermList,(getPartCode[#[[3]]]==40)&];
pCodeKin={20,2010,210,220,120};
pCodePot={2000,100,200,300,400,2100};
lKin = Select[TermList,(FreeQ[pCodeKin,getPartCode[#[[3]]]]==False)&];
lPot = Select[TermList,(FreeQ[pCodePot,getPartCode[#[[3]]]]==False)&];
DynamicStatusAddTerms[terms]="checking charge conservation";
rVVV=Plus@@( MakeTerms/@lVVV);
rVVVV=Plus@@(MakeTerms/@lVVVV);
rPot=Plus@@(MakeTerms/@lPot);
rKin=Plus@@(MakeTerms/@lKin);
DynamicStatusAddTerms[terms]="All Done";
If[Overwrite===True, 
Return[{rVVV+rVVVV+rPot+rKin,rVVV,rVVVV,rPot,rKin}];,
Return[{0,rVVV,rVVVV,rPot,rKin}];
];


];


(*
CreateTermList[terms_]:=Block[{i,j,pos,temp,fields,list={},coup,numerical,potenz,structure},
If[Head[terms]===Plus,
For[i=1,i\[LessEqual]Length[terms],
fields = List@@Select[List@@terms[[i]],(Head[#]\[Equal]Dot)&][[1]];
(* numerical = Cases[List@@terms[[i]],x_?NumberQ,2]; *)
numerical = Cases[List@@terms[[i]],x_?NumberQ,1]; (* CHECK *)
structure=Flatten[{Cases[List@@terms[[i]],x_Delta],Cases[List@@terms[[i]],x_epsTensor],Cases[List@@terms[[i]],x_CG]}];
If[numerical==={},numerical=1;,numerical=numerical[[1]]];
If[structure=!={},numerical=Times@@structure numerical];
coupling=List@@Select[List@@terms[[i]],(Head[#]=!=Dot && Head[#]\[Equal]Symbol)&];
Pcoupling=List@@Select[List@@terms[[i]],(NumericQ[#]===False&&Head[#]===Power)&,99];
If[Pcoupling==={},potenz=1;,potenz=Pcoupling[[1,2]]; coupling=Join[coupling,{Pcoupling[[1,1]]}];];
If[coupling==={}, coupling=NoName;,coupling=coupling[[1]];];
list=Join[list,{{numerical,coupling,fields,potenz}}];
i++;];,
fields = List@@Select[List@@terms,(Head[#]\[Equal]Dot)&][[1]];
(* numerical = Cases[List@@terms,x_?NumberQ,2]; *)
numerical = Cases[List@@terms,x_?NumberQ,1]; (* CHECK *)
structure=Flatten[{Cases[List@@terms,x_Delta],Cases[List@@terms,x_epsTensor],Cases[List@@terms,x_CG]}];
If[numerical==={},numerical=1;,numerical=numerical[[1]]];
If[structure=!={},numerical=Times@@structure numerical];
coupling=List@@Select[List@@terms,(Head[#]=!=Dot && Head[#]\[Equal]Symbol && Head[#]=!=Delta && Head[#]=!=epsTensor && Head[Head[#]]=!=CG)&];
Pcoupling=List@@Select[List@@terms,(NumericQ[#]===False&&Head[#]===Power)&,99];
If[Pcoupling==={},potenz=1;,potenz=Pcoupling[[1,2]]; coupling=Join[coupling,{Pcoupling[[1,1]]}];];
If[coupling==={}, coupling=NoName;,coupling=coupling[[1]];];
list=Join[list,{{numerical,coupling,fields,potenz}}];
];
Return[list];
]; *)

CreateTermList[terms_]:=Block[{i,j,pos,temp,fields,list={},coup,coeff,potenz,structure,coupling},
If[FreeQ[FullForm[terms],Dot]==False,
fields=Select[List@@terms,(Head[#]==Dot)&][[1]];,
fields=Select[List@@terms,(FreeQ[SFields,#]==False || FreeQ[FFields,#]==False )&];
];

coupling=List@@Select[List@@terms,(Head[#]=!=Dot&&(Head[#]==Symbol || MemberQ[{B,T,L},Head[#]])&&Head[#]=!=Delta&&Head[#]=!=epsTensor&&Head[#]=!=Sig && Head[#]=!=Lam &&Head[Head[#]]=!=CG  && FreeQ[SFields,#]  && FreeQ[FFields,#]) &];
coupling=Select[coupling,FreeQ[Transpose[parameters][[1]],#]&];
If[coupling==={},
LagInput::NoCoupling="No coupling defined for term in Lagrangian involving ``";
Message[LagInput::NoCoupling,terms];,
coupling=coupling[[1]];
];
coeff=Simplify[terms/coupling/fields];
fields=List@@fields;
potenz=Select[Cases[List@@terms,x_Power,3],(FreeQ[#,coupling]==False)&]/.Power[a_,b_]->b/.{a_Integer}->a/.{}->1;
list={coeff,coupling,fields,potenz};
Return[list];
];


MakeTerms[entry_]:=Block[{i,i1,j,temp,part,particles,particleNr=1,fermNr=1,head,newParticle,pos,invFields={},coup,fields,withHead={},head2, headDer,lorIndex,Fincluded=False,IndStructure},
fields=entry[[3]];
(* coup = genTest[entry[[2]],fields,False]^entry[[4]]; *)
AdditionalParametersLagrange = Join[AdditionalParametersLagrange,{entry[[2]]}];
particles=1;

For[i=1,i<=Length[fields],
If[Head[fields[[i]]]===Der,
headDer=DER;
newParticle=fields[[i,1]];
lorIndex=fields[[i,2]];,
newParticle=fields;
headDer=Evaluate;
newParticle=fields[[i]];
];
If[Head[newParticle]===conj, head=conj;,head=Evaluate;];
newParticle = RE[newParticle];

If[FreeQ[SFields,newParticle]==False,
pos = Position[SFields,newParticle][[1,1]];
head2=getHead[SFields[[pos]]];
invFields = Join[invFields,{Fields[[pos,3]]}];
withHead = Join[withHead,{head2[head[Fields[[pos,3]]]]}];
];

If[FreeQ[FFields,newParticle]==False,
pos = Position[FFields,newParticle][[1,1]];
head2=getHead[FFields[[pos]]];
invFields = Join[invFields,{Fields[[pos,3]]}];
withHead = Join[withHead,{head2[head[Fields[[pos,3]]]]}];
];

If[getTypeOld[newParticle]===F,
Fincluded=True;
particles*=((head[RE[getFull2[newParticle][fermNr]/.subGC[particleNr]]]) sum[genf[particleNr],getGenStart[newParticle],GetGen[newParticle]] /. DER[a_]->Der[a,lorIndex]);
particleNr++;
fermNr++;,
particles*=((headDer[head[RE[getFull2[newParticle]/.subGC[particleNr]]]]) sum[genf[particleNr],getGenStart[newParticle],GetGen[newParticle]]/. DER[a_]->Der[a,lorIndex]);
particleNr++;
];
i++;];

PrintDebug["Checking charge Conservation of ", Dot@@entry[[3]]];

SA`LagTermsNonSUSY=Join[SA`LagTermsNonSUSY,{withHead}];

CheckLagTermChargeConservation[withHead];

For[i1=1,i1<=Length[Global],
SA`CheckGlobalLagLevel=True;
If[CheckChargeConservationGlobal[List@@withHead,Global[[i1,2]]]==False,
Message[Lagrange::ViolationGlobal ,Global[[i1,2]]];
];
SA`CheckGlobalLagLevel=False;
i1++;];


If[FreeQ[entry[[1]],Delta] && FreeQ[entry[[1]],epsTensor] && FreeQ[entry[[1]],CG],
temp=MakeIndexStructure[withHead] /. subFieldsOne;,
temp=entry[[1]]/. subFieldsOne;
];
coup = genTest[entry[[2]],fields,False]^entry[[4]];

If[SupersymmetricModel===False,
GenerateCGCsForBrokenGroups[temp,fields,particles,invFields,withHead]; 
If[Head[ContractionRGE[entry[[2]]]]===ContractionRGE,
IndStructure=MakeIndexStructureRGEnonSUSY[withHead] /.subCGCBroken;,
IndStructure=ContractionRGE[coup]
];
ExtractSymmetryOfParametersNS[entry[[2]],withHead,IndStructure];
(*
IndStructure=temp*GetNormalizationFactors[invFields];
*)

(* adapted sign *)
Switch[Plus@@(getPartCode[getBlank[#]&/@fields]),
2100,
SA`FFSlist=Join[SA`FFSlist,{{List@@(particles /. sum[a__]->1 /.A_[b_Integer]->A),-coup  IndStructure entry[[1]],False}}];,
2000,
SA`FFlist=Join[SA`FFlist,{{List@@(particles /. sum[a__]->1/.A_[b_Integer]->A),-coup  IndStructure entry[[1]],False}}];,
200,
SA`SSlist=Join[SA`SSlist,{{List@@(particles /. sum[a__]->1),-coup  IndStructure (entry[[1]]/. {Delta[a__]->1,epsTensor[__]->1}),LagInputIncludeHC}}];,
300,
SA`SSSlist=Join[SA`SSSlist,{{List@@(particles /. sum[a__]->1),-coup  IndStructure (entry[[1]]/. {Delta[a__]->1,epsTensor[__]->1}),LagInputIncludeHC}}];,
400,
SA`SSSSlist=Join[SA`SSSSlist,{{List@@(particles /. sum[a__]->1),-coup  IndStructure (entry[[1]]/. {Delta[a__]->1,epsTensor[__]->1}),LagInputIncludeHC}}];
];
];

If[FreeQ[entry[[1]],Delta] && FreeQ[entry[[1]],epsTensor] && FreeQ[entry[[1]],CG],
temp=entry[[1]]*coup*particles*MakeIndexStructure[withHead] /. subFieldsOne;,
temp=entry[[1]]*coup*particles /. subFieldsOne;
];

SA`LagrangianContractions=Join[SA`LagrangianContractions,{{Dot@@withHead,temp /. sum[a__]->1}}];
If[invFields=!={},
temp = SumOverExpandedIndizes[temp,invFields];
];


If[Fincluded==True,
temp = temp + (temp /. { A_?(#=!=SU&)[1]->A[2],A_?(#=!=SU&)[2]->A[1]});
];


If[LagInputIncludeHC,
temp = temp + (conj[temp]  /. {A_?(#=!=SU&)[1]->A[2],A_?(#=!=SU&)[2]->A[1],A_?(#=!=SU&)[3]->A[4],A_?(#=!=SU&)[4]->A[3]}  /. CG[a_,b_]:>CG[a,Reverse/@b])
];



Return[temp /. conj[gamma[a_]]->gamma[a]];



];

getHead[x_]:=If[Head[x]===conj, Return[conj];,Return[Evaluate];];

DC[fieldIN_,ind1_,ind2_,lor_]:=Block[{},
If[Head[fieldIN]===conj,
head=conj;field=RE[fieldIN];,
head=Evaluate; field=fieldIN;
];
pos=Position[SFields,field];
If[pos==={},
pos=Position[FFields,field];
];
pos=pos[[1,1]];
Return[(Deri[head[part[getFull2[field],2]],lorentz/.subGC[lor]]+ head[part[getFull2[field],2]]KovariantDerivative[pos,ind1,ind2,lor])];
];

FS[i_Integer,i1_,i2_]:= (Deri[part[SGauge[[i]],i1],lorentz /. subGC[i2]] -Deri[part[SGauge[[i]],i2],lorentz /.subGC[i1]]);
FS[i_Symbol,i1_,i2_]:= (Deri[part[SGauge[[Position[SGauge,i][[1,1]]]],i1],lorentz /. subGC[i2]] -Deri[part[SGauge[[Position[SGauge,i][[1,1]]]],i2],lorentz /.subGC[i1]]);

del[x_Times,y_]:=del[x[[1]],y] Times@@ Take[List@@x,{2,-1}]+x[[1]] del[ Times@@ Take[List@@x,{2,-1}],y];
del[x_Dot,y_]:=del[x[[1]],y].Dot@@ Take[List@@x,{2,-1}]+x[[1]].del[Dot@@ Take[List@@x,{2,-1}],y];
del[x_sum,y_]:=0;
del[x_Rational,y_]:=0;
del[x_Integer,y_]:=0;
del[x_Complex,y_]:=0;
del[x_Sin,y_]:=0;
del[x_Cos,y_]:=0;
del[x_Tan,y_]:=0;
del[x_ArcTan,y_]:=0;
del[x_ArcSin,y_]:=0;
del[x_ArcCos,y_]:=0;
del[x_,y_]:=Deri[x,y] /; (FreeQ[SFields,RE[x]]===False || FreeQ[FFields,RE[y]]===False);


