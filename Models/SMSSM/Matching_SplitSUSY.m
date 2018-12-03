$NameUV="SingletSplitSUSY";

$ParametrisationUV = {
   vd -> v Cos[ArcTan[TanBeta]],
   vu -> v Sin[ArcTan[TanBeta]],
   v -> epsUV,
   Yu[a_,b_] :> Delta[3,b] Delta[a,3] Yu[a,b],
   Yd[a_,b_] :> Delta[3,b] Delta[a,3] Yd[a,b],
   Ye[a_,b_] :> Delta[3,b] Delta[a,3] Ye[a,b],
   T[Yu][a__] :> epsUV,
   T[Yd][a__] :> epsUV,
   T[Ye][a__] :> epsUV,
   mq2[a__] :> Delta[a] MSUSY^2,
   mu2[a__] :> Delta[a] MSUSY^2,
   md2[a__] :> Delta[a] MSUSY^2,
   me2[a__] :> Delta[a] MSUSY^2,
   ml2[a__] :> Delta[a] MSUSY^2,
   conj[x_] -> x,
   MassB->epsUV,
   MassWB->epsUV,
   MassG->epsUV,
   vS->epsUV,
   B[MS]->epsUV^2,
   MS->epsUV,
   B[\[Mu]]->MSUSY^2 TanBeta/(1+TanBeta^2),
   L[L1]->epsUV^3,
   L1->epsUV^2
   };
   
$SimplificationsMatching={
  conj[x_] -> x,
  UVscaleQ->MSUSY
};
   
$AssumptionsMatching={
 TanBeta>0, 
 MSUSY>0
}


$SolveTadpolesUV = {mHd2, mHu2, ms2};

$MatchingConditions = {
   LambdaH -> -1/3 hh[indD].hh[indD].hh[indD].hh[indD],
   LambdaS -> -1/3 hh[indS].hh[indS].hh[indS].hh[indS],
   LambdaSH -> -hh[indD].hh[indD].hh[indS].hh[indS],
   KappaS -> -1/(3 Sqrt[2]) hh[indS].hh[indS].hh[indS],
   KappaSH -> -1/Sqrt[2] hh[indS].hh[indD].hh[indD]
   };

ExtractExternalIndices[input_]:=Block[{},
If[(ZH[1,3]/.input)==1,
 indS=1;
 indD=2;,
 indS=2;,
 indH=1;
];
];

