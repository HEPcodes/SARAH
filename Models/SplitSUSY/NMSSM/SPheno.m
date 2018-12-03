MINPAR={
        {1,  BinoNMSSM},
        {2,  WinoNMSSM},
        {3,  GluinoNMSSM},
        {4,  TanBeta},
        {5,  LambdaNMSSM},
        {6,  KappaNMSSM},
        {7,  TLambdaNMSSM},
        {8,  TKappaNMSSM},
        {9,  vSNMSSM},
        {10, m0},
        {11, Azero},
        {12, MAinput},
        {13, Muinput}
};

RenormalizationScaleFirstGuess = 100^2;
RenormalizationScale = v^2;

RealParameters = {TanBeta,m0,MAinput,vSNMSSM}; 


UseParameterAsGUTscale = {m0};

BoundaryHighScale={
 {\[Lambda], LambdaNMSSM},
 {\[Kappa], KappaNMSSM}, 
 {T[\[Lambda]], TLambdaNMSSM},
 {T[\[Kappa]], TKappaNMSSM}, 
 {T[Ye], Azero*Ye},
 {T[Yd], Azero*Yd},
 {T[Yu], Azero*Yu},
 {mq2, DIAGONAL m0^2},
 {ml2, DIAGONAL m0^2},
 {md2, DIAGONAL m0^2},
 {mu2, DIAGONAL m0^2},
 {me2, DIAGONAL m0^2},
 {MassB,  BinoNMSSM},
 {MassWB, WinoNMSSM}, 
 {MassG,  GluinoNMSSM}
};

BoundaryRenScale={
  {vS,        vSNMSSM}
}


ParametersToSolveTadpoles = {m2,mSr2};

DEFINITION[MatchingConditions]=Default[OHDM];

ListDecayParticles = Automatic;
ListDecayParticles3B = Automatic;

DefaultInputValues ={m0 -> 10^6, m12 -> 10^6, Azero -> 0, Muinput->10^6, MAinput->10^6, TanBeta->2};


(*---------------------------------------------------*)
(* information for matching to MSSM *)
(*---------------------------------------------------*)

MatchingToModel= {"SMSSM"};
MatchingScale = {m0};

IncludeParticlesInThresholds={
 {hh,Ah,Hpm,Su,Sd,Se,Sv}
};


  
AssumptionsMatchingScale={
 {
  {vd,epsUV*Cos[ArcTan[TanBeta]]}, 
  {vu,epsUV*Sin[ArcTan[TanBeta]]},
  {vS,epsUV},
  {B[MS], epsUV^2},
  {MS, epsUV},
  {L[L1], epsUV^3 },
  {L1, epsUV^2},
  {\[Mu], epsUV^2}
 }
};

(*
BoundaryMatchingScaleUp={
 {
 {Yu, Sqrt[1+TanBeta^2]/TanBeta*Yu},
 {Yd, Sqrt[1+TanBeta^2]*Yd},
 {Ye, Sqrt[1+TanBeta^2]*Ye},
 {vS,        vSNMSSM},
 {B[\[Mu]], MAinput^2/(TanBeta + 1/TanBeta)},
 {\[Mu],Muinput}
}
};

BoundaryMatchingScaleDown={
 {
 {Yu, TanBeta/Sqrt[1+TanBeta^2]*Yu},
 {Yd, 1/Sqrt[1+TanBeta^2]*Yd},
 {Ye, 1/Sqrt[1+TanBeta^2]*Ye},
 {MG,        MassB },
 {MW,        MassWB},
 {MB,        MassG},
 {vS,        vSNMSSM},
 {Ysr,       \[Kappa]},
 {Ysi,       \[Kappa]},
 {Ysu,       \[Lambda] Sin[ArcTan[TanBeta]]},
 {Ysd,       \[Lambda] Cos[ArcTan[TanBeta]]},
 {mSi2,      mSr2},
 {g1u,Sqrt[3/5] g1 Sin[ArcTan[TanBeta]]   },
 {g2u,g2 Sin[ArcTan[TanBeta]]   }, 
 {g1d,Sqrt[3/5] g1 Cos[ArcTan[TanBeta]]   },
 {g2d,g2 Cos[ArcTan[TanBeta]]   } ,
 {LambdaH , -1/3 EFTcoupNLO[hh[indD].hh[indD].hh[indD].hh[indD]]},
 {LambdaSr , -1/3 EFTcoupNLO[hh[indS].hh[indS].hh[indS].hh[indS]]},
 {LambdaSi, LambdaSr},
 {LambdaSri, LambdaSr},
 {LambdaSHr, -hh[indD].EFTcoupNLO[hh[indD].hh[indS].hh[indS]]},
 {LambdaSHi, LambdaSHr},
 {KappaSr  , -1/(3 Sqrt[2]) EFTcoupNLO[hh[indS].hh[indS].hh[indS]]},
 {KappaSri, KappaSr},
 {KappaSi, KappaSr},
 {KappaSHr , -1/Sqrt[2] EFTcoupNLO[hh[indS].hh[indD].hh[indD]]    }
 }
};


indS=1;
indD=2;
*)

BoundaryMatchingScaleUp={
 {
 {Yu, Sqrt[1+TanBeta^2]/TanBeta*Yu},
 {Yd, Sqrt[1+TanBeta^2]*Yd},
 {Ye, Sqrt[1+TanBeta^2]*Ye},
 {B[\[Mu]], MAinput^2/(TanBeta + 1/TanBeta)},
 {\[Mu],Muinput}
}
};

BoundaryMatchingScaleDown={
 {
 {MG,        MassB },
 {MW,        MassWB},
 {MB,        MassG}, 
 {mSi2,      mSr2},
 
 (* new Yukawa couplings *)
 {Ysr,  -1/Sqrt[2]*(EFTcoupNLO[Chi[FsR].Chi[FsR].hh[phiS]][PL])},
 {Ysi,       Ysr},
 {Ysu,  1/(2 Sqrt[2])*(EFTcoupNLO[Chi[FHu0].Chi[FsR].hh[phiU]][PL]/(ZN[FHu0,4] ZN[FsR,5]))},
 {Ysd,  1/(2 Sqrt[2])*(EFTcoupNLO[Chi[FHd0].Chi[FsR].hh[phiU]][PL]/(ZN[FHd0,3] ZN[FsR,5]))}, 
 {g1u, -2*(EFTcoupNLO[Chi[fB].Chi[FHu0].hh[phiu]][PL]/(ZN[fB,1] ZN[FHu0,4]))},
 {g2u,  2*(EFTcoupNLO[Chi[fW0].Chi[FHu0].hh[phiu]][PL]/(ZN[fW0,2] ZN[FHu0,4]))}, 
 {g1d,  2*(EFTcoupNLO[Chi[fB].Chi[FHd0].hh[phiu]][PL]/(ZN[fB,1] ZN[FHd0,3]))},
 {g2d, -2*(EFTcoupNLO[Chi[fW0].Chi[FHd0].hh[phiu]][PL]/(ZN[fW0,2] ZN[FHd0,3]))},
 
  (* new scalar couplings *)
 {LambdaH,  -1/3 EFTcoupNLO[hh[phiu].hh[phiu].hh[phiu].hh[phiu]]},
 {LambdaSr,  -1/3 EFTcoupNLO[hh[phiS].hh[phiS].hh[phiS].hh[phiS]]},
 {LambdaSi,  LambdaSr},
 {LambdaSri, LambdaSr},
 {LambdaSHr, -hh[phiu].EFTcoupNLO[hh[phiu].hh[phiS].hh[phiS]]},
 {LambdaSHi, LambdaSHr},
 {KappaSr,   -1/(3 Sqrt[2]) EFTcoupNLO[hh[phiS].hh[phiS].hh[phiS]]},
 {KappaSri,  KappaSr},
 {KappaSHr,  -1/Sqrt[2] EFTcoupNLO[hh[phiS].hh[phiu].hh[phiu]]    }
 }
};



ParametersToSolveTadpoleMatchingScale={
 {mHd2,mHu2,ms2}
};
 


