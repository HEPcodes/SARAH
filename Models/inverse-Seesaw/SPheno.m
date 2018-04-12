MINPAR={{1,m0},
        {2,m12},
        {3,TanBeta},
        {4,SignumMu},
        {5,Azero}};

RealParameters = {TanBeta,m0};

ParametersToSolveTadpoles = {\[Mu],B[\[Mu]]};

RenormalizationScaleFirstGuess = m0^2 + 4 m12^2;
RenormalizationScale = MSu[1]*MSu[6];

ConditionGUTscale = g1 == g2;

BoundaryHighScale={
{g1, Sqrt[(g1^2 + g2^2)/2]},
{g2, g1},
{T[Ye], Azero*Ye},
{T[Yd], Azero*Yd},
{T[Yu], Azero*Yu},
{T[Yv], Azero*Yv},
{mq2, DIAGONAL m0^2},
{ml2, DIAGONAL m0^2},
{md2, DIAGONAL m0^2},
{mu2, DIAGONAL m0^2},
{me2, DIAGONAL m0^2},
{mv2, DIAGONAL m0^2},
{ms2, DIAGONAL m0^2},
{mHd2, m0^2},
{mHu2, m0^2},
{MassB, m12},
{MassWB,m12},
{MassG,m12}
};


BoundarySUSYScale={
{MUS, LHInput[MUS]},
{Mv, LHInput[Mv]},
{Yv,LHInput[Yv]},
{B[MUS], LHInput[B[MUS]]},
{B[Mv], LHInput[B[Mv]]}
};


InitializationValues = {
 {\[Mu],0},
 {B[\[Mu]],0}

};

QuadruplePrecision = {Fv};



BoundaryLowScaleInput={
 {vd,Sqrt[4 mz2/(g1^2+g2^2)]*Cos[ArcTan[TanBeta]]},
 {vu,Sqrt[4 mz2/(g1^2+g2^2)]*Sin[ArcTan[TanBeta]]}
};

DEFINITION[MatchingConditions]=Default[THDMII];

UseHiggs2LoopMSSM = True;

ListDecayParticles = Automatic;
ListDecayParticles3B = Automatic;

(* Example for input values *)
(* Not including neutrino fit! *)
DefaultInputValues = {m0 -> 1000, m12 -> 1500, TanBeta->20, SignumMu ->1, Azero -> -1500, Yv[a_,a_]->0.1, Mv[a_,a_]->1000, MUS[a_,a_]->0.001};
