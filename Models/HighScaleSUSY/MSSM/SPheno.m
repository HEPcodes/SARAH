MINPAR={{ 1, m0},
        { 2, m12},
        { 3, TanBeta},
        { 5, Azero},
        {23, Muinput},
        {24, MAinput}     
        };
        
RealParameters = {TanBeta,m0,MAinput};        

RenormalizationScaleFirstGuess = 100^2;
RenormalizationScale = v^2;


UseParameterAsGUTscale = {m0};

BoundaryHighScale={
 {T[Ye], Azero*Ye},
 {T[Yd], Azero*Yd},
 {T[Yu], Azero*Yu},
 {mq2, DIAGONAL m0^2},
 {ml2, DIAGONAL m0^2},
 {md2, DIAGONAL m0^2},
 {mu2, DIAGONAL m0^2},
 {me2, DIAGONAL m0^2},
 {MassB, m12},
 {MassWB,m12},
 {MassG,m12},
 {\[Mu], Muinput},
 {B[\[Mu]], MAinput^2/(TanBeta + 1/TanBeta)}
};


ParametersToSolveTadpoles = {m2};

DEFINITION[MatchingConditions]=Default[OHDM];

ListDecayParticles = Automatic;
ListDecayParticles3B = Automatic;

DefaultInputValues ={m0 -> 10^6, m12 -> 10^6, Azero -> 0, Muinput->10^6, MAinput->10^6, TanBeta->2};


(*---------------------------------------------------*)
(* information for matching to MSSM *)
(*---------------------------------------------------*)

MatchingToModel= {"MSSM"};
MatchingScale = {m0};

IncludeParticlesInThresholds={
 {hh,Ah,Hpm,Su,Sd,Se,Sv,Chi,Cha}
};

AssumptionsMatchingScale={
 {
  {vd,epsUV*Cos[ArcTan[TanBeta]]}, 
  {vu,epsUV*Sin[ArcTan[TanBeta]]}
 }
};


BoundaryMatchingScaleUp={
{
  {Yu, Sqrt[1+TanBeta^2]/TanBeta*Yu},
  {Yd, Sqrt[1+TanBeta^2]*Yd},
  {Ye, Sqrt[1+TanBeta^2]*Ye}
 }
};

BoundaryMatchingScaleDown={
{
  {\[Lambda], -1/3 EFTcoupNLO[hh[1].hh[1].hh[1].hh[1]]}
 }
};

ParametersToSolveTadpoleMatchingScale={
 {mHd2,mHu2}
};
 


