OnlyLowEnergySPheno = True;

MINPAR={
        {1,LHInput},
        {2,LTInput},
        {3,LT2Input},
        {4,LTHInput},
        {5,LTH2Input},        
        {6,KInput},

        {10,vTInput}

       };

ParametersToSolveTadpoles = {mu2,MT};

DEFINITION[MatchingConditions]= {
 {Ye, YeSM},
 {Yd, YdSM},
 {Yu, YuSM},
 {g1, g1SM},
 {g2, g2SM},
 {g3, g3SM},
 {v,  vSM}
 };


BoundaryLowScaleInput={
  {vT, vTInput},
  {LH, LHInput},
  {LT, LTInput},
  {LT2, LT2Input},
  {LHT, LTHInput},
  {LHT2, LTH2Input},
  {KHT, KInput}
};


ListDecayParticles = {Fu,Fe,Fd,hh,Ah,Hpm};
ListDecayParticles3B = {{Fu,"Fu.f90"},{Fe,"Fe.f90"},{Fd,"Fd.f90"}};

DefaultInputValues ={LHInput -> 0.1, LTInput -> 0.2 LT2Input -> 0.1, LTHInput -> 0.2, LTH2Input->0.2, KInput->10., vTInput->3 };


RenConditionsDecays={
{dCosTW, 1/2*Cos[ThetaW] * (PiVWp/(MVWp^2) - PiVZ/(mVZ^2)) },
{dSinTW, -dCosTW/Tan[ThetaW]},
{dg2, 1/2*g2*(derPiVPheavy0 + PiVPlightMZ/MVZ^2 - (-(PiVWp/MVWp^2) + PiVZ/MVZ^2)/Tan[ThetaW]^2 + (2*PiVZVP*Tan[ThetaW])/MVZ^2)  },
{dg1, dg2*Tan[ThetaW]+g2*dSinTW/Cos[ThetaW]- dCosTW*g2*Tan[ThetaW]/Cos[ThetaW]}
};
