Off[General::spell]

Model`Name = "BLSSM";
Model`NameLaTeX ="B-L-SSM";
Model`Authors = "L.Basso, F.Staub";
Model`Date = "2012-09-01";

(* 2013-09-01: changing to new conventions for Superfields, Superpotential and global symmetries *)



(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2],MParity};
MpM = {-1,-1,1};
MpP = {1,1,-1};

(* Vector Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,  False, MpM};
Gauge[[2]]={WB, SU[2], left,        g2,  True,  MpM};
Gauge[[3]]={G,  SU[3], color,       g3,  False, MpM};
Gauge[[4]]={Bp,  U[1], BminusL,     gBL, False, MpM};

(* Chiral Superfields *)

SuperFields[[1]] = {q, 3, {uL,  dL},    1/6, 2, 3,  1/6, MpM};  
SuperFields[[2]] = {l, 3, {vL,  eL},   -1/2, 2, 1, -1/2, MpM};
SuperFields[[3]] = {Hd,1, {Hd0, Hdm},  -1/2, 2, 1,    0, MpP};
SuperFields[[4]] = {Hu,1, {Hup, Hu0},   1/2, 2, 1,    0, MpP};

SuperFields[[5]] = {d, 3, conj[dR],   1/3, 1, -3, -1/6, MpM};
SuperFields[[6]] = {u, 3, conj[uR],  -2/3, 1, -3, -1/6, MpM};
SuperFields[[7]] = {e, 3, conj[eR],     1, 1,  1,  1/2, MpM};
SuperFields[[8]] = {vR,3, conj[vR],     0, 1,  1,  1/2, MpM};

SuperFields[[9]]  = {C1, 1, C10,  0, 1, 1, -1, MpP};
SuperFields[[10]] = {C2, 1, C20,  0, 1, 1,  1, MpP};


(*------------------------------------------------------*)
(* Superpotential *)
(*------------------------------------------------------*)

SuperPotential = Yu u.q.Hu - Yd d.q.Hd - Ye e.l.Hd + \[Mu] Hu.Hd + Yv vR.l.Hu - MuP C1.C2 + Yn vR.C1.vR;


(*----------------------------------------------*)
(*   DEFINITION                                  *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};


DEFINITION[GaugeES][DiracSpinors]={
  Bino ->{fB, conj[fB]},
  Wino -> {fWB, conj[fWB]},
  Glu -> {fG, conj[fG]},
  H0 -> {FHd0, conj[FHu0]},
  HC -> {FHdm, conj[FHup]},
  Fd1 -> {FdL, 0},
  Fd2 -> {0, FdR},
  Fu1 -> {FuL, 0},
  Fu2 -> {0, FuR},
  Fe1 -> {FeL, 0},
  Fe2 -> {0, FeR},
  Fv1 -> {FvL, 0},
  Fv2 -> {0, FvR},
  FC -> {FC10, conj[FC20]},
  FB -> {fBp, conj[fBp]}
};



(*--- Gauge Sector ---- *)



DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3],VBp},{VP,VZ,VZp},ZZ},
  {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW},
  {{fWB[1],fWB[2],fWB[3]},{fWm,fWp,fW0},ZfW}
};


       
(*--- VEVs ---- *)
        
          	
DEFINITION[EWSB][VEVs]= 
{    {SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid,1/Sqrt[2]}},
     {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu,1/Sqrt[2]}},
     {SvL, {0, 0}, {sigmaL, \[ImaginaryI]/Sqrt[2]},{phiL,1/Sqrt[2]}},
   {SvR, {0, 0}, {sigmaR, \[ImaginaryI]/Sqrt[2]},{phiR,1/Sqrt[2]}},
{SC10, {x1, 1/Sqrt[2]}, {sigma1, \[ImaginaryI]/Sqrt[2]},{phi1, 1/Sqrt[2]}},
{SC20, {x2, 1/Sqrt[2]}, {sigma2, \[ImaginaryI]/Sqrt[2]},{phi2, 1/Sqrt[2]}}     
};
 

 
(*--- Matter Sector ---- *)
 
 
DEFINITION[EWSB][MatterSector]= 
{    {{SdL, SdR}, {Sd, ZD}},
     {{SuL, SuR}, {Su, ZU}},
     {{SeL, SeR}, {Se, ZE}},
     {{sigmaL,sigmaR}, {SvIm, ZVI}},
     {{phiL,phiR}, {SvRe, ZVR}}, 
     {{phid, phiu,phi1, phi2}, {hh, ZH}}, 
     {{sigmad, sigmau,sigma1,sigma2}, {Ah, ZA}},
     {{SHdm,conj[SHup]},{Hpm,ZP}},
     {{fB, fW0, FHd0, FHu0,fBp,FC10,FC20}, {L0, ZN}}, 
     {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}},
     {{FvL,conj[FvR]},{Fvm,UV}},
     {{{FeL},{conj[FeR]}},{{FEL,ZEL},{FER,ZER}}},
     {{{FdL},{conj[FdR]}},{{FDL,ZDL},{FDR,ZDR}}},
     {{{FuL},{conj[FuR]}},{{FUL,ZUL},{FUR,ZUR}}}                    \
       }; 



DEFINITION[EWSB][Phases]= 
{    {fG, PhaseGlu}
    }; 

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  FDL, conj[FDR]},
 Fe ->{  FEL, conj[FER]},
 Fu ->{  FUL, conj[FUR]},
 Fv ->{  Fvm, conj[Fvm]},
 Chi ->{ L0, conj[L0]},
 Cha ->{ Lm, conj[Lp]},
 Glu ->{ fG, conj[fG]}
};	
