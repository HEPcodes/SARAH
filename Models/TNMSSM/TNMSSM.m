Off[General::spell]

Model`Name = "TNMSSM";
Model`NameLaTeX ="Triplet extended NMSSM";
Model`Authors = "F.Staub";
Model`Date = "2015-11-16";


(* 2013-09-01: changing to new conventions for Superfields, Superpotential and global symmetries *)
(* 2015-11-16: changed SPheno.m *)

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global symmetries *)


Global[[1]]={Z[2],RParity}; 
Global[[2]]={Z[3],Z3}; 

(* Charges under Global symmetries *)
RpM = {-1,-1,1}; RpP = {1,1,-1};
Z3t = Exp[2*Pi*\[ImaginaryI]/3];

(* Vector Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False, RpM, 1};
Gauge[[2]]={WB, SU[2], left,        g2,True,  RpM, 1};
Gauge[[3]]={G,  SU[3], color,       g3,False, RpM, 1};


(* Chiral Superfields *)

SuperFields[[1]] = {q, 3, {uL,  dL},    1/6, 2, 3, RpM, Z3t};  
SuperFields[[2]] = {l, 3, {vL,  eL},   -1/2, 2, 1, RpM, Z3t};
SuperFields[[3]] = {Hd,1, {Hd0, Hdm},  -1/2, 2, 1, RpP, Z3t};
SuperFields[[4]] = {Hu,1, {Hup, Hu0},   1/2, 2, 1, RpP, Z3t};

SuperFields[[5]] = {d, 3, conj[dR],   1/3, 1, -3, RpM, Z3t};
SuperFields[[6]] = {u, 3, conj[uR],  -2/3, 1, -3, RpM, Z3t};
SuperFields[[7]] = {e, 3, conj[eR],     1, 1,  1, RpM, Z3t};

SuperFields[[8]] = {s, 1, sR,     0, 1,  1, RpP, Z3t};

SuperFields[[9]] = {t,   1, {{Tp/Sqrt[2],-Tpp},{T0, -Tp/Sqrt[2]}},   1, 3, 1, RpP, Z3t};
SuperFields[[10]] = {tb, 1, {{Tm/Sqrt[2],-T0b},{Tmm, -Tm/Sqrt[2]}}, -1, 3, 1, RpP, Z3t};


(*------------------------------------------------------*)
(* Superpotential *)
(*------------------------------------------------------*)

SuperPotential = Yu u.q.Hu - Yd d.q.Hd - Ye e.l.Hd + \[Lambda] Hu.Hd.s + \[Kappa]/3 s.s.s \
                 + LamT s.tb.t + XiU Hu.tb.Hu + XiD Hd.t.Hd;



(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)


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
  Fv -> {FvL, 0},
  FS -> {FsR, conj[FsR]},
  FT -> {FT0, conj[FT0b]},
  FTC -> {FTm, conj[FTp]},
  FTCC -> {FTmm, conj[FTpp]}
};


(* ----- After EWSB ----- *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW},
  {{fWB[1],fWB[2],fWB[3]},{fWm,fWp,fW0},ZfW}
};      
        
 
DEFINITION[EWSB][VEVs]= 
{    {SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid,1/Sqrt[2]}},
     {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu,1/Sqrt[2]}},
     {SsR, {vS, 1/Sqrt[2]}, {sigmaS, \[ImaginaryI]/Sqrt[2]},{phiS,1/Sqrt[2]}},
     {ST0, {vT, 1/Sqrt[2]}, {sigmaT, \[ImaginaryI]/Sqrt[2]},{phiT,1/Sqrt[2]}},
     {ST0b, {vTb, 1/Sqrt[2]}, {sigmaTb, \[ImaginaryI]/Sqrt[2]},{phiTb,1/Sqrt[2]}}     };


 
DEFINITION[EWSB][MatterSector]= 

{    {{SdL, SdR}, {Sd, ZD}},
     {{SvL}, {Sv, ZV}},
     {{SuL, SuR}, {Su, ZU}},
     {{SeL, SeR}, {Se, ZE}},
     {{phid, phiu, phiS,phiT,phiTb}, {hh, ZH}},
     {{sigmad, sigmau,sigmaS,sigmaT,sigmaTb}, {Ah, ZA}},
     {{SHdm,conj[SHup],STm,conj[STp]},{Hpm,ZP}},
     {{STmm,conj[STpp]},{Hmm,ZMM}},
     {{fB, fW0, FHd0, FHu0,FsR,FT0,FT0b}, {L0, ZN}}, 
     {{{fWm, FHdm,FTm}, {fWp, FHup,FTp}}, {{Lm,UM}, {Lp,UP}}},
     {{{FeL},{conj[FeR]}},{{FEL,ZEL},{FER,ZER}}},
     {{{FdL},{conj[FdR]}},{{FDL,ZDL},{FDR,ZDR}}},
     {{{FuL},{conj[FuR]}},{{FUL,ZUL},{FUR,ZUR}}}            
   }; 

DEFINITION[EWSB][Phases]= 
{    {fG, PhaseGlu},
     {FTmm, PhaseTmm}
    }; 


DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  FDL, conj[FDR]},
 Fe ->{  FEL, conj[FER]},
 Fu ->{  FUL, conj[FUR]},
 Fv ->{  FvL, 0},
 Chi ->{ L0, conj[L0]},
 Cha ->{ Lm, conj[Lp]},
 Glu ->{ fG, conj[fG]},
 DCha -> {FTmm, conj[FTpp]}
};
