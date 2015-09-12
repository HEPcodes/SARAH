
ParticleDefinitions[GaugeES] = {
      {SdL0,  {  Description -> "Left Down-Squarks"}},
      {SdR0,  { Description -> "Right Down-Squarks"}},
      {SuL0,  { Description -> "Left Up-Squarks"}}, 
      {SuR0,  { Description -> "Right Up-Squarks" }},                
      {SeL0,  { Description -> "Left Selectron"}}, 
      {SeR0,  { Description -> "Right Selectron"}},                         
      {SvL0,  { Description -> "Left Sneutrino"}},                         
      {SHd0, { Description -> "Neutral Down-Higgs"}},                         
      {SHdm, { Description -> "Charged Down-Higgs"}}, 
      {SHu0, { Description -> "Neutral Up-Higgs"}},                         
      {SHup, { Description -> "Charged Up-Higgs"}},                                                   
      {VB,   { Description -> "B-Boson"}},                                                   
      {VG,   { Description -> "Gluon"}},          
      {VWB,  { Description -> "W-Bosons"}},          
      {gB,   { Description -> "B-Boson Ghost"}},                                                   
      {gG,   { Description -> "Gluon Ghost" }},          
      {gWB,  { Description -> "W-Boson Ghost"}},    
      {Glu,  { Description -> "Gluino"}},
      {Wino, { Description -> "Wino"}},
      {Bino, { Description -> "Bino"}},                                                                                           
      {H0,   { Description -> "Neutral Higgsinos"}},
      {HC,   { Description -> "Charged Higgsinos"}},
      
      {Fd1,  { Description -> "Dirac Left Down-Quark"}},
      {Fd2,  { Description -> "Dirac Right Down-Quark"}},
      {Fu1,  { Description -> "Dirac Left Up-Quark"}},
      {Fu2,  { Description -> "Dirac Right Up-Quark"}},
      {Fe1,  { Description -> "Dirac Left Electron"}},
      {Fe2,  { Description -> "Dirac Right Electron"}},
      {Fv,   { Description -> "Neutrinos" }},
      
      
      {SbL, {LaTeX ->"\\tilde{b}_L"}},
      {SbR, {LaTeX ->"\\tilde{b}_R"}},
      {SsL, {LaTeX ->"\\tilde{s}_L"}},
      {SsR, {LaTeX ->"\\tilde{s}_R"}},
      {SdL, {LaTeX ->"\\tilde{d}_L"}},
      {SdR, {LaTeX ->"\\tilde{d}_R"}},
      {StL, {LaTeX ->"\\tilde{t}_L"}},
      {StR, {LaTeX ->"\\tilde{t}_R"}},
      {ScL, {LaTeX ->"\\tilde{c}_L"}},
      {ScR, {LaTeX ->"\\tilde{c}_R"}},
      {SuL, {LaTeX ->"\\tilde{u}_L"}},
      {SuR, {LaTeX ->"\\tilde{u}_R"}},
      {SeL, {LaTeX ->"\\tilde{e}_L"}},
      {SeR, {LaTeX ->"\\tilde{e}_R"}},
      {SmL, {LaTeX ->"\\tilde{\\mu}_L"}},
      {SmR, {LaTeX ->"\\tilde{\\mu}_R"}},
      {StauL, {LaTeX ->"\\tilde{\\tau}_L"}},
      {StauR, {LaTeX ->"\\tilde{\\tau}_R"}}
      
      
      
      
      
      };
      
      
      
      
  ParticleDefinitions[EWSB] = {
             
                                             
      
      {Sd ,  {  Description->"Down Squark",
                PDG ->  {1000001,2000015},
                Mass -> LesHouches,
                FeynArtsNr -> 141,   
                LaTeX -> "\\tilde{d}",
                OutputName -> "sd" }},
                
       {Ss ,  {  Description->"Strange Squark",
                PDG ->  {1000003,2000003},
                Mass -> LesHouches,
                FeynArtsNr -> 142,   
                LaTeX -> "\\tilde{s}",
                OutputName -> "ss" }},        
               
       {Sb ,  {  Description->"Bottom Squark",
                PDG ->  {1000005,2000005},
                Mass -> LesHouches,
                FeynArtsNr -> 143,   
                LaTeX -> "\\tilde{b}",
                OutputName -> "sb" }},            
              
           
      {Su ,  {  Description->"Up Squark",
                PDG ->  {1000002,2000002},
                Mass -> LesHouches,
                FeynArtsNr -> 131,   
                LaTeX -> "\\tilde{u}",
                OutputName -> "su" }}, 
                
     {Sc ,  {  Description->"Charmed Squark",
                PDG ->  {1000004,2000004},
                Mass -> LesHouches,
                FeynArtsNr -> 132,   
                LaTeX -> "\\tilde{c}",
                OutputName -> "sc" }},    
      
      {St ,  {  Description->"Top Squark",
                PDG ->  {1000006,2000006},
                Mass -> LesHouches,
                FeynArtsNr -> 133,   
                LaTeX -> "\\tilde{t}",
                OutputName -> "st" }},                        
      
      
      {Se ,  { Description -> "Sleptons"}}, 
      {Sv ,  { Description -> "Sneutrinos"}}, 
                
                                          
      {hh ,  { Description -> "Higgs"}}, 
      {Ah ,  { Description -> "Pseudo-Scalar Higgs"}},                       
      {Hpm,  { Description -> "Charged Higgs"}},                                                   
      {VP,   { Description -> "Photon"}}, 
      {VZ,   { Description -> "Z-Boson" }}, 
      {VG,   { Description -> "Gluon" }},          
      {VWm,  { Description -> "W-Boson" }},         
      {gP,   { Description -> "Photon Ghost"}},                                                   
      {gWm,  { Description -> "Negative W-Boson Ghost"}}, 
      {gWmC, { Description -> "Positive W-Boson Ghost" }}, 
      {gZ,   { Description -> "Z-Boson Ghost" }},    
      {gG,   { Description -> "Gluon Ghost" }},         
      {Glu,  { Description -> "Gluino" }},
      {Chi,  { Description -> "Neutralinos"}},
      {Cha,  { Description -> "Charginos"}}  ,          
                               
                 
      {Fd,    { Description->"Down Quark",
                 PDG -> {1},
                 Mass -> 0.0035,
                 FeynArtsNr -> 41,
                 LaTeX -> "d",
                 OutputName -> "d" }},   
                 
     {Fs,    { Description->"Strange Quark",
                 PDG -> {3},
                 Mass -> 0.104,
                 FeynArtsNr -> 42,
                 LaTeX -> "s",
                 OutputName -> "s" }},   
                 
      {Fb,    { Description->"Bottom Quark",
                 PDG -> {5},
                 Mass -> 4.20,
                 FeynArtsNr -> 43,
                 LaTeX -> "b",
                 OutputName -> "b" }},                               
                 
      {Fu,    { Description->"Up Quark",
                 PDG -> {2},
                 Width -> 0, 
                 Mass -> 0.0015,
                 FeynArtsNr -> 31,
                 LaTeX -> "u",
                 OutputName -> "u" }},   
                 
       {Fc,    { Description->"Charmed Quark",
                  PDG -> {4},
                 Width ->  0, 
                 Mass -> 1.27,
                 FeynArtsNr -> 32,
                 LaTeX -> "c",
                 OutputName -> "c" }}, 
                 
      {Ft,    { Description->"Top Quark",
                 PDG -> {6},
                 Width -> 1.59, 
                 Mass -> 171.2,
                 FeynArtsNr -> 33,
                 LaTeX -> "t",
                 OutputName -> "t" }},                            
                 
       {Fe,   { Description -> "Leptons" }},                              
       {Fv,   { Description -> "Neutrinos" }}    
                                                                                                          
     
        };    
        
        
        
 WeylFermionAndIndermediate = {
       {FHd0, { Description -> "Neutral Down-Higgsino"}},      
       {FHu0, { Description -> "Neutral Up-Higgsino" }}, 
       {FHdm, { Description -> "Charged Down-Higgsino"}},                                                                                                             
       {FHup, { Description -> "Charged Up-Higgsino"}},
       {L0,   { Description -> "Neutralino Weyl-Spinor"}},
       {Lm,   { Description -> "Negative Chargino Weyl-Spinor"}},
       {Lp,   { Description -> "Positive Chargino Weyl-Spinor"}}, 
       {fG,   { Description ->"Gluino Weyl-Spinor"}},
       {fWB,  { Description ->"Wino Weyl-Spinor"}},
       {fW0,  { Description ->"Neutral Wino" }},
       {fWm,  { Description ->"Negative Wino"}},                 
       {fWp,  { Description ->"Positive Wino"}},                 
       {fB,   { Description ->"Bino Weyl-Spinor"}},    
       {phid, { Description -> "Scalar Down"}},                                                                       
       {phiu, { Description -> "Scalar Up"}}, 
       {sigmad,   { Description -> "Pseudo Scalar Down"}},                                                                                        
       {sigmau,   { Description -> "Pseudo Scalar Up" }},
       {SHd,  { Description -> "Down-Higgs"}},
       {SHu,  { Description -> "Up-Higgs"}},
       {Sl,   { Description -> "Left Slepton" }},
       {Sq,   { Description -> "Left Squark" }},
       {FdL0,  { Description -> "Left Down-Quark" }},
       {FdR0,  { Description -> "Right Down-Quark" }},              
       {FuL0,  { Description -> "Left Up-Quark" }},
       {FuR0,  { Description -> "Right Up-Quark" }},
       {FeL0,  { Description -> "Left Electron" }},
       {FeR0,  { Description -> "Right Electron" }},              

       {FeL,  { LaTeX -> "e_L" }},
       {FeR,  { LaTeX -> "e_R" }},
       {FdL,  { LaTeX -> "d_L" }},
       {FdR,  { LaTeX -> "d_R" }},              
       {FuL,  { LaTeX -> "u_L" }},
       {FuR,  { LaTeX -> "u_R" }},
       
       {FvL0,  { LaTeX -> "\\nu" }},
       {FveL,  { LaTeX -> "\\nu_e" }},
       {FvmL,  { LaTeX -> "\\nu_{\\mu}" }},
       {FvtL,  { LaTeX -> "\\nu_{\\tau}" }},


       {FmL,  { LaTeX -> "\\mu_L" }},
       {FmR,  { LaTeX -> "\\mu_R" }},
       {FsL,  { LaTeX -> "s_L" }},
       {FsR,  { LaTeX -> "s_R" }},              
       {FcL,  { LaTeX -> "c_L" }},
       {FcR,  { LaTeX -> "c_R" }},

       {FtauL,  { LaTeX -> "\\tau_L" }},
       {FtauR,  { LaTeX -> "\\tau_R" }},
       {FbL,  { LaTeX -> "b_L" }},
       {FbR,  { LaTeX -> "b_R" }},              
       {FtL,  { LaTeX -> "t_L" }},
       {FtR,  { LaTeX -> "t_R" }},




       {FHd,  { Description -> "Down-Higgsino" }},
       {FHu,  { Description -> "Up-Higgsino" }},
       {Fl,   { Description -> "Left Leptons"}},
       {Fq,   { Description -> "Left Quarks"}},  
       {FvL,  { Description -> "Left Neutrino"}},                  

       {e,    { Description -> "Right Electron Superfield" }},
       {d,    { Description -> "Right Down-Quark Superfield" }},                 
       {q,    { Description -> "Left Quark Superfield" }},                 
       {u,    { Description -> "Right Up-Quark Superfield" }},                 
       {l,    { Description -> "left Lepton Superfield" }},  
       {Hd,   { Description -> "Down-Higgs Superfield"}},                                         
       {Hu,   { Description -> "Up-Higgs Superfield" }},                 
       {G,    { Description -> "Gluon Superfield" }},                 
       {B,    { Description -> "B Superfield" }},                 
       {WB,   { Description -> "W Superfield" }}


        };       


