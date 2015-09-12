ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    

{AlphaS,    {Description -> "Alpha Strong"}},	
{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

{e,         { Description -> "electric charge"}}, 

{Yu,        { Description -> "Up-Yukawa-Coupling"   }}, 
{Yd,        { Description -> "Down-Yukawa-Coupling"}},
{Ye,        { Description -> "Lepton-Yukawa-Coupling"}}, 

{T[Ye],     { Description -> "Trilinear-Lepton-Coupling"}},
{T[Yd],     { Description -> "Trilinear-Down-Coupling"}}, 
{T[Yu],     { Description -> "Trilinear-Up-Coupling"}}, 

{\[Mu],     { Description -> "Mu-parameter"}}, 
{B[\[Mu]],  { Description -> "Bmu-parameter"}},        

{mq2,       { Description -> "Softbreaking left Squark Mass"}},
{me2,       { Description -> "Softbreaking right Slepton Mass"}},
{ml2,       { Description -> "Softbreaking left Slepton Mass"}},
{mu2,       { Description -> "Softbreaking right Up-Squark Mass"}},
{md2,       { Description -> "Softbreaking right Down-Squark Mass"}},
{mHd2,      { Description -> "Softbreaking Down-Higgs Mass"}}, 
{mHu2,      { Description -> "Softbreaking Up-Higgs Mass"}}, 

{MassB,     { Description -> "Bino Mass parameter"}},
{MassWB,    { Description -> "Wino Mass parameter"}},
{MassG,     { Description -> "Gluino Mass parameter"}},
                           
{vd,        { Description -> "Down-VEV"}}, 
{vu,        { Description -> "Up-VEV"}},       
{v,         { Description -> "EW-VEV"}},
             
{\[Beta],   { Description -> "Pseudo Scalar mixing angle"  }},             
{TanBeta,   { Description -> "Tan Beta" }},              
{\[Alpha],  { Description -> "Scalar mixing angle" }},                          
             
    
{ZD,        { Description -> "Down-Squark-Mixing-Matrix" }},             
{ZU,        { Description -> "Up-Squark-Mixing-Matrix"}},             
{ZE,        { Description -> "Slepton-Mixing-Matrix"}}, 
{ZV,        { Description->"Sneutrino Mixing-Matrix"}},   
{ZH,        { Description->"Scalar-Mixing-Matrix", 
              Real -> False,
				  Dependence->None,
  				  DependenceOptional->None,
  				  DependenceNum->None,
              LesHouches -> NMHMIX      }},
{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix", 
               Real ->False,
				  Dependence->None,
  				  DependenceOptional->None,
  				  DependenceNum->None,
               LesHouches -> NMAMIX     }},
{ZP,        { Description->"Charged-Mixing-Matrix"}},                      
                                          
 
{ZN,        { Description->"Neutralino Mixing-Matrix",
              LesHouches -> NMNMIX   }}, 

{UP,        { Description->"Chargino-plus Mixing-Matrix"}}, 
{UM,        { Description->"Chargino-minus Mixing-Matrix"}}, 

{ZEL,       { Description ->"Left-Lepton-Mixing-Matrix"}},
{ZER,       { Description ->"Right-Lepton-Mixing-Matrix" }},                          
{ZDL,       { Description ->"Left-Down-Mixing-Matrix"}},                       
{ZDR,       { Description ->"Right-Down-Mixing-Matrix"}},              
{ZUL,       { Description ->"Left-Up-Mixing-Matrix"}},                        
{ZUR,       { Description ->"Right-Up-Mixing-Matrix"}},           
              
{ThetaW,    { Description -> "Weinberg-Angle"}},                 
{ThetaWp,  { Description -> "Theta'",
			 Dependence -> None,
 			 DependenceNum -> None  }},

{ZZ, {Description ->   "Photon-Z-Z' Mixing Matrix"}},
{ZW, {Description -> "W Mixing Matrix" }},
{ZfW, {Description ->    "Wino Mixing Matrix"}},
                        
{PhaseGlu,  { Description -> "Gluino-Phase" }},   
                                                                      
{\[Kappa],   {Description -> "Singlet Self-Interaction"}},                               
{T[\[Kappa]],  { Description -> "Softbreaking Singlet Self-Interaction" }}, 
{\[Lambda],   { Description -> "Singlet-Higgs-Interaction"   }},                               
{T[\[Lambda]],  {Description -> "Softbreaking Singlet-Higgs-Interaction"}},        
             
{ms2,       { Description -> "Softbreaking Singlet Mass" }},


{ms12,       { LaTeX -> "m^2_{s_1}",
              OutputName -> ms12,
              LesHouches -> {MSOFT,200} }},

{ms22,       { LaTeX -> "m^2_{s_2}",
              OutputName -> ms22,
              LesHouches -> {MSOFT,201} }},


{ms32,       { LaTeX -> "m^2_{s_3}",
              OutputName -> ms32,
              LesHouches -> {MSOFT,202} }},


{vS,        { Description -> "Singlet-VEV"}},


{gp,       {Description -> "U(1)' Gauge Coupling"}},


{Qq,        { LaTeX -> "Q_q",
              OutputName -> Qq,
              Real -> True,
			  LesHouches -> {XCharge,1}}},             
{Ql,        { LaTeX -> "Q_l",
              OutputName -> Ql,
              Real -> True,
			  LesHouches -> {XCharge,2}}},
{QHu,        { LaTeX -> "Q_{H_u}",
              OutputName -> QHu,
              Real -> True,
			  LesHouches -> {XCharge,3}}},
{QHd,        { LaTeX -> "Q_{H_d}",
              OutputName -> QHd,
              Real -> True,
			  LesHouches -> {XCharge,4}}},
{Qd,        { LaTeX -> "Q_d",
              OutputName -> Qd,
              Real -> True,
			  LesHouches -> {XCharge,5}}},
{Qu,        { LaTeX -> "Q_u",
              OutputName -> Qu,
              Real -> True,
			  LesHouches -> {XCharge,6}}},             
{Qe,        { LaTeX -> "Q_e",
              OutputName -> Qe,
              Real -> True,
			  LesHouches -> {XCharge,7}}},
{Qs,        { LaTeX -> "Q_s",
              OutputName -> Qs,
              Real -> True,
			  LesHouches -> {XCharge,8}}},
{Qs1,     {LaTeX -> "Q_1",
              OutputName -> Qs1,
                 Real ->True,
			  LesHouches -> {XCharge,9} }},
{Qs2,     {LaTeX -> "Q_2",
              OutputName -> Qs2,
                 Real ->True,
			  LesHouches -> {XCharge,10} }},
{Qs3,     {LaTeX -> "Q_3",
              OutputName -> Qs3,
                 Real ->True,
			  LesHouches -> {XCharge,11} }},
                 
{MassU,        { LaTeX -> "M_U",
                 OutputName -> M4,
                 LesHouches -> {MSOFT,4} }},
{ThetaZ,    { LaTeX -> "\\Theta_Z"}},

{vS1,      { LaTeX -> "v_1",
              OutputName -> vS1,
              LesHouches -> {HMIX,200}}},
{vS2,      { LaTeX -> "v_2",
              OutputName -> vS2,
              LesHouches -> {HMIX,201}}},
{vS3,      { LaTeX -> "v_3",
              OutputName -> vS3,
              LesHouches -> {HMIX,202}}}

 }; 
 


