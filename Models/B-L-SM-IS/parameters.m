ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},

{g11p,        {Description -> "Mixed Gauge Coupling 2"}},
{g1p1,        {Description -> "Mixed Gauge Coupling 1"}},
{g1p,       {   Description -> "B-L-Coupling"}},

{MZp,       {   Description -> "Z' mass"}},


{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    
{AlphaS,    {Description -> "Alpha Strong"}},	
{e,         { Description -> "electric charge"}},
{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},
 

{Yu,        { Description -> "Up-Yukawa-Coupling",
			 DependenceNum ->  Sqrt[2]/vH* {{Mass[Fu,1],0,0},
             									{0, Mass[Fu,2],0},
             									{0, 0, Mass[Fu,3]}}}}, 
             									
{Yd,        { Description -> "Down-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/vH* {{Mass[Fd,1],0,0},
             									{0, Mass[Fd,2],0},
             									{0, 0, Mass[Fd,3]}}}},
             									
{Ye,        { Description -> "Lepton-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/vH* {{Mass[Fe,1],0,0},
             									{0, Mass[Fe,2],0},
             									{0, 0, Mass[Fe,3]}}}}, 
                                                                            
                                                                           
{Mu,         { Description -> "SM Mu Parameter"}},                                        
{\[Lambda],  { Description -> "SM Higgs Selfcouplings"}},
{vH,          { Description -> "EW-VEV",
               DependenceSPheno -> None }},
{vX,      {  LaTeX -> "x",
             Dependence ->  None, 
             OutputName -> vX,
             Real -> True,
             LesHouches -> {BL,43} }},

{ThetaW,    { Description -> "Weinberg-Angle"}},
{ThetaWp,  { Description -> "Theta'", DependenceNum -> None }},

{ZZ, {Description ->   "Photon-Z-Z' Mixing Matrix"}},
{ZW, {Description -> "W Mixing Matrix"}},

{Ys, {OutputName -> YS,
      LaTeX -> "Y_S",
      LesHouches -> YS}},

{Yx, {OutputName -> YX,
      LaTeX -> "Y_x",
      LesHouches -> YX}},

{Yv, {OutputName -> Yv,
      LaTeX -> "Y_\\nu",
      LesHouches -> Yv}},

{L1, {OutputName -> lam1,
      LaTeX -> "\\lambda_1",
      LesHouches -> {BL,1}}},


{L2, {OutputName -> lam3,
      LaTeX -> "\\lambda_3",
      LesHouches -> {BL,3}}},

{L3, {OutputName -> lam5,
      LaTeX -> "\\lambda_5",
      LesHouches -> {BL,5}}},

{Mu1, {OutputName -> Mu1,
      LaTeX -> "\\mu_1",
      LesHouches -> {BL,10}}},

{Mu2, {OutputName -> Mu2,
      LaTeX -> "\\mu_2",
      LesHouches -> {BL,11}}},

{mH2, {OutputName -> mH2,
      LaTeX -> "\\m^2_H",
      LesHouches -> {BL,20}}},

{mchi2, {OutputName -> mX2,
      LaTeX -> "\\m^2_\\chi",
      LesHouches -> {BL,21}}},

{MUS, {OutputName -> muS,
      LaTeX -> "\\mu_s",
      LesHouches -> MUS}},

{MUS1, {OutputName -> muS1,
      LaTeX -> "\\mu_s1",
      LesHouches -> MUS1}},


      
{Vu,        {Description ->"Left-Up-Mixing-Matrix"}},
{Vd,        {Description ->"Left-Down-Mixing-Matrix"}},
{Uu,        {Description ->"Right-Up-Mixing-Matrix"}},
{Ud,        {Description ->"Right-Down-Mixing-Matrix"}}, 
{Ve,        {Description ->"Left-Lepton-Mixing-Matrix"}},
{Ue,        {Description ->"Right-Lepton-Mixing-Matrix"}},

{SS,        {OutputName -> SS,
	    Description->"S matrix", 
            Dependence -> None,
             DependenceOptional -> None,
               DependenceNum -> None,
             LesHouches ->SSMIX
       }},



{ZM,	    {Description -> "Neutrino-Mixing-Matrix"}},

{ZH,        { Description->"Scalar-Mixing-Matrix", 
               Dependence -> None,
               DependenceOptional -> None,
               DependenceNum -> None   }},
{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix", 
                Dependence -> None,
               DependenceOptional -> None,
               DependenceNum -> None   }}



 }; 
 

