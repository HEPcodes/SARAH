ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    
{AlphaS,    {Description -> "Alpha Strong"}},	
{e,         { Description -> "electric charge"}}, 
{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},


{Yu,        { Description -> "Up-Yukawa-Coupling",
              Form -> Diagonal,
              DependenceNum ->  Sqrt[2]/vu  {{Mass[Fu,1],0,0},
             							    {0, Mass[Fc,1],0},
             								{0, 0, Mass[Ft,1]}} }}, 
{Yd,        { Description -> "Down-Yukawa-Coupling",
              Form -> Diagonal, 
             DependenceNum ->  Sqrt[2]/vd {{Mass[Fd,1],0,0},
             								{0, Mass[Fs,1],0},
             								{0, 0, Mass[Fb,1]}} }},
             								
{Ye,        { Description -> "Lepton-Yukawa-Coupling"}}, 

{T[Ye],     { Description -> "Trilinear-Lepton-Coupling"}},
{T[Yd],     { Description -> "Trilinear-Down-Coupling",
              Form -> Diagonal}}, 
{T[Yu],     { Description -> "Trilinear-Up-Coupling",
              Form -> Diagonal}}, 

{\[Mu],     { Description -> "Mu-parameter",
              OutputName -> "Mue"}}, 
{B[\[Mu]],  { Description -> "Bmu-parameter"}},        

{mq2,       { Description -> "Softbreaking left Squark Mass",
              Form -> Diagonal}},
{me2,       { Description -> "Softbreaking right Slepton Mass",
              Form -> Diagonal}},
{ml2,       { Description -> "Softbreaking left Slepton Mass",
              Form -> Diagonal}},
{mu2,       { Description -> "Softbreaking right Up-Squark Mass",
              Form -> Diagonal}},
{md2,       { Description -> "Softbreaking right Down-Squark Mass",
              Form -> Diagonal}},
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
             
    
{ZH,        { Description->"Scalar-Mixing-Matrix"}},
{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix"}},
{ZP,        { Description->"Charged-Mixing-Matrix"}},                      
                                          
 
{ZN,        { Description->"Neutralino Mixing-Matrix" }}, 
{UP,        { Description->"Chargino-plus Mixing-Matrix"}}, 
{UM,        { Description->"Chargino-minus Mixing-Matrix"}}, 

             
{ThetaW,    { Description -> "Weinberg-Angle"}},                           
{PhaseGlu,  { Description -> "Gluino-Phase" }},   

{ZZ, {Description ->   "Photon-Z Mixing Matrix"}},
{ZW, {Description -> "W Mixing Matrix" }},
{ZfW, {Description ->    "Wino Mixing Matrix"}},
       
       
{ZD,        { LaTeX -> "Z^D",
             OutputName -> ZD, 
             LesHouches -> sdownmix			  }},
             
{ZS,        { LaTeX -> "Z^S",
             OutputName -> ZS, 
             LesHouches -> sstrmix			  }},             
             
{ZB,        { LaTeX -> "Z^B",
             OutputName -> ZB, 
             LesHouches -> sbotmix			  }},
     
         
             
{ZU,        { LaTeX -> "Z^U",
             OutputName -> ZU, 
             LesHouches ->  supmix			       }},
             
         
             
{ZC,        { LaTeX -> "Z^C",
             OutputName -> ZC, 
             LesHouches ->  scharmmix			       }},


         
             
{ZT,        { LaTeX -> "Z^T",
             OutputName -> ZT, 
             LesHouches ->  stopmix			       }},
             
             
{ZE,        { LaTeX -> "Z^E",
             OutputName -> ZE, 
             LesHouches -> selemix			        }}, 
             
{ZM,        { LaTeX -> "Z^\\mu",
             OutputName -> ZM, 
             LesHouches -> smumix			        }}, 

{ZTau,        { LaTeX -> "Z^\\tau",
             OutputName -> ZTau, 
             LesHouches -> staumix			        }},

{WOp,       {Description -> "Weinberg Operator"}},
{Q[WOp],       {Description -> "Soft Breaking Weinberg Operator"}},

{ZEL,       { Description ->"Left-Lepton-Mixing-Matrix"}},
{ZER,       { Description ->"Right-Lepton-Mixing-Matrix" }},        

{UV,       { Description ->"Neutrino-Mixing-Matrix"}}
             
             
                            

 }; 
 

