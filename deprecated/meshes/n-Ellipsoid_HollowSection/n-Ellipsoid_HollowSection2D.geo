Include "../Parameters.geo";
//Geom_AR = GEOM_AR_2; MeshLevel = 0; MeshType = TRI; MeshCurving = CURVED; Extended = EXTENSION_DISABLED; PDEName = POISSON;

// Modifiable Parameters
AddRadial = 0;
SwapTriOrientation = 0;



// Geometry Specification
rIn  = 0.5;
rOut = 1.0;

aIn  = 0.5;
aOut = 1.0;

bIn  = 0.5;

If (MeshCurving == CURVED)
	If (Geom_AR == GEOM_AR_3)
		bOut = 3.0*aOut;
	ElseIf (Geom_AR == GEOM_AR_2)
		bOut = 2.0*aOut;
	ElseIf (Geom_AR == GEOM_AR_1)
		bOut = 1.0*aOut;
	EndIf

	t   = Pi/4.0;
	r   = rIn;
	r_e = Sqrt(1.0/((Cos(t)/aOut)^2.0+(Sin(t)/bOut)^2.0));

	Point(1) = {aIn,0,0,lc};
	Point(2) = {aOut,0,0,lc};
	Point(3) = {0,bIn,0,lc};
	Point(4) = {0,bOut,0,lc};
	Point(5) = {r*Cos(t),r*Sin(t),0,lc};
	Point(6) = {r_e*Cos(t),r_e*Sin(t),0,lc};
	Point(7) = {0,0,0,lc};

	Line(1001)    = {1,2};
	Line(1002)    = {3,4};
	Line(1003)    = {5,6};
	Ellipse(1004) = {5,7,3,1}; // start, centre, point on major-axis, end
	Ellipse(1005) = {5,7,3,3};
	Ellipse(1006) = {6,7,4,2};
	Ellipse(1007) = {6,7,4,4};
ElseIf (MeshCurving == TOBECURVED)
	bOut = aOut;

	Point(1) = {aIn,0,0,lc};
	Point(2) = {aOut,0,0,lc};
	Point(3) = {0,bIn,0,lc};
	Point(4) = {0,bOut,0,lc};
	Point(5) = {aIn,bIn,0,lc};
	Point(6) = {aOut,bOut,0,lc};
	Point(7) = {0,0,0,lc};

	Line(1001) = {1,2};
	Line(1002) = {3,4};
	Line(1003) = {5,6};
	Line(1004) = {5,1};
	Line(1005) = {5,3};
	Line(1006) = {6,2};
	Line(1007) = {6,4};
EndIf

Transfinite Line {1001:1003} = 1*2^(MeshLevel+AddRadial)+1 Using Progression 1;
Transfinite Line {1004:1007} = 1*2^(MeshLevel)+1 Using Progression 1;

Line Loop (4001) = {1001,-1006,-1003,1004};
Line Loop (4002) = {-1005,1003,1007,-1002};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

If (MeshCurving == CURVED)
	Transfinite Surface{4001} Right;
	Transfinite Surface{4002};
ElseIf (MeshCurving == TOBECURVED)
	If (SwapTriOrientation)
		Transfinite Surface{4001} Right;
		Transfinite Surface{4002};
	Else
		Transfinite Surface{4001};
		Transfinite Surface{4002} Right;
	EndIf
EndIf


If (MeshType == QUAD)
	Recombine Surface{4001,4002};
ElseIf (MeshType == MIXED2D)
	Recombine Surface{4001};
EndIf



// Physical Parameters for '.msh' file
BC_Straight =   BC_STEP_SC;
BC_Curved   = 2*BC_STEP_SC;

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};

If (Extended == EXTENSION_ENABLED)
	l = 2*rIn;

	Point(8)  = {aIn,-l,0,lc};
	Point(9)  = {aOut,-l,0,lc};

	Point(10) = {-aIn,0,0,lc};
	Point(11) = {-aOut,0,0,lc};
	Point(12) = {-r*Cos(t),r*Sin(t),0,lc};
	Point(13) = {-r_e*Cos(t),r_e*Sin(t),0,lc};
	Point(14)  = {-aIn,-l,0,lc};
	Point(15)  = {-aOut,-l,0,lc};

	Line(1008) = {8,9};
	Line(1009) = {1,8};
	Line(1010) = {2,9};


	Line(1011)    = {10,11};
	Line(1012)    = {12,13};
	Ellipse(1013) = {12,7,3,10};
	Ellipse(1014) = {12,7,3,3};
	Ellipse(1015) = {13,7,4,11};
	Ellipse(1016) = {13,7,4,4};

	Line(1017) = {14,15};
	Line(1018) = {10,14};
	Line(1019) = {11,15};

	Transfinite Line {1008,1011,1012,1017} = 1*2^(MeshLevel+AddRadial)+1 Using Progression 1;
	Transfinite Line {1013:1016}           = 1*2^(MeshLevel)+1 Using Progression 1;
	Transfinite Line {1009,1010,1018,1019} = 1*2^(MeshLevel)+1 Using Progression 1;

	Line Loop (4003) = {-1011,-1013,1012,1015};
	Line Loop (4004) = {1014,1002,-1016,-1012};
	Line Loop (4005) = {1008,-1010,-1001,1009};
	Line Loop (4006) = {-1017,-1018,1011,1019};

	Plane Surface(4003) = {4003};
	Plane Surface(4004) = {4004};
	Plane Surface(4005) = {4005};
	Plane Surface(4006) = {4006};

	Transfinite Surface{4004};
	Transfinite Surface{4003,4005,4006} Right;

	If (MeshType == QUAD)
		Recombine Surface{4003:4006};
	ElseIf (MeshType == MIXED2D)
		Recombine Surface{4002,4005};
	EndIf

	Physical Surface(9403) = {4003};
	Physical Surface(9404) = {4004};
	Physical Surface(9405) = {4005};
	Physical Surface(9406) = {4006};

	// Physical Parameters for '.msh' file
	If (PDEName == POISSON)
		Physical Line(10011) = {1008,1017};           // Straight Dirichlet
		Physical Line(10012) = {1009:1010,1018:1019}; // Straight Neumann
		Physical Line(20012) = {1004:1007,1013:1016}; // Curved   Neumann
	ElseIf (PDEName == EULER) // Euler
		Physical Line(10004) = {1017};                // Straight Total Temperature/Pressure
		Physical Line(10003) = {1008};                // Straight Back Pressure
		Physical Line(10002) = {1009:1010,1018:1019}; // Straight SlipWall
		Physical Line(20002) = {1004:1007,1013:1016}; // Curved   SlipWall
	EndIf


	// Visualization in gmsh

	Color Black{ Surface{4001:4006}; }
	Geometry.Color.Points = Black;
Else
	// Physical Parameters for '.msh' file
	If (PDEName == POISSON) // Poisson
		Physical Line(10011) = {1002}; // Straight Dirichlet
		Physical Line(10012) = {1001}; // Straight Neumann
		Physical Line(20011) = {1004:1005}; // Curved Dirichlet
		Physical Line(20012) = {1006:1007}; // Curved Neumann
	ElseIf (PDEName == EULER) // Euler
		Physical Line(10006) = {1002};                // Straight Supersonic Outflow
		Physical Line(10005) = {1001};                // Straight Supersonic Inflow
		Physical Line(20002) = {1004:1007};           // Curved   SlipWall
	EndIf

	// Visualization in gmsh

	Color Black{ Surface{4001:4002}; }
	Geometry.Color.Points = Black;
EndIf
