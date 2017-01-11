template <int dim>
class MODEL{

public:

	MODEL(Dune::ParameterTree& config){

		modelName = config.get<std::string>("modelName","myOutput");

		char path[1024], *MakeDir;
 		int length1 = strlen("mkdir ");
		int length2 = strlen(path);
		MakeDir = new char [length1 + length2 + 1];
		strcpy(MakeDir, "mkdir ");
		strcat(MakeDir, modelName.c_str());
		system(MakeDir);

		verb_InitialGuess = config.get<bool>("verb_InitialGuess",true);
		verb_timeIntegration = config.get<bool>("verb_timeIntegration",true);

		adaptive = config.get<bool>("adaptive_timeStepping",true);

		// === Mesh & Geometry Setup

		geometryType = config.get<std::string>("geometryType","donut");

		numLayers = config.get<int>("numLayers", 1);
		ply_th = config.get<double>("plyThickness", 0.22);
		    
		InternalRadius = config.get<double>("InternalRadius",15.0);

		L[0] = M_PI * (2.0 * InternalRadius); // Circumference of Internal Radius
		L[1] = numLayers * ply_th;

		N[0] = config.get<int>("Nx",100);
		N[1] = config.get<int>("Nz",4);

		periodic[0] = config.get<bool>("Periodic_x",true);
		periodic[1] = config.get<bool>("Periodic_z",false);

		overlap = config.get<int>("grid_overlap",1);

		ply_percentage.resize(3);

		ply_percentage[0] = config.get<double>("pp0", 0.25);
		ply_percentage[1] = config.get<double>("pp90", 0.25);
		ply_percentage[2] = config.get<double>("pp45", 0.50);

		a_int = config.get<double>("alpha_interface", 1.0 / 17.0);
		a_ply = config.get<double>("alpha_interface", 16.0 / 17.0);

		// ===

		E1 = config.get<double>("Elastic_E1",120.0);
	   	Gint = config.get<double>("Elastic_Gint",1.0e-4);
	   	Gply = config.get<double>("Elastic_Gply",1.0e-3);

	   	// Consolidation Model
	    a = config.get<double>("Consol_a",1.0);
	    Vf0 = config.get<double>("Consol_Vf0",0.5);
		Vfmax = config.get<double>("Consol_Vfmax",0.68);
	    beta = config.get<double>("Consol_beta",350.0); //
	    Ef  = config.get<double>("Consol_Ef",230.0); //

	    // Permeability Model
	    rf = config.get<double>("Perm_rf",4.0e-3);
	    k1 = config.get<double>("Perm_k1",0.7);
	    k2 = config.get<double>("Perm_k2",0.2);

	    // Cure Model

        Ea = config.get<double>("Cure_Ea",66500); 
    	Aa = config.get<double>("Cure_Aa",1.53e5); 
    	m = config.get<double>("Cure_m",0.813); 
    	n = config.get<double>("Cure_n",2.74); 

    	// Visocisty Model

    	muinfty = config.get<double>("Vis_muinfty",3.45e-19);
        Emu = config.get<double>("Vis_Emu",76536);
        A = config.get<double>("Vis_A",3.8); 
        B = config.get<double>("Vis_B",2.5);
        ag = config.get<double>("Vis_ag",0.47); 

	    // Solver Options

        Newton_AssemblyThreshold = config.get<double>("Newton_AssemblyThreshold",0.0);
	    Newton_Verb_Level = config.get<int>("Newton_Verb_Level",1);
	    Newton_Reduction = config.get<double>("Newton_Reduction",1e-4);
	    Newton_MinLinearReduction = config.get<double>("Newton_MinLinearReduction",2e-5);
        Newton_MaxIterations = config.get<int>("Newton_MaxIterations",125);
        Newton_LineSearchMaxIterations = config.get<int>("Newton_LineSearchMaxIterations",10);

        timeStepper_Verb = config.get<int>("timeStepper_Verb", 1);

        adapt_m = config.get<double>("adaptive_m",2);

        dt = config.get<double>("initial_timeStep",0.1);
        dtmax = config.get<double>("adaptive_maxtimeStep",20.0);

        t1 = config.get<double>("cure_t1",1667.98617611);
		t2 = config.get<double>("cure_t2",1800.00);
		t3 = config.get<double>("cure_t3",1499.99352005);
		t4 = config.get<double>("cure_t4",7200.00);

		tend = t1 + t2 + t3 + t4;

        initial_Alpha = config.get<double>("initial_Alpha",1e-4);
        initial_Temp = config.get<double>("initial_Temp",30.0);








	}

	Dune::FieldVector<double,dim> inline getL(){ return L; }
	Dune::array<int,dim> inline getN(){ return N; }
	std::bitset<dim> inline getPeriodic(){ return periodic; }
	int inline getOverlap(){ return overlap; }
	bool inline getVerb_IG(){return verb_InitialGuess;}
	std::string inline getOutputFolder(){return modelName;}
	int inline getNewtonVerbLevel(){return Newton_Verb_Level; }
	double inline  getNewtonReduction(){return Newton_Reduction; }
	int inline getNewtonMaxIterations(){return Newton_MaxIterations;}
	double inline getNewtonLineReduction(){return Newton_MinLinearReduction;}
	double inline getNewtonAssemblyThreshold(){return Newton_AssemblyThreshold;}
	int inline getNewtonLineSearchMaxIterations(){return Newton_LineSearchMaxIterations;}
	double inline get_Vf0(){return Vf0;}
	bool inline getOSMVerb(){return timeStepper_Verb;}
	double inline get_dt(){return dt;}

	double inline get_initialAlpha(){return initial_Alpha;}
	double inline get_initialTemp(){return initial_Temp;}

	double inline get_tend(){return tend;}

	int inline getVerb_TI(){return verb_timeIntegration;}

	int inline getAdaptM(){return adapt_m;}

	int inline isAdaptive(){return adaptive;}

	double inline get_dtmax(){return dtmax;}


	// Function for geometry transformation of rectangular grid to part geometry (to be defined by the user)
	
	Dune::FieldVector<double,dim> inline evaluateGeometryTransformation(const Dune::FieldVector<double,2>& x){
	
		Dune::FieldVector<double,dim> y(0.0);

		if (geometryType == "donut"){

			double phi = x[0] / InternalRadius;

			y[0] = (InternalRadius + x[1]) * std::sin(phi);
			y[1] = (InternalRadius + x[1]) * std::cos(phi);

		}

		return y;
	}

	// Function to mark boundary nodes and evaluate their prescribed values at the boundary

	bool inline isDirichlet(Dune::FieldVector<double,2>& x, int dof){

		// isDirichlet returns bool if degree of freedom (dof) at node with coords (x) lies on boundary

		bool Dirich = false; // Initialise to default value

		if (dof < 4){

			double distance = std::abs(std::sqrt(x[0] * x[0] + x[1] * x[1]) - InternalRadius);

			if (distance < 1.0e-6){
				Dirich = true;
			}

		}
		else{
			Dirich = true;

		}

		return Dirich; // return bool

	} // end isDirichlet

	double inline evaluateDirichlet(const Dune::FieldVector<double,2>& x, int dof,double time){
		
		// Returns the value at a boundary node (at x) for degree of freedom dof

		double value = 0.0;

		return value;
	}


	// Tensors

	Dune::FieldMatrix<double,6,6> inline getElasticTensor(double Vf, double ROT){

		Dune::FieldMatrix<double,6,6> C(0.0);
		
		C[0][0] = (ply_percentage[0]  + 0.25 * ply_percentage[2]) * E1; // Rule of Mixtures

		double As = (3.0 * M_PI * Ef / std::pow(beta,4));
  		double term1 = 1.0 / Vf - 1.0 / Vfmax;
    	C[1][1] = -As * ( (-Vf / Vf0) / (std::pow(term1,4)) + 4.0 * ((1.0 - Vf / Vf0) / (std::pow(term1,5)))) * (1 / (Vf * Vf));

		C[0][1] = 0.0; C[1][0] = 0.0; // Fiber bed - Poisson Effects Neglected

		C[2][2] = a_int / Gint + a_ply / Gply;	C[2][2] = 1.0 / C[2][2];

		C[3][3] = a_int * Gint + a_ply * Gply;

		double tply = a_ply * ply_th;

		C[4][4] = C[0][0] * tply * tply * (C[3][3] - C[2][2])/(12 * C[2][2]);

		return C;
	}

	Dune::FieldMatrix<double,2,2> inline getPermTensor(double Vf,double mu, double ROT){

		Dune::FieldMatrix<double,2,2> K(0.0);

		K[0][0] = (1.0 / (4.0 * k1)) * ( rf * rf ) * std::pow(1 - Vf,3) / (Vf * Vf);

		double tmp_d = std::sqrt(Vfmax / Vf);

		K[1][1] = ((rf * rf) / (4.0 * k2)) * std::pow(tmp_d - 1.0,3) / (tmp_d + 1.0);

		K /= mu;

		return K;

	}

	double inline getVisocity(double Temp, double degree_of_cure){
	/*
		 * Viscosity Model :- mu = muinfty * exp(Emu / R * Temp) * (ag / (ag - degree_of_cure))^(A + B * degree_of_cure)
		 * ------
		 * Kenny, J.M., “Integration of Process Models with Control and Optimization of Polymer
		 * Composites Fabrication”, Proceedings of the Third Conference on Computer Aided Design
		 * in Composites Materials Technology, 1992, pp. 530-544.
		 * ----
		 * Comments: Good Model for 8552 Model tends to predict a lower minimum visocity at rates greater than 5K/min,
		 * this not a particular problem as autoclave heating rates rarely exceed 4oC
	*/

		double mu = muinfty * std::exp (Emu / (Gas_Constant * Temp) ) * std::pow(ag / (ag - degree_of_cure),A + B * degree_of_cure); //

		return mu;

	}

	void inline cureModel(double a0, double T0, double dt){

	// solve the cure model using Rugge-Kutta method

	double c1 = g(T0, a0);

	double c2 = g(T0 + 0.5 * dt * temperature_rampRate, a0 + 0.5 * dt * c1);

	double c3 = g(T0 + 0.5 * dt * temperature_rampRate, a0 + 0.5 * dt * c2);

	double c4 = g(T0 + dt * temperature_rampRate, a0 + dt * c3);

	alpha_temp = a0 + (dt / 6.0) * (c1 + 2.0 * c2 + 2.0 * c3 + c4);

}

double inline g(double T,double a){
	double da = Aa * std::exp(-Ea / (Gas_Constant * T)) * std::pow(a,m) * std::pow (1-a,n);
	return da;
}


	double inline get_new_temp(double tt,double Dt){

		double t = tt + Dt;

		if (t < t1){
			temp_temp = initial_Temp  + t * (2.7778 / 60.0);
		}
		else if(t < t2 + t1){
			temp_temp = 107.2222;
		}
		else if (t < t3 + t2 + t1){
			temp_temp = 107.2222 + (t - (t2 + t1)) * (2.7778 / 60.0);

		}
		else{
			temp_temp = 176.6667;
		}

		temp_temp += absolute_zero; // Convert to Kelvins



		return temp_temp;
	}


	// Setting Values

	void inline setTime(double input){time = input;}
	void inline setPressure(double input){pressure = input;}

	void inline set_alpha(double vv){alpha.push_back(vv);}
	double inline get_alpha(){return alpha[current_step];}

	void inline set_alpha_temp(double a0){	alpha_temp = a0; }
	double inline get_alpha_temp(){ return alpha_temp;}

	void inline setStep(int input){current_step = input;}

	void inline set_temp(double vv){temp.push_back(vv);}

	void inline setCure(){alpha.push_back(alpha_temp); }

	void inline recordTime(double currentTime){timeVector.push_back(currentTime);}





	





private:

	int current_step;

	double time, pressure;

	double dt, tend, t1,t2,t3,t4 ;

	std::string modelName;

	bool verb_InitialGuess, verb_timeIntegration;

	// Solver Parameters

	int Newton_Verb_Level, Newton_MaxIterations, Newton_LineSearchMaxIterations ;

	double Newton_Reduction, Newton_MinLinearReduction, Newton_AssemblyThreshold;

	int timeStepper_Verb;

	int adapt_m;

	// Mesh & Geometry Parameters

	std::string geometryType;

	Dune::FieldVector<double,dim> L;
	std::bitset<dim> periodic;
	Dune::array<int,dim> N;

	int overlap;

	int numLayers;

	double ply_th;

	double InternalRadius;

	// Composite Parameters

	std::vector<double> ply_percentage;

	double E1, Gint, Gply, a_ply, a_int;

	double a, Vf0, Vfmax, beta, Ef;
    

    double Gas_Constant = 8.3144598; // J/ (mol * K)
    double absolute_zero = 273.15;
	double atmos = -1.01324e-4;

	double rf, k1, k2; // Permebaility Model

	double alpha_temp, Ea, Aa, m, n; // Cure Model

	double muinfty, Emu, A , B, ag;
        
	std::vector<double> alpha, temp, timeVector; // Vector containing Cure State and Temperature

	double temp_temp;

	double initial_Alpha,initial_Temp;

	double temperature_rampRate = 0.0;

	bool adaptive;

	double dtmax;


};
