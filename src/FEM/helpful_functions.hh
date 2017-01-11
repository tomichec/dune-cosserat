
Dune::FieldMatrix<double,6,22> getB(std::vector<Dune::FieldVector<double,2> >& dp2, std::vector<Dune::FieldVector<double,2> >& dp1, std::vector<double> & p1){
     
Dune::FieldMatrix<double,6,22> B(0.0);


        for (int i = 0; i < 9; i++){
            B[0][i] = dp2[i][0]; // E11
            B[1][i + 9] = dp2[i][1]; // E22

            B[2][i] = dp2[i][1];
            B[3][i + 9] = dp2[i][0];

            if (i < 4){
              B[2][i + 18] = p1[i];
              B[3][i + 18] = -p1[i];
              B[4][i + 18] = dp1[i][0];
              B[5][i + 18] = dp1[i][1];
            }
        }


        return B;
}


Dune::FieldVector<double,6> rotateTensor3(Dune::FieldVector<double,6>& e, Dune::FieldVector<double,3> & ROT, int tran){

	// For now only rotate around x2 as all that is required

	Dune::FieldMatrix<double,6,6> R(0.0);

	if (tran == 1){ ROT *= -1.0; }

	double c = std::cos(ROT[1]); double s = std::sin(ROT[1]);

	R[0][0] = c * c; R[0][2] = s * s; R[0][4] = 2 * c * s;
	R[1][1] = 1.0;
	R[2][0] = s * s; R[2][2] = c * c; R[2][4] = -2 * c * s;
	R[3][3] = c; R[3][5] = -s;
	R[4][0] = -c * s; R[4][2] = c * s; R[4][4] = c * c - s * s;
	R[5][3] = s; R[5][5] = c;

	Dune::FieldVector<double,6> e_r;

	for (int i = 0; i < 6; i ++){
		for (int j = 0; j < 6; j++){
			if (R[i][j] < 1e-8){R[i][j] = 0.0;}
		}
	}

	R.mv(e,e_r);

	return e_r;

}

Dune::FieldVector<double,6> rotateTensor3(Dune::FieldVector<double,6>& e, double ROT, int tran){

	// For now only rotate around x2 as all that is required

	Dune::FieldMatrix<double,6,6> R(0.0);

	if (tran == 1){ ROT *= -1.0; }

	double c = std::cos(ROT); double s = std::sin(ROT);

	R[0][0] = c * c; R[0][1] = s * s; R[0][2] = c * s; R[0][3] = c * s;
	R[1][0] = s * s; R[1][1] = c * c; R[1][2] = -c * s; R[1][3] = -c * s;
	R[2][0] = -c * s; R[2][1] = c * s; R[2][2] = c * c; R[2][3] = -s * s;
	R[3][0] = -c * s; R[3][1] = c * s; R[3][2] = -s * s; R[3][3] = c * c;
	R[4][4] = c; R[4][5] = s;
	R[5][4] = -s; R[5][5] = c;

	Dune::FieldVector<double,6> e_r;

	for (int i = 0; i < 6; i ++){
		for (int j = 0; j < 6; j++){
			if (R[i][j] < 1e-8){R[i][j] = 0.0;}
		}
	}

	R.mv(e,e_r);

	return e_r;

}





Dune::FieldVector<double,3> rotateTensor2(Dune::FieldVector<double,3>& e, Dune::FieldVector<double,3> & ROT, int tran){

	// For now only rotate around x2 as all that is required

	Dune::FieldMatrix<double,3,3> T(0.0);

	if (tran == 1){ ROT *= -1.0; }

	double c = std::cos(ROT[1]); double s = std::sin(ROT[1]);

	T[0][0] = c; T[2][2] = c; T[0][2] = s; T[2][0] = -s; T[1][1] = 1.0;

	for (int i = 0; i < 3; i ++){
		for (int j = 0; j < 3; j++){
			if (T[i][j] < 1e-8){T[i][j] = 0.0;}
		}
	}


	Dune::FieldVector<double,3> e_r;

	T.mv(e,e_r);

	return e_r;

}

Dune::FieldVector<double,2> rotateTensor2(Dune::FieldVector<double,2>& e, double ROT, int tran){

	Dune::FieldMatrix<double,2,2> T(0.0);

	if (tran == 1){ ROT *= -1.0; }

	double c = std::cos(ROT); double s = std::sin(ROT);

	T[0][0] = c; T[1][1] = c; T[0][1] = s; T[1][0] = -s;

	for (int i = 0; i < 2; i ++){
		for (int j = 0; j < 2; j++){
			if (T[i][j] < 1e-8){T[i][j] = 0.0;}
		}
	}


	Dune::FieldVector<double,2> e_r;

	T.mv(e,e_r);

	return e_r;

}
