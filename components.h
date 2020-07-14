float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

void calculateTau(Matrix &C, mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    zeroes(C,12,3);
    C.at(0).at(0) = 3*(pow(y1,2)+2*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)+3*(2*z1+z2+z3+z4));     
    C.at(1).at(0) = pow(y1,2)+y1*(2*y2+y3+y4)+pow(3*y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)+3*(z1+2*z2+z3+z4);
    C.at(2).at(0) = pow(y1,2)+y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+pow(3*y3,2)+2*y3*y4+pow(y4,2)+3*(z1+z2+2*z3+z4);     
    C.at(3).at(0) = pow(y1,2)+y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*(pow(y4,2)+z1+z2+z3+2*z4);

    C.at(4).at(1) = 3*(pow(y1,2)+2*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)+3*(2*z1+z2+z3+z4));
    C.at(5).at(1) = pow(y1,2)+y1*(2*y2+y3+y4)+pow(3*y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)+3*(z1+2*z2+z3+z4);
    C.at(6).at(1) = pow(y1,2)+y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+pow(3*y3,2)+2*y3*y4+pow(y4,2)+3*(z1+z2+2*z3+z4);     
    C.at(7).at(1) = pow(y1,2)+y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*(pow(y4,2)+z1+z2+z3+2*z4);

    C.at(8).at(2) = 3*(pow(y1,2)+2*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)+3*(2*z1+z2+z3+z4));
    C.at(9).at(2) = pow(y1,2)+y1*(2*y2+y3+y4)+pow(3*y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)+3*(z1+2*z2+z3+z4);
    C.at(10).at(2) = pow(y1,2)+y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+pow(3*y3,2)+2*y3*y4+pow(y4,2)+3*(z1+z2+2*z3+z4);     
    C.at(11).at(2) = pow(y1,2)+y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*(pow(y4,2)+z1+z2+z3+2*z4);

    C.at(4).at(0) = 0;      C.at(0).at(1) = 0;      C.at(0).at(2) = 0;
    C.at(5).at(0) = 0;      C.at(1).at(1) = 0;      C.at(1).at(2) = 0;
    C.at(6).at(0) = 0;      C.at(2).at(1) = 0;      C.at(2).at(2) = 0;
    C.at(7).at(0) = 0;      C.at(3).at(1) = 0;      C.at(3).at(2) = 0;
    C.at(8).at(0) = 0;      C.at(8).at(1) = 0;      C.at(4).at(2) = 0;
    C.at(9).at(0) = 0;      C.at(9).at(1) = 0;      C.at(5).at(2) = 0;
    C.at(10).at(0) = 0;     C.at(10).at(1) = 0;     C.at(6).at(2) = 0;
    C.at(11).at(0) = 0;     C.at(11).at(1) = 0;     C.at(7).at(2) = 0;

}

float calculatePe(mesh m, int i){
     element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    return x1 + x2 + x3 + x4 +z1 +z2 + z3 +z4;
}

void calculatePI(Matrix &C, mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    zeroes(C,12,3);
    C.at(0).at(0) = 3*pow(x1,2)+2*x1*(x2+x3+x4)+pow(x2,2)+x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2);   
    C.at(1).at(0) = pow(x1,2)+x1*(2*x2+x3+x4)+3*pow(x2,2)+2*x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2);
    C.at(2).at(0) = pow(x1,2)+x1*(x2+2*x3+x4)+pow(x2,2)+x2*(2*x3+x4)+3*pow(x3,2)+2*x3*x4+pow(x4,2);    
    C.at(3).at(0) = pow(x1,2)+x1*(x2+x3+2*x4)+pow(x2,2)+x2*(x3+2*x4)+pow(x3,2)+2*x3*x4+3*pow(x4,2);

    C.at(4).at(1) = 3*pow(x1,2)+2*x1*(x2+x3+x4)+pow(x2,2)+x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2);   
    C.at(5).at(1) = pow(x1,2)+x1*(2*x2+x3+x4)+3*pow(x2,2)+2*x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2);
    C.at(6).at(1) = pow(x1,2)+x1*(x2+2*x3+x4)+pow(x2,2)+x2*(2*x3+x4)+3*pow(x3,2)+2*x3*x4+pow(x4,2);     
    C.at(7).at(1) = pow(x1,2)+x1*(x2+x3+2*x4)+pow(x2,2)+x2*(x3+2*x4)+pow(x3,2)+2*x3*x4+3*pow(x4,2);

    C.at(8).at(2) = 3*pow(x1,2)+2*x1*(x2+x3+x4)+pow(x2,2)+x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2);   
    C.at(9).at(2) = pow(x1,2)+x1*(2*x2+x3+x4)+3*pow(x2,2)+2*x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2);
    C.at(10).at(2) = pow(x1,2)+x1*(x2+2*x3+x4)+pow(x2,2)+x2*(2*x3+x4)+3*pow(x3,2)+2*x3*x4+pow(x4,2);    
    C.at(11).at(2) = pow(x1,2)+x1*(x2+x3+2*x4)+pow(x2,2)+x2*(x3+2*x4)+pow(x3,2)+2*x3*x4+3*pow(x4,2);

    C.at(4).at(0) = 0;      C.at(0).at(1) = 0;      C.at(0).at(2) = 0;
    C.at(5).at(0) = 0;      C.at(1).at(1) = 0;      C.at(1).at(2) = 0;
    C.at(6).at(0) = 0;      C.at(2).at(1) = 0;      C.at(2).at(2) = 0;
    C.at(7).at(0) = 0;      C.at(3).at(1) = 0;      C.at(3).at(2) = 0;
    C.at(8).at(0) = 0;      C.at(8).at(1) = 0;      C.at(4).at(2) = 0;
    C.at(9).at(0) = 0;      C.at(9).at(1) = 0;      C.at(5).at(2) = 0;
    C.at(10).at(0) = 0;     C.at(10).at(1) = 0;     C.at(6).at(2) = 0;
    C.at(11).at(0) = 0;     C.at(11).at(1) = 0;     C.at(7).at(2) = 0;


}

void calculateTrara(Matrix &C, mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();
    
    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();
    
    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    zeroes(C,3,12);

    C.at(0).at(0) = 2*x1*(3*z1+z2+z3+z4)+x2*(2*z1+2*z2+z3+z4)+x3*(2*z1+z2+2*z3+z4)+x4*(2*z1+z2+z3+2*z4)+6*(2*y1+y2+y3+y4);
    C.at(0).at(1) = x1*(2*z1+2*z2+z3+z4)+2*x2*(z1+3*z2+z3+z4)+x3*(z1+2*z2+2*z3+z4)+x4*(z1+2*z2+z3+2*z4)+6*(y1+2*y2+y3+y4);
    C.at(0).at(2) = x1*(2*z1+z2+2*z3+z4)+x2*(z1+2*z2+2*z3+z4)+2*x3*(z1+z2+3*z3+z4)+x4*(z1+z2+2*(z3+z4))+6*(y1+y2+2*y3+y4);
    C.at(0).at(3) = x1*(2*z1+z2+z3+2*z4)+x2*(z1+2*z2+z3+2*z4)+x3*(z1+z2+2*(z3+z4))+2*(x4*(z1+z2+z3+3*z4)+3*(y1+y2+y3+2*y4));

    C.at(1).at(4) = 2*x1*(3*z1+z2+z3+z4)+x2*(2*z1+2*z2+z3+z4)+x3*(2*z1+z2+2*z3+z4)+x4*(2*z1+z2+z3+2*z4)+6*(2*y1+y2+y3+y4);
    C.at(1).at(5) = x1*(2*z1+2*z2+z3+z4)+2*x2*(z1+3*z2+z3+z4)+x3*(z1+2*z2+2*z3+z4)+x4*(z1+2*z2+z3+2*z4)+6*(y1+2*y2+y3+y4);
    C.at(1).at(6) = x1*(2*z1+z2+2*z3+z4)+x2*(z1+2*z2+2*z3+z4)+2*x3*(z1+z2+3*z3+z4)+x4*(z1+z2+2*(z3+z4))+6*(y1+y2+2*y3+y4);
    C.at(1).at(7) = x1*(2*z1+z2+z3+2*z4)+x2*(z1+2*z2+z3+2*z4)+x3*(z1+z2+2*(z3+z4))+2*(x4*(z1+z2+z3+3*z4)+3*(y1+y2+y3+2*y4));

    C.at(2).at(8) = 2*x1*(3*z1+z2+z3+z4)+x2*(2*z1+2*z2+z3+z4)+x3*(2*z1+z2+2*z3+z4)+x4*(2*z1+z2+z3+2*z4)+6*(2*y1+y2+y3+y4);
    C.at(2).at(9) = x1*(2*z1+2*z2+z3+z4)+2*x2*(z1+3*z2+z3+z4)+x3*(z1+2*z2+2*z3+z4)+x4*(z1+2*z2+z3+2*z4)+6*(y1+2*y2+y3+y4);
    C.at(2).at(10) = x1*(2*z1+z2+2*z3+z4)+x2*(z1+2*z2+2*z3+z4)+2*x3*(z1+z2+3*z3+z4)+x4*(z1+z2+2*(z3+z4))+6*(y1+y2+2*y3+y4);
    C.at(2).at(11) = x1*(2*z1+z2+z3+2*z4)+x2*(z1+2*z2+z3+2*z4)+x3*(z1+z2+2*(z3+z4))+2*(x4*(z1+z2+z3+3*z4)+3*(y1+y2+y3+2*y4));

    C.at(0).at(4) = 0;
    C.at(0).at(5) = 0;
    C.at(0).at(6) = 0;
    C.at(0).at(7) = 0;
    C.at(0).at(8) = 0;
    C.at(0).at(9) = 0;
    C.at(0).at(10) = 0;
    C.at(0).at(11) = 0;

    C.at(1).at(0) = 0;
    C.at(1).at(1) = 0;
    C.at(1).at(2) = 0;
    C.at(1).at(3) = 0;
    C.at(1).at(8) = 0;
    C.at(1).at(9) = 0;
    C.at(1).at(10) = 0;
    C.at(1).at(11) = 0;

    C.at(2).at(0) = 0;
    C.at(2).at(1) = 0;
    C.at(2).at(2) = 0;
    C.at(2).at(3) = 0;
    C.at(2).at(4) = 0;
    C.at(2).at(5) = 0;
    C.at(2).at(6) = 0;
    C.at(2).at(7) = 0;

}





float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixT,matrixP,matrixL,matrixR,matrixTau,matrixPI,matrixTrara;
    float u_bar,nu,rho,Ve,J,Determinant;

    /* [ T+P  L ]
       [  R   0 ]
    */

    //Matrix T
    Matrix Alpha, Beta;

    
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(360*Determinant);
    
    calculateTau(matrixTau,m,e);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_a, productMatrixMatrix(matrixTau,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixT);


    //Matrix P
    Matrix Alpha_t,Beta_t;
    float pe = calculatePe(m,e);    
    float real_k = (float) (pe)/24*(Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixP);

    
    //Matrix L
    Matrix Omega;
    
    
    float real_g = (float) (J)/(360*Determinant);
    calculatePI(matrixPI,m,e);
    calculateOmega(Omega);
    productRealMatrix(real_g,productMatrixMatrix(matrixPI,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixL);

    //Matrix R
    Matrix Omega_t;
    float real_d = (float)(J/(720*Determinant));

    calculateTrara(matrixTrara,m,e);
    transpose(Omega, Omega_t);
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,matrixTrara,3,3,12),4,3,12),matrixR);


    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixT,matrixP,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixL);
    ubicarSubMatriz(M,12,15,0,11,matrixR);

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

void calculateAtilde(Matrix &C, mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    zeroes(C,4,1);

    C.at(0).at(0) = pow(3*z1,2)+2*z1*(z2+z3+z4)+ pow(z2,2)+z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2);
    C.at(1).at(0) = pow(z1,2)+z1*(2*z2+z3+z4)+3*pow(z2,2)+2*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2);
    C.at(2).at(0) =  pow(z1,2)+z1*(z2+2*z3+z4)+pow(z2,2)+z2*(2*z3+z4)+pow(3*z3,2)+2*z3*z4+pow(z4,2);
    C.at(3).at(0) =  pow(z1,2)+z1*(z2+z3+2*z4)+pow(z2,2)+z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+pow(3*z4,2);

}


Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix, matrixTtilt, respuesta;

    calculateF(f, m);
    calculateAtilde(matrixTtilt,m,e);
    calculateGamma(g_matrix);


    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,16);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    float real = (36*J/360);
    productRealMatrix(real, matrixTtilt,respuesta);

    
    b.at(12) = real;
    b.at(13) = real;
    b.at(14) = real;
    b.at(15) = real;
    
    
    return b;
}
