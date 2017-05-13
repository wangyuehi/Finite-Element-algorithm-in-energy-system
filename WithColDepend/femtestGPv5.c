
#include "math.h"
#include "stdio.h"
#include <stdlib.h>
#include "femtestGPv1.h"
#include <windows.h>  
#include <time.h> //time_t time()  clock_t clock()  
#include <Mmsystem.h>             //timeGetTime()  
#pragma comment(lib, "Winmm.lib")   //timeGetTime()  

static double pi = 3.1416;

/*initial info*/

static int NNO = 534;      // number of nodes---- - Number of nodes changes in each time step due to supplementary nodes but their contribution will be assembled to other nodes so, the matrix size will not be affects by the supplementary nodes
static double NEL = 1008;  // number of elements---- - which are fixed in moving structures
static int ap = 0;         // antiperiodicity index for reversqrtg the current direction in each position reset
static double vel = 26.82; // speed(m / s)
static double L = 0.2667;  // machine depth(m)---- - used in force calculation section only

//moving band upper nodes
int mbu[49] = { 24, 25, 26, 27, 28, 29, 30, 31, 110, 111
, 112, 113, 114, 115, 116, 117, 196, 197, 198, 199
, 200, 201, 202, 203, 282, 283, 284, 285, 286, 287
, 288, 289, 368, 369, 370, 371, 372, 373, 374, 375
, 454, 455, 456, 457, 458, 459, 460, 461, 521
};

//moving band downer nodes
int mbd[49] = { 16, 17, 18, 19, 20, 21, 22, 23, 102, 103, 104, 105, 106, 107, 108, 109, 188, 189, 190, 191, 192, 193, 194, 195, 274, 275, 276, 277, 278, 279, 280, 281, 360, 361, 362, 363, 364, 365, 366, 367, 446, 447, 448, 449, 450, 451, 452, 453, 520 };

//moving band meshes
int remesh[96] = { 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165
, 166, 167, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331
, 332, 333, 334, 335, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497
, 498, 499, 500, 501, 502, 503, 656, 657, 658, 659, 660, 661, 662, 663
, 664, 665, 666, 667, 668, 669, 670, 671, 824, 825, 826, 827, 828, 829
, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 992, 993, 994, 995
, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007 };
/*{ 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008 };*/

int n_supp = 0;
int N_m = 49;            // indicates the number of moving band nodes on upper and lower sides
/*mesh*/
int  mesh[1008][3] = { 0 };
double X[600] = { 0 };
double Y[600] = { 0 };

/*time dependency modeling*/

double dt = 1e-3;                       //time step
int dtrev = 1000;
double T = 0.01; //simulation time
double A[600][1000] = { 0 };    //the solution matrix
double At[600] = { 0 }; //temporary At used to expand the result before saving to A
double F_thrust[1000] = { 0 };
double F_Normal[1000] = { 0 };

/*meshes speed and conductivity*/
double v[1008] = { 0 };
double Sigma[1008] = { 0 };

double J[1008] = { 0 };    //Current density, only available in primary mover windings


double  nu[1008] = { 0 };
double dnoodB2[1008] = { 0 };

double An[600][20] = { 0 };

int N1, N2, N3;

double DET;

double a2, a3;

double B[1200] = { 0 };

double SS[700][700] = { 0 }; //stifness matrix
double SJ[700][700] = { 0 }; // Jakoobian Matrix
double b[800] = { 0 };    //source term

int	L_P[70];
int R_P[70];

int I;
int res;

double res1[600] = { 0 };

double den[516] = { 0 };
int ns = 516;
double res2 = 0.0;
// SS, SJ and b are now 516 * 516
double frac[600] = { 0 };
/* Naive solver code*/
double num[530][530] = { 0 };
double res3[600] = { 0 };
double res4[600] = { 0 };
int levelOf[516] = { 0 };
int ancester[516] = { 0 };



void Initialize(){

	int i, j;


	/*initial solution*/
Initloop1:for (i = 0; i < 600; i++){
	A[i][0] = 0.001;
	A[i][1] = 0.f;
	A[i][2] = 0.f;
}

		  for (i = 0; i < 600; i++){
		  Initloop2:for (j = 1; j < 20; j++){
			  A[i][j] = 0.0;
		  }
		  }


}
void MeshSpeedAndConductivity(){
	/*meshes speed and conductivity*/
	int i;


MeshSpeedAndConductivity_loop1:for (i = 51; i < 992; i++){
	if (i == 152)
		i = 219;
	else if (i == 320)
		i = 387;
	else if (i == 488)
		i = 555;
	else if (i == 656)
		i = 723;
	else if (i == 824)
		i = 891;

	v[i] = vel;
}



						   MeshSpeedAndConductivity_loop2:for (i = 10; i < 866; i++){
							   if (i == 26)
								   i = 178;
							   else if (i == 194)
								   i = 346;
							   else if (i == 362)
								   i = 514;
							   else if (i == 530)
								   i = 682;
							   else if (i == 698)
								   i = 850;

							   Sigma[i] = 32500000;
						   }


}
void MovingBandModel(){
	int i, j;
	/*skip the gaps between nodes*/
MovingBandModel_label0:for (i = 24; i < 534; i++){



	if (i == 86)
		i = 110;
	if (i == 172)
		i = 196;
	if (i == 258)
		i = 282;
	if (i == 344)
		i = 368;
	if (i == 430)
		i = 454;
	if (i == 516)
		i = 521;

	X[i] = X[i] + vel*dt;
}
					   /* restarting the position of the rotor when it reaches the end of stator */

					   if (X[24] > 0.1524){
						   ap = ap + 1; // the signal for reversing the current direction after each position reset to satisfy the antiperiodicity properties
					   MovingBandModel_loop2:for (i = 24; i < 534; i++){

						   /*skip the gaps between nodes*/

						   if (i == 86)
							   i = 110;
						   else if (i == 172)
							   i = 196;
						   else if (i == 258)
							   i = 282;
						   else if (i == 344)
							   i = 368;
						   else if (i == 430)
							   i = 454;
						   else if (i == 516)
							   i = 521;

						   X[i] = X[i] - 0.1524;
					   }
					   }


					   /* finding the appropriate elements in the airgap(excluded the supplementary nodes) */

					   n_supp = 0;
					   N_m = 49;

				   MovingBandModel_loop3:for (i = 0; i < 49 - 1; i++){

					   if (X[24] >= X[mbd[i]] && X[24] < X[mbd[i + 1]]){
						   //	printf("%d, Xmbdi = %f, X[mbdi+1] = %f\n", mesh[165][1], X[mbd[i]], X[mbd[i+1]]);
					   MovingBandModel_innerLoop0:	for (j = 0; j < N_m - n_supp - 1; j++){



						   mesh[remesh[2 * j]][0] = mbu[j];
						   mesh[remesh[2 * j]][1] = mbd[j + n_supp + 1];
						   mesh[remesh[2 * j]][2] = mbd[j + n_supp];

						   mesh[remesh[2 * j + 1]][0] = mbu[j];
						   mesh[remesh[2 * j + 1]][1] = mbu[j + 1];
						   mesh[remesh[2 * j + 1]][2] = mbd[j + n_supp + 1];

					   }
													break;
					   }// as soon as finding appropriate mesh the for loop should be ended
					   else{
						   n_supp = n_supp + 1; // indicates the number of supplementary nodes
					   }

				   }
										 //printf("n_supp = %d", n_supp);
										 /* defining the supplementary nodes positions using out nodes due to movement */
									 MovingBandModel_OutLoop1:	for (i = N_m - n_supp; i < N_m; i++){
										 //printf("i = %d,i + NNO + n_supp - N_m = %d\n", i, i + NNO + n_supp - N_m);
										 X[i + NNO + n_supp - N_m] = X[mbu[i - 1]] - 0.1524;
										 Y[i + NNO + n_supp - N_m] = 0.0318;

									 }


																/* defining the new elements related to the supplementary nodes */
																if (n_supp > 0){

																MovingBandModel_OutLoop2:	for (int i = 0; i < n_supp - 1; i++){
																	mesh[remesh[2 * (i + N_m - n_supp - 1)]][0] = NNO + i;
																	mesh[remesh[2 * (i + N_m - n_supp - 1)]][1] = mbd[i + 1];
																	mesh[remesh[2 * (i + N_m - n_supp - 1)]][2] = mbd[i];


																	mesh[remesh[2 * (i + N_m - n_supp - 1) + 1]][0] = NNO + i;
																	mesh[remesh[2 * (i + N_m - n_supp - 1) + 1]][1] = NNO + i + 1;
																	mesh[remesh[2 * (i + N_m - n_supp - 1) + 1]][2] = mbd[i + 1];
																}

																							/*should this be subtracted by 1 */
																							mesh[remesh[2 * (N_m - 1) - 2]][0] = NNO + n_supp - 1;
																							mesh[remesh[2 * (N_m - 1) - 2]][1] = mbd[n_supp];
																							mesh[remesh[2 * (N_m - 1) - 2]][2] = mbd[n_supp - 1];

																							mesh[remesh[2 * (N_m - 1) - 1]][0] = NNO + n_supp - 1;
																							mesh[remesh[2 * (N_m - 1) - 1]][1] = mbu[0];
																							mesh[remesh[2 * (N_m - 1) - 1]][2] = mbd[n_supp];

																}

}
int isEven(int exp)
{
	if (exp % 2 == 1)
		return -1;
	else
		return 1;


}
void SpecifyingSource(double t){
	/*------------------------Specifying the Source---------------------------- -*/
	int i;

	double Acu = 1.122e-3 / 2;                      // copper cross section area of a slot divided by 2 for each winding(m2)
	int f = 92;                                //supply frequency(Hz)

	double Amp = sqrt(2) * 3 * 1541 / (4 * Acu);          //supply current density amplitude(A)  J = NI / (c*A) where c is the number of parallel circuits

	double Ja = Amp*sinf(2 * pi*f*(t + dt));           //Current density of phase a
	double Jb = Amp*sinf(2 * pi*f*(t + dt) - 4 * pi / 3);
	double Jc = Amp*sinf(2 * pi*f*(t + dt) - 2 * pi / 3);


	/* for considering the reverse current when the mover reaches the end and reset to the begin*/
SpecifyingSource_loop1:for (i = 0; i < NEL; i++){
	J[i] = 0;
}
				   SpecifyingSource_loop2:for (i = 113; i < 320; i++){
					   if (i == 152)
						   i = 304;
					   J[i] = isEven(ap)*Ja;                // first slot up down + second slot up

				   }

									  SpecifyingSource_loop3:	for (i = 953; i < 976; i++){
										  J[i] = -1 * isEven(ap)*Ja;
									  }


															SpecifyingSource_loop4:for (i = 617; i < 992; i++){
																if (i == 640)
																	i = 785;
																else if (i == 824)
																	i = 976;

																J[i] = isEven(ap)*Jb;

															}

																			   SpecifyingSource_loop5:for (i = 281; i < 656; i++){
																				   if (i == 304)
																					   i = 449;
																				   else if (i == 488)
																					   i = 640;

																				   J[i] = -1 * isEven(ap)*Jc;

																			   }


}
void SpecifyingMaterials(int I){
	/*indicates the material permieability of each element, 1 for air, copper and the aluminium sheet regions*/
	int i, j;
	double detrev, pirev;
	pirev = 1 / (4 * pi*1e-7);
SpecifyingMaterials_loop1:
	for (i = 0; i < 1008; i++){
		nu[i] = pirev;
	}
SpecifyingMaterials_loop2:for (i = 0; i < 1008; i++){
	dnoodB2[i] = 0;
}
					  SpecifyingMaterials_loop3:for (i = 0; i < 850; i++)

					  {

						  if (i == 10)
							  i = 168;
						  else if (i == 178)
							  i = 336;
						  else if (i == 346)
							  i = 504;
						  else if (i == 514)
							  i = 672;
						  else if (i == 682)
							  i = 840;

						  N1 = mesh[i][0];
						  N2 = mesh[i][1];
						  N3 = mesh[i][2];

						  DET = X[N2] * Y[N3] + X[N1] * Y[N2] + X[N3] * Y[N1] - X[N1] * Y[N3] - X[N3] * Y[N2] - X[N2] * Y[N1];
						  detrev = 1 / (DET);
						  //printf("DET = %f\n", DET);
						  a2 = (An[N2][I] * Y[N3] - An[N3][I] * Y[N2] - An[N1][I] * (Y[N3] - Y[N2]) + Y[N1] * (An[N3][I] - An[N2][I])) *detrev;
						  a3 = (An[N3][I] * X[N2] - An[N2][I] * X[N3] - X[N1] * (An[N3][I] - An[N2][I]) + An[N1][I] * (X[N3] - X[N2])) *detrev;

						  B[i] = sqrt(a2*a2 + a3 *a3);

						  //preventing B from having zero value
						  if (B[i] == 0){
							  B[i] = 1e-12;
						  }
						  //magnetic saturation curve "1010 steel"
						  if (B[i] < 0.6){
							  nu[i] = 135.28;
							  dnoodB2[i] = 0;
						  }
						  else{
							  nu[i] = 79.57*((B[i] - 0.6) + 1.7);
							  dnoodB2[i] = 79.57;
						  }


					  }

											SpecifyingMaterials_loop4:for (i = 51; i < 939; i++){

												if (i == 99)
													i = 219;
												else if (i == 267)
													i = 387;
												else if (i == 435)
													i = 555;
												else if (i == 603)
													i = 723;
												else if (i == 771)
													i = 891;


												N1 = mesh[i][0];                                                       // Global node number of that element
												N2 = mesh[i][1];
												N3 = mesh[i][2];

												DET = X[N2] * Y[N3] + X[N1] * Y[N2] + X[N3] * Y[N1] - X[N1] * Y[N3] - X[N3] * Y[N2] - X[N2] * Y[N1];

												detrev = 1 / (DET);

												a2 = (An[N2][I] * Y[N3] - An[N3][I] * Y[N2] - An[N1][I] * (Y[N3] - Y[N2]) + Y[N1] * (An[N3][I] - An[N2][I])) *detrev;
												a3 = (An[N3][I] * X[N2] - An[N2][I] * X[N3] - X[N1] * (An[N3][I] - An[N2][I]) + An[N1][I] * (X[N3] - X[N2])) *detrev;
												B[i] = sqrt(a2 * a2 + a3 * a3);

												/*preventing B from having zero value*/
												if (B[i] == 0)
													B[i] = 1e-12;

												/*magnetic saturation curve "M19"*/
												if (B[i] < 0.6){
													nu[i] = 135.28;
													dnoodB2[i] = 0;
												}
												else{
													nu[i] = 79.57*(((B[i] - 0.6)*(B[i] - 0.6)) + 1.7);
													dnoodB2[i] = 79.57;
												}




											}


}
void ComputeS1(double S1[3][3], double S3[3][3], double COEFF1, double COEFF3, double Q1, double Q2, double Q3, double R1, double R2, double R3){


	S1[0][0] = COEFF1*(Q1*Q1 + R1*R1);
	S1[0][1] = COEFF1*(Q1*Q2 + R1*R2);
	S1[0][2] = COEFF1*(Q1*Q3 + R1*R3);
	S1[1][0] = S1[0][1];
	S1[1][1] = COEFF1*(Q2*Q2 + R2*R2);
	S1[1][2] = COEFF1*(Q2*Q3 + R2*R3);
	S1[2][0] = S1[0][2];
	S1[2][1] = S1[1][2];
	S1[2][2] = COEFF1*(Q3*Q3 + R3*R3);

	S3[0][0] = COEFF3*Q1;
	S3[0][1] = COEFF3*Q2;
	S3[0][2] = COEFF3*Q3;
	S3[1][0] = S3[0][0];
	S3[1][1] = S3[0][1];
	S3[1][2] = S3[0][2];
	S3[2][0] = S3[0][0];
	S3[2][1] = S3[0][1];
	S3[2][2] = S3[0][2];

}
void ComputeS2(double S2[3][3], double COEFF2){


	S2[0][0] = COEFF2;
	S2[0][1] = COEFF2*0.5;
	S2[0][2] = COEFF2*0.5;
	S2[1][0] = S2[0][1];
	S2[1][1] = COEFF2;
	S2[1][2] = COEFF2*0.5;
	S2[2][0] = S2[0][2];
	S2[2][1] = S2[1][2];
	S2[2][2] = COEFF2;


}
void ComputeS(double S1[3][3], double S2[3][3], double S3[3][3], double COEFF1, double COEFF2, double COEFF3, double Q1, double Q2, double Q3, double R1, double R2, double R3){
	ComputeS1(S1, S3, COEFF1, COEFF3, Q1, Q2, Q3, R1, R2, R3);
	ComputeS2(S2, COEFF2);

}
void MatrixEquation(double t, int I){
	int i, j, i1, i2, k, n;
	double Q1, Q2, Q3, R1, R2, R3, COEFF1, COEFF2, COEFF3, COEFF4;
	double S1[3][3];
	double S2[3][3];
	double J1[3][3];
	double S3[3][3];
	double b1[3];
	double b2[3];
	double tmp1, tmp2, tmp3;


	for (i = 0; i<3; i++){
	MatrixEquation_loopi1:for (j = 0; j<3; j++){
		S1[i][j] = 0;
		S2[i][j] = 0;
		S3[i][j] = 0;
		J1[i][j] = 0;
	}
						  b1[i] = 0;
						  b2[i] = 0;
	}
	for (i = 0; i < 700; i++)
	{
	MatrixEquation_loop1_1:for (j = 0; j < 700; j++){
		SS[i][j] = 0;
		SJ[i][j] = 0;
	}

	}
MatrixEquation_loopinitB:for (i = 0; i < 800; i++){
	b[i] = 0;
}
						 /*---------------------- - stiffness matrix formation------------------------ -*/
						 //NEL
						 for (i = 0; i < 1008; i++){



							 N1 = mesh[i][0]; // Global node number of that element
							 N2 = mesh[i][1];
							 N3 = mesh[i][2];


							 Q1 = Y[N2] - Y[N3];
							 Q2 = Y[N3] - Y[N1];
							 Q3 = Y[N1] - Y[N2];
							 R1 = X[N3] - X[N2];
							 R2 = X[N1] - X[N3];
							 R3 = X[N2] - X[N1];


							 tmp1 = X[N2] - X[N1];
							 tmp1 *= Y[N3];

							 tmp2 = X[N1] - X[N3];
							 tmp2 *= Y[N2];

							 tmp3 = X[N3] - X[N2];
							 tmp3 *= Y[N1];
							 DET = tmp1 + tmp2 + tmp3;

							 COEFF1 = nu[i] * 1 / ((DET * 2));
							 COEFF2 = Sigma[i] * DET * 1 / ((12 * dt));
							 COEFF3 = Sigma[i] * v[i] * 0.166667;
							 COEFF4 = J[i] * DET *0.166667;
							 /*if (i > 987 && i < 991){
							 printf("i %d,COE4 %f,Ji %f,DET %f\n", i, COEFF4, J[i], DET);
							 }*/
							 /*----stiffness matrix for each element corresponding to the main term------*/

							 ComputeS(S1, S2, S3, COEFF1, COEFF2, COEFF3, Q1, Q2, Q3, R1, R2, R3);

							 /*S1[0][0] = COEFF1*(Q1*Q1 + R1*R1);
							 S1[0][1] = COEFF1*(Q1*Q2 + R1*R2);
							 S1[0][2] = COEFF1*(Q1*Q3 + R1*R3);
							 S1[1][0] = S1[0][1];
							 S1[1][1] = COEFF1*(Q2*Q2 + R2*R2);
							 S1[1][2] = COEFF1*(Q2*Q3 + R2*R3);
							 S1[2][0] = S1[0][2];
							 S1[2][1] = S1[1][2];
							 S1[2][2] = COEFF1*(Q3*Q3 + R3*R3);*/



							 for (k = 0; k < 3; k++){
							 MatrixEquation_loop2_1_1:for (n = 0; n < 3; n++){
								 J1[k][n] = S1[k][n] + 4 / (nu[i] * nu[i] * DET) *(dnoodB2[i]) * (S1[k][0] * An[N1][I] + S1[k][1] * An[N2][I] + S1[k][2] * An[N3][I]) *
									 (S1[n][0] * An[N1][I] + S1[n][1] * An[N2][I] + S1[n][2] * An[N3][I]);

							 }
							 }

							 /*--stiffness matrix for each element corresponding to the conducting part--*/



							 /*S2[0][0] = COEFF2;
							 S2[0][1] = COEFF2*0.5;
							 S2[0][2] = COEFF2*0.5;
							 S2[1][0] = S2[0][1];
							 S2[1][1] = COEFF2;
							 S2[1][2] = COEFF2*0.5;
							 S2[2][0] = S2[0][2];
							 S2[2][1] = S2[1][2];
							 S2[2][2] = COEFF2;*/

							 /*----stiffness matrix for each element corresponding to the moving part----*/



							 /*	S3[0][0] = COEFF3*Q1;
							 S3[0][1] = COEFF3*Q2;
							 S3[0][2] = COEFF3*Q3;
							 S3[1][0] = S3[0][0];
							 S3[1][1] = S3[0][1];
							 S3[1][2] = S3[0][2];
							 S3[2][0] = S3[0][0];
							 S3[2][1] = S3[0][1];
							 S3[2][2] = S3[0][2];*/

							 /*------ - source term corresponding to the eddy current current source------ -*/

							 b1[0] = COEFF2*(A[N1][(int)(t / dt + 1) - 1] + 0.5*A[N2][(int)(t / dt + 1) - 1] + 0.5*A[N3][(int)(t / dt + 1) - 1]);
							 b1[1] = COEFF2*(0.5*A[N1][(int)(t / dt + 1) - 1] + A[N2][(int)(t / dt + 1) - 1] + 0.5*A[N3][(int)(t / dt + 1) - 1]);
							 b1[2] = COEFF2*(0.5*A[N1][(int)(t / dt + 1) - 1] + 0.5*A[N2][(int)(t / dt + 1) - 1] + A[N3][(int)(t / dt + 1) - 1]);

							 /*-------- - source term corresponding to the external current source-------- -*/


							 b2[0] = COEFF4;
							 b2[1] = COEFF4;
							 b2[2] = COEFF4;


							 /*----------------construction of the global stiffness matrix-------------- -*/

							 for (i1 = 0; i1 < 3; i1++){
								 b[mesh[i][i1]] = b[mesh[i][i1]] + b1[i1] + b2[i1]; // source term


							 MatrixEquation_loop2_2_1:for (i2 = 0; i2 < 3; i2++){

								 SS[mesh[i][i1]][mesh[i][i2]] = SS[mesh[i][i1]][mesh[i][i2]] + S1[i1][i2] + S2[i1][i2] + S3[i1][i2]; //inserting local stiffness matrix to global stiffness matrix

							 }
							 }


							 /*----------------construction of the global Jakoonian matrix-------------- -*/
							 //the local jakoobian matrices J2 and J3 are equal to S2 and S3

							 for (i1 = 0; i1 < 3; i1++){
							 MatrixEquation_loop2_3_1:for (i2 = 0; i2 < 3; i2++){

								 SJ[mesh[i][i1]][mesh[i][i2]] = SJ[mesh[i][i1]][mesh[i][i2]] + J1[i1][i2] + S2[i1][i2] + S3[i1][i2]; 	//inserting local stiffness matrix to global Jakoobian matrix
							 }
							 }



						 }


}
int ismember(int x, int n, int* A){

	int i = 0;
ismember_Loop:
	for (; i < n; i++){
		if (A[i] == x)
			return i;
	}
	return -1;

}
void rowDeleting(int row, int col, int target){
	int i, j;
RowDeleting_Loop: for (i = target; i < 700 - 1; i++){
RowDeleting_innerLoop:	for (j = 0; j < 700; j++){
	SS[i][j] = SS[i + 1][j];
	SJ[i][j] = SJ[i + 1][j];
}
}
}
void rowDeleting1(int row, int col, double a[600][20], int target){
	int i, j;
RowDeleting1_Loop: for (i = target; i < 599; i++){
RowDeleting1_innerLoop: for (j = 0; j < 20; j++){
	a[i][j] = a[i + 1][j];
}
}
}
void colDeleting(int row, int col, int target){
	int i, j;
ColDeleting_Loop: for (i = 0; i < 700; i++){
ColDeleting_innerLoop:	for (j = target; j < 699; j++){
	SS[i][j] = SS[i][j + 1];
	SJ[i][j] = SJ[i][j + 1];
}
}
}
void vectorDeleting(int len, double a[600], int target){

	int i;
VectorDeleting_Loop:for (i = target; i < 799; i++){
	a[i] = a[i + 1];


}

}
void DeletingSubGroup(int index){
	rowDeleting(700, 700, index);
	colDeleting(700, 700, index);

}
void BounderyConditions(){

	/*anti periodicity boundary condition on the right hand side of the machine*/


	int l, i, j, m, n;
	double sstmp1, sstmp2, sstmp3, sstmp4, sjtmp1, sjtmp2, sjtmp3, sjtmp4, diff1, diff2;
	int  k1 = 0;
	int fix1[] = { 0, 2, 4, 12, 16, 24, 32, 42, 47, 54, 58, 62, 66, 70, 74, 78, 82, 84 };
	int index[] = { 0, 1, 86, 87, 172, 173, 258, 259, 344, 345, 430, 431, 516, 84, 85, 170, 171, 256, 257, 342, 343, 428, 429, 514, 515, 533 };

BounderyConditions_fix:	for (l = 0; l < 18; l++){

	L_P[l] = fix1[l];

	R_P[l] = l + 516;

}


						//for ( l = 18; l < 18 + n_supp; l++){
						//
						//
						//	}



						//BounderyConditions_1_0:	for ( l = 0; l < 18; l++){




						//	}

					BounderyConditions_Loop1_1:for (l = 18; l < 18 + n_supp; l++){
						L_P[l] = NNO + l - 18;
						R_P[l] = mbu[N_m - n_supp + l - 19];// node number on the anitiperiodicity boundary right
					}
											   // for stifness matrix




										   BounderyConditions_loop2_1:for (i = 0; i < 18 + n_supp; i++){
										   BounderyConditions_loop2_2:for (j = 0; j < 18 + n_supp; j++){

											   SS[R_P[i]][R_P[j]] = SS[R_P[i]][R_P[j]] + SS[L_P[i]][L_P[j]];
											   SJ[R_P[i]][R_P[j]] = SJ[R_P[i]][R_P[j]] + SJ[L_P[i]][L_P[j]];

										   }

										   }
																	  //					BounderyConditions_oloop4:   for (i = 0; i < 18 + n_supp; i++){
																	  //					BounderyConditions_iloop4:for (j = 0; j < 18 + n_supp; j++){
																	  //
																	  //
																	  //
																	  //
																	  //
																	  //					}
																	  //
																	  //					}

																  BounderyConditions_oLoop3:for (i = 0; i < 18 + n_supp; i++){

																  BounderyConditions_iloop3:for (j = 0; j < NNO + n_supp; j++){
																	  if (j == 0){
																		  b[R_P[i]] = b[R_P[i]] - b[L_P[i]];
																	  }
																	  sstmp1 = SS[R_P[i]][j];
																	  sstmp2 = SS[L_P[i]][j];
																	  sstmp3 = SS[j][R_P[i]];
																	  sstmp4 = SS[j][L_P[i]];


																	  if (ismember(j, 18 + n_supp, L_P) == -1){
																		  diff1 = sstmp1 - sstmp2;
																		  diff2 = sstmp3 - sstmp4;
																		  SS[R_P[i]][j] = diff1;
																		  SS[j][R_P[i]] = diff2;
																	  }
																  }

																  }



																						BounderyConditions_oLoop5:for (i = 0; i < 18 + n_supp; i++){

																						BounderyConditions_iLoop5:for (j = 0; j < NNO + n_supp; j++){
																							sjtmp1 = SJ[R_P[i]][j];
																							sjtmp2 = SJ[L_P[i]][j];
																							sjtmp3 = SJ[j][R_P[i]];
																							sjtmp4 = SJ[j][L_P[i]];
																							if (ismember(j, 18 + n_supp, L_P) == -1){
																								SJ[R_P[i]][j] = sjtmp1 - sjtmp2;
																								SJ[j][R_P[i]] = sjtmp3 - sjtmp4;
																							}
																						}

																						}


																												  //BEL = [1, 2, 6, 5, 12, 11, 20, 28, 27, remesh(96 - 2 * n_supp + 1):remesh(96), remesh(1), remesh(2), 40, 39, 52, 53, 54, 60, 61, 62, 63, 68, 69, 70, 74, 75, 76, 81, 82, 83, 84, 89, 90, 91, 92, 97, 96] % element numbers on the boundary that have at least one node on the boundary

																												  // Drichlet boundary condition(A = 0) on top and bottom sides of the geometry
																												  // in order to make BC nodes to have zero potential we have to make zero the corresponding nodes potential of 1 - An 2 - SJ\(SS*An(:, I) - b) which results in :


																												  /* { 1, 2, 87, 88, 173, 174, 259, 260, 345, 346, 431, 432, 517, 85, 86, 171, 172, 257, 258, 343, 344, 429, 430, 515, 516, 534 }; // node numbers at which A = 0*/
																												  for (m = 0; m < 26; m++){

																													  for (n = 0; n < 700; n++){
																														  SS[index[m]][n] = 0;
																														  SJ[index[m]][n] = 0;
																														  if (n == 699){
																															  SS[index[m]][index[m]] = 1;
																															  SJ[index[m]][index[m]] = 1;
																															  b[index[m]] = 0;


																														  }
																													  }

																												  }


																												  for (i = 0; i < 18; i++){


																													  //DeletingGroup(fix1[i]-k1);
																													  DeletingSubGroup(fix1[i] - k1);


																													  vectorDeleting(800, b, fix1[i] - k1);

																													  rowDeleting1(600, 20, An, fix1[i] - k1);

																													  k1++;
																												  }


																											  BounderyConditions_oLoop8:for (i = NNO; i < NNO + n_supp; i++){

																												  //DeletingGroup( i - k1);
																												  DeletingSubGroup(i - k1);

																												  vectorDeleting(800, b, i - k1);

																												  rowDeleting1(600, 20, An, i - k1);
																												  k1++;
																											  }


																																		//for ( i = NNO; i < NNO + n_supp; i++){
																																		//	rowDeleting(700, 700, i - k1);
																																		//	colDeleting(700, 700,i - k1);
																																		//	vectorDeleting(800, b, i - k1);
																																		//	//rowDeleting(700, 700, SJ, i - k1);
																																		//	//colDeleting(700, 700, SJ, i - k1);
																																		//	rowDeleting1(600, 20, An, i - k1);

																																		//	k1++;
																																		//}


}
void minus(double a[], double b[], double res[], int len){
minus_lopp:for (int i = 0; i < len; i++){
	res[i] = a[i] - b[i];

}


}
void norm(double* a, int len, double* res){
	double temp = 0.f;
norm_loop:for (int i = 0; i < len; i++){
	temp += a[i] * a[i];

}

		  *res = sqrt(temp);



}
void translate1Down(int destend, int deststart, int col, double a[600][20]){
	int l = 0;
	int n = 0;
translate1Down_loop:for (l = destend; l >deststart - 1; l--){
translate1Down_iloop:for (n = 0; n < col; n++){

	a[l][n] = a[l - 1][n];

}

}

}

void GPLU(double A1[700][700], double BB[516], double x1[600]) {
	/* solve A1*x1 = BB */

	int n = 516; /* size of A1 */

	int s, i, j, k, m, in, ii, len, l, levelindex, first;

	int levelnum = 0; /* how many levels */
	int levellength[516] = { 0 }; /*length of every level*/
	int level[516][516] = { 0 };/*columns of every level*/
	int ancester[516] = { 0 };
	int nonZeroIndexLU[516][460] = { 0 };/*non zero indices of matrix A1*/
	double As[516][516] = { 0 };
	double B[516] = {0};
	int levelOf[516] = { 0 }; /*column dependency level*/
	double y1[516] = { 0 };

	s = 0;
	for (i = 0; i < n; i++)
		B[i] = BB[i];
	for (i = 0; i < n; i++) {
	GPLU_loop0: for (j = 0; j < 460; j++) {
		nonZeroIndexLU[i][j] = -1;
	}
	}
	for (i = 0; i < n; i++) {
	GPLU_loop1: for (j = 0; j < n; j++) {
		As[i][j] = A1[i][j];
	}
	}

	/* LU decomposition using GP */

	/*find new nonzeros (fill in) and column dependencies*/
	for (i = 0; i < n; i++) {
		s = 0;
		first = 1;
	GPLU_First: for (j = 0; j < n; j++) {
		if (As[j][i] != 0) {
			As[j][i] = 1;

			nonZeroIndexLU[i][s] = j;
			s++;
			if (i > 0) {
				if (j < i) {
					/*compute dependency*/
					if (levelOf[j] + 1 > levelOf[i]) {
						ancester[i] = j;
						levelOf[i] = levelOf[j] + 1;

						if (first == 1) {
							levelnum++;
							first = 0;
						}
					}

					/*find fill in*/
				GPLU_loop2: for (int m = 0; m < 460; m++) {
					in = nonZeroIndexLU[j][m];
					if (in != -1) {
						As[in][i] = 1;
					}
					else {
						break;
					}

				}

				}
			}
		}

	}
				level[levelOf[i]][levellength[levelOf[i]]++] = i;
	}

	/*L and U calculation, in column dependency order. Columns in same level should compute in parellel*/

GPLU_LU: for (l = 0; l < levelnum; l++) {
	len = levellength[l];
	for (levelindex = 0; levelindex < len; levelindex++) {
		k = level[l][levelindex];
	GPLU_loop3: for (i = 0; i < n; i++) {
		x1[i] = A1[i][k];
	}

			GPLU_outX1: for (m = 0; m < 460; m++) {
				j = nonZeroIndexLU[k][m];
				if (j > -1) {
					if (j >= k) {
						continue;
					}
					else {
					GPLU_loop4: for (ii = 0; ii < 460; ii++) {
						in = nonZeroIndexLU[j][ii];
						if (in != -1) {
							if (in <= j) {
								continue;
							}
							else {
								x1[in] -= As[in][j] * x1[j];
							}
						}
						else {
							break;
						}

					}
					}
				}
				else {
					break;
				}
			}
					GPLU_loop5: /*for (i = 0; i < 460; i++) {
								j = nonZeroIndexLU[k][i];
								if (j > -1)
								if (j <= k) {
								As[j][k] = x1[j];
								}
								else {
								As[j][k] = x1[j] / As[k][k];
								}
								else
								break;

								}*/
						for (i = 0; i < 516; i++){
							if (i <= k)
								As[i][k] = x1[i];
							else
								As[i][k] = x1[i] / As[k][k];
						}
	}
}



/* solve Ax = b using L and U*/

/*-----backward substitution----*/
/* ----Ly = BB---- */
GPLU_back0out: for (k = 0; k < n; k++) {
y1[k] = B[k];

GPLU_backloopi1: for (i = 0; i < 460; i++) {
j = nonZeroIndexLU[k][i];
if (j == -1) {
	break;
}
if (j < k) {
	continue;
}
else {
	if (j == k) {
		B[k] -= y1[k];
	}
	else {
		B[j] = B[j] - y1[k] * As[j][k];
	}
}
}
}

			   /*----Ux = y----*/
		   GPLU_init2: for (i = 0; i < 516; i++) {
			   x1[i] = 0;
		   }
				   GPLU_back1out: for (k = n - 1; k >= 0; k--) {
					   x1[k] = y1[k] / As[k][k];

				   GPLU_backloopi2: for (i = 0; i < 460; i++) {
					   j = nonZeroIndexLU[k][i];

					   if (j == -1) {
						   break;
					   }
					   else {
						   if (j <= k && j >= 0) {
							   y1[j] -= x1[k] * As[j][k];
						   }
					   }
				   }

				   }

}



int NaiveSolver(int I, double* nares){
	int i = 0;
	int j = 0;
	int k = 0;
	int m = 0;
	int n = 0;
	int is = 0;
	int js = 0;
	int l = 0;
	double norm1, norm2;
	int k3 = 0;
	int fix1[] = { 0, 2, 4, 12, 16, 24, 32, 42, 47, 54, 58, 62, 66, 70, 74, 78, 82, 84 };
	// expanding the An to the NNO + N_supp matrix required for the next iteration
	int k2 = 1;



	ns = 516;
	res2 = 0.f;
	for (i = 0; i < 530; i++){
		for (j = 0; j < 530; j++){
			num[i][j] = 0;
		}
	}
	for (i = 0; i < 600; i++){
		res1[i] = 0;
		res3[i] = 0;
		res4[i] = 0;
		frac[i] = 0;
		if (i < 516)
			den[i] = 0;
	}




	for (m = 0; m < 516; m++){
	NaiveSolver_loop1:for (n = 0; n < 516; n++){
		num[m][n] = SJ[m][n];

	}
	}

	//double res1[516] = {0};
NaiveSolver_loop3:for (i = 0; i < 516; i++){
NaiveSolver_loop3_1:for (j = 0; j < 516; j++){


	res1[i] += SS[i][j] * An[j][I];

}
}



				  //double den[516] = {0};
				  minus(res1, b, den, 516);


				 


			  //NaiveSolver_loop5_1:for (m = 0; m < 516; m++){

				 // num[m][516] = den[m];

			  //}
				/*  if (I == 1){
					  printf("Ssssssssssssssssssssssssssssss\n");
					  for (i = 0; i < 516; i++)
						  printf("%f\t",An[i][I]);
				  }*/
			
								  /*-----------------test GP------------------------*/
				 
				//  printf("\n");
				  //if (I == 1){
					
					 // for (i = 0; i < 20; i++){
						//  //for (j = 0; j < 20; j++){
						//  printf("%f\t", b[i]);
						//  //}
					 // }
					 // for (i = 0; i < 20; i++){
						//  //for (j = 0; j < 20; j++){
						//  printf("%f\t",res1[i]);
						//  //}
					 // }
					
				  //}

				  GPLU(SJ, den, frac, levelOf, ancester);
								  /*------------------------------------------------*/
				 /* if (I == 1){
					  for (i = 0; i < 20; i++){
						  printf("%f\t", frac[i]);
					  }
					 
						  return;
				  }
				 */
				  
																				NaiveSolver_an:for (m = 0; m < 600; m++){
																					An[m][I + 1] = An[m][I] - frac[m];
																					//if (I ==0)
																					//printf("%f - %f\t", An[m][I], frac[m]);

																				}



																							 //  printf("%f, %f", SJ[70][0], SJ[499][0]);

																							   /*Convergency criteria*/


																						   NaiveSolver_cinv1:for (m = 0; m < 600; m++){
																							   res3[m] = An[m][I + 1] - An[m][I];
																							   res4[m] = An[m][I];
																						   }
																											 //if (I == 1){
																											 //	for (int n = 0; n < 600; n++){
																											 //		printf("%d %f %f\n", n, An[n][I ] , An[n][I+1]);
																											 //	}
																											 //}
																											 norm1 = 0.f;
																											 norm(res3, 600, &norm1);

																											 norm2 = 0.f;
																											 norm(res4, 600, &norm2);
																											 *nares = norm1 / norm2;

																											 //	printf("%f,%f,%f\n", norm1, norm2, norm1 / norm2);


																											 if (norm1 / norm2 < 5e-3)   {//convergency condition
#ifndef __SYNTHESIS__
																												 printf(" convergent!\n");
#endif
																												 return -1; //when the solution is found, the for loop will break with a 516 * 1 An
																											 }
																											 if (I == 4){
#ifndef __SYNTHESIS__
																												 printf(" non convergent!\n");
#endif

																												 return 0;

																											 }




																											 for (m = 0; m < 18; m++){
																												 /*	for (int l = 515+k2; l >fix1[m] ; l--){
																												 for (int n = 0; n < 3; n++){

																												 An[l ][n] = An[l-1][n];

																												 }

																												 }*/
																												 translate1Down(515 + k2, fix1[m] + 1, 20, An);


																											 NaiveSolver_conv2:	for (l = 0; l < 3; l++){
																												 An[fix1[m]][l] = 0;
																											 }

																																k2++;
																											 }


																										 NaiveSolver_transAnLoop2:for (m = NNO; m < NNO + n_supp; m++){
																										 NaiveSolver_transAninLoop:	for (l = 515 + k2; l >m; l--){
																										 NaiveSolver_Trans2:	for (n = 0; n < 3; n++){
																											 An[l][n] = An[l - 1][n];
																										 }
																										 }

																																NaiveSolver_An0:for (l = 0; l < 3; l++){
																																	An[m][l] = 0;
																																}

																																				k2++;
																										 }

																																  // these two loops should be seperate as it is, should not be combined


																															  NaiveSolver_Trans3:	for (m = 0; m < 18; m++){
																															  NaiveSolver_trans3in:	for (n = 0; n < 3; n++){
																																  An[fix1[m]][n] = -An[R_P[k3]][n];
																															  }
																																					k3++;
																															  }

																																				NaiveSolver_transAnLop4:for (m = NNO; m < NNO + n_supp; m++){
																																				NaiveSolver_trans4:	for (n = 0; n < 3; n++){

																																					An[m][n] = -An[R_P[k3]][n];
																																				}
																																									k3++;

																																				}




																																										return 0;
}

int NonLinearIteration(double t, double* n){
	int i = 0;
NonLinearIteration_loop1:for (i = 0; i < NNO + n_supp; i++){
	An[i][0] = 0.001;       // the first iteration, initial solution of the corresponding time step
}


					 NonLinearIteration_loop2:for (I = 0; I < 20; I++){

						 SpecifyingMaterials(I);

						 MatrixEquation(t, I);

						 BounderyConditions();

						 res = NaiveSolver(I, n);
						// printf("%d\n", I);
						//
						// if (I == 1)break;

						 if (res == -1){
#ifndef __SYNTHESIS__
							  printf("I = %d\n", I);
#endif
							 return -1;

						 }


					 }
											  return 0;
}
void ExpandAt(int I, double t){


	int i = 0;
	int ddivres1 = 0;
	int fix1[] = { 0, 2, 4, 12, 16, 24, 32, 42, 47, 54, 58, 62, 66, 70, 74, 78, 82, 84 };
	// expanding the An to the NNO + N_supp matrix required for the next iteration
	int k2 = 1;
	int m, l;
	int k3 = 1;


ExpandAt_loop1:for (i = 0; i < 516; i++){

	At[i] = An[i][I + 1];


}


			   // expanding the At to the NNO + N_supp matrix required for the next time step

		   ExpandAt_out2:for (m = 0; m < 18; m++){
		   ExpandAt_inner2:for (l = 515 + (k2++); l > fix1[m]; l--){

			   At[l] = At[l - 1];

		   }
		   }
					 ExpandAt_out2_5:for (m = 0; m < 18; m++){


						 At[fix1[m]] = 0;



					 }


								 ExpandAt_outloop3:for (m = NNO; m < NNO + n_supp; m++){
								 ExpandAt_inLoop3:for (l = 515 + (k2++); l >m; l--){
									 At[l] = At[l - 1];
								 }



								 }

											   ExpandAt_outloop3_1:for (m = NNO; m < NNO + n_supp; m++){


												   At[m] = 0;



											   }
																   // these two loops should be seperate as it is, should not be combined



															   ExpandAt_loop3_5:for (m = 0; m < 18; m++){

																   At[fix1[m]] = -At[R_P[k3]];

																   k3++;
															   }
																			ExpandAt_Loop4:for (m = NNO; m < NNO + n_supp; m++){

																				At[m] = -At[R_P[k3]];

																				k3++;

																			}

																						   ddivres1 = (int)(t * dtrev + 2) - 1;
																					   ExpandAt_loop5:for (i = 0; i < NNO + n_supp; i++){
																						   A[i][ddivres1] = At[i];//saving expanded At to A as the final solution of each time step

																					   }


}
void ForceCalculation(double t){


	int i, j;
	int ddivres1;
	double detrev;
	int ii[] = { 26, 27, 27, 28, 29, 30, 30, 31, 32, 33, 33, 34, 35, 36, 36, 37 };/*{ 27, 28, 28, 29, 30, 31, 31, 32, 33, 34, 34, 35, 36, 37, 37, 38 };*/
	int add[] = { 0, 168, 336, 504, 672, 840 };


ForceCalculation_loop1:for (i = 0; i < 16; i++){// the element numbers corresponding to the force calculation band
ForceCalculation_loop2:for (j = 0; j < 6; j++){

	N1 = mesh[ii[i] + add[j]][0];                                                      // Global node number of that element
	N2 = mesh[ii[i] + add[j]][1];
	N3 = mesh[ii[i] + add[j]][2];
	DET = X[N2] * Y[N3] + X[N1] * Y[N2] + X[N3] * Y[N1] - X[N1] * Y[N3] - X[N3] * Y[N2] - X[N2] * Y[N1];
	detrev = 1 / DET;
	ddivres1 = (int)(t * dtrev + 2) - 1;
	//B = curl[A] = dA / dy i - dA / dx j][ where A = a1 + a2x + a3y][ so dA / dy = a3 and dA / dx = a2
	a2 = (A[N2][ddivres1] * Y[N3] - A[N3][ddivres1] * Y[N2] - A[N1][ddivres1] * (Y[N3] - Y[N2]) + Y[N1] * (A[N3][ddivres1] - A[N2][ddivres1])) *detrev;
	a3 = (A[N3][ddivres1] * X[N2] - A[N2][ddivres1] * X[N3] - X[N1] * (A[N3][ddivres1] - A[N2][ddivres1]) + A[N1][ddivres1] * (X[N3] - X[N2])) *detrev;
	F_thrust[ddivres1] = F_thrust[ddivres1] + 21 * L*a2*a3 / (4 * pi*1e-7) * 0.0254 / 16;
	F_Normal[ddivres1] = F_Normal[ddivres1] + 21 * L*(a2 * a2 - a3 * a3) / (2 * 4 * pi*1e-7) * 0.0254 / 16;
}
}

}

void LoadData(int m[1008][3], double x[600], double y[600]){
	int i = 0;
	int j = 0;


	for (i = 0; i<1008; i++){

		for (j = 0; j<3; j++){
			mesh[i][j] = m[i][j];

		}
	}

	for (i = 0; i<534; i++){
		X[i] = x[i];
		Y[i] = y[i];

	}


}
void readFile(int mesh[1008][3], double X[600], double Y[600]){
	FILE *fp;
	int i, j;

	fp = fopen("mesh.txt", "r");
	for (i = 0; i < 1008; i++)
		for (j = 0; j < 3; j++)
			fscanf(fp, "%d", &mesh[i][j]);

	fp = fopen("X.txt", "r");
	for (i = 0; i < 534; i++)
		fscanf(fp, "%lf", &X[i]);

	fp = fopen("Y.txt", "r");
	for (i = 0; i < 534; i++)
		fscanf(fp, "%lf", &Y[i]);

}
void TOP(double naivearray[1000], int m[1008][3], double x[600], double y[600]){
	double t;
	int i = 0;
	int j = 0;
	time_t timeBegin, timeEnd;
	
	LoadData(m, x, y);
	Initialize();
	//readFile(mesh, X, Y);


	MeshSpeedAndConductivity();



TopLoop:for (t = 0; t < T; t += dt)
{


	MovingBandModel();

	SpecifyingSource(t);


	timeBegin = clock();
	//dwBegin = timeGetTime();

	int conv = NonLinearIteration(t, &naivearray[(int)(t * dtrev)]);
	

	//break;

	//dwEnd = timeGetTime();
	timeEnd = clock();
	printf("%f ms\n", (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC);




	ExpandAt(I, t);


	ForceCalculation(t);
#ifndef __SYNTHESIS__
	//printf("t = %f\n", t);
#endif
}



}

