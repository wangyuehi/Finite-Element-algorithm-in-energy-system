#ifndef _FEMTESTGPV1_H
#define _FEMTESTGPV1_H

#include <stdio.h>

void readFile(int mesh[1008][3], double X[600], double Y[600]);
void Initialize();
void MeshSpeedAndConductivity();
void MovingBandModel();
void SpecifyingSource(double t);
void SpecifyingMaterials(int I);
void MatrixEquation(double t, int I);
int ismember(int x, int n, int* A);
void rowDeleting(int row, int col, int target);
void rowDeleting1(int row, int col, double a[600][20], int target);
void colDeleting(int row, int col, int target);
void vectorDeleting(int len, double a[600], int target);
void BounderyConditions();
void multOneCol(double** A, double** b, int n, double* c, int len);
void minus(double a[], double b[], double res[], int len);
void norm(double* a, int len, double* res);
void translate1Down(int destend, int deststart, int col, double a[600][20]);
int NaiveSolver(int I, double* nares);
int NonLinearIteration(double t, double* n);
void ExpandAt(int I, double t);
void ForceCalculation(double t, double F_thrust1test);
void TOP(float naivearray[1000], int m[1008][3], double x[600], double y[600]);

#endif