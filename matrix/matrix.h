#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED
double** beolvas(char* hely);
double* soronkent(char* sor);
void cout(double** matrix);
double** sum(double** M1, double** M2);
double** sub(double** M1, double** M2);
double** mult(double** M1, double** M2);
double** fmem(double** M, int oszlop);
double** almatrix(double** M, int oszlop, int sor);
double det(double** A);
void freem(double** M);
void mulc(double** M, double lambda);
double** kulonsorra(double** M, int sor);
double** transponalt(double** M);
void sorcsere(double** M, int i, int j);
double** adj(double** M);
double** inverse(double** M);
void letisztaz(double** M);
double** deepcpy(double** M);
double** GJE(double** M);
void fkiir(char* hely, double** mit);
void Mmalloc(double** M,int i, int j);
#endif // MATRIX_H_INCLUDED
