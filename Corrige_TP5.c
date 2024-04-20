#include <stdio.h>
#include <stdlib.h>
#define N 4 //permet de changer plus rapidement la dimension si besoin (remplace N par 4 après compilation)

double A1[N] = {1,3,4,1};
double A2[N] = {2,4,3,4};
double A3[N] = {1,1,2,5};
double A4[N] = {1,4,5,2};
double* A[N] = {A1, A2, A3, A4};
double b[N] = {1,2,3,3};



void affiche_mat(double** A){
    for (int i = 0; i<N;i++) {
        for (int j = 0; j<N;j++) printf("%.2f ", A[i][j]);
        printf("\n");
    }
    printf("\n");
}

void affiche_vect(double* x){
    for (int j = 0; j<N;j++) printf("%.2f ",x[j]);
        printf("\n");
        printf("\n");
}

//question 1
void mult_vec(double** A, double* x, double* y){
	for (int i = 0; i<N;i++){
		y[i] = 0;
		for (int j= 0;j<N;j++){
			y[i] += A[i][j]*x[j];
		}
	}
}

//question 2
void mult_mat(double** A, double** B, double** C){
	for (int i = 0; i<N;i++){
		for (int j = 0;j<N;j++){
			C[i][j] = 0;
			for (int k=0;k<N;k++) C[i][j] += A[i][k]*B[k][j];
		}
	}
}

//question 3
void echange(double** A, double* b, int k, int l){
	double* p = A[k];
	double x = b[k];
	A[k] = A[l];
	b[k] = b[l];
	A[l] = p;
	b[l] = x;
}

//question 4
void ajout_scal(double** A, double* b, double s, int k, int l){
	for (int i=0;i<N;i++){
		A[k][i] += s*A[l][i];
	}
	b[k] += s*b[l];
}


//question 5
double fabs(double x){
	if (x<0) return -x;
	return x;
}

int max_ligne(double** A, int i){
	int max=i;
	for (int j=i+1;j<N;j++){
		if (fabs(A[j][i]) > fabs(A[max][i])) max = j;
	}
	return max;
}

//question 6
void pivot(double** A, double* b){
	for (int i = 0; i<N;i++){
		int k = max_ligne(A, i);
		echange(A,b,i,k);//après l'échange, A[i][i] est bien le max de la colonne i
		for (int j=i+1;j<N;j++){
			ajout_scal(A, b, -A[j][i]/A[i][i], j ,i);
		}
	}
}

//question 7
void gauss(double** A, double* b, double* y){
	pivot(A,b); //on transforme notre systeme pour que A soit triangulaire
	//a partir d'une matrice triangulaire, on peut facilement calculer la solution
	for (int i=N-1;i>=0;i--){
		y[i] = b[i];
		for (int j=i+1;j<N;j++){
			y[i] += -A[i][j]*y[j];
		}
		y[i] = y[i]/A[i][i];
	}
}

//question 8
void init_col(double* x, int i){//cette fonction calcule la i-eme colonne de l'identité dans X
	for (int j=0;j<N;j++) x[j] = 0;
	x[i] = 1;
}

void copy_mat(double** A, double** Cop){//cette fonction copie A dans Cop
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++) Cop[i][j]=A[i][j];
    }
}

void inv_gauss(double** A, double** A1){
	double col[N];
	double y[N];
	double* Cop[N]; for(int i=0;i<N;i++) Cop[i] = malloc(N*sizeof(double));
	for (int i=0;i<N;i++){
		//pour chaque colonne c de I, on calcule la solution de Ax = c
		copy_mat(A,Cop); //comme on veut reutiliser A en l'état apres gauss, on doit la copier dans une nouvelle matrice
		init_col(col, i);
		gauss(Cop,col,y); //on calcule dans y la solution de Ax = c
		for (int j =0;j<N;j++){ //en fait, la solution est egale à la ieme colonne de A^-1.
			A1[j][i] = y[j];
		}
	}
}


int main(void){
	double* A1[N];
	for (int i =0;i<N;i++) A1[i] = malloc(N*sizeof(double));
	double* A2[N];
	for (int i =0;i<N;i++) A2[i] = malloc(N*sizeof(double));
	inv_gauss(A, A1);
	affiche_mat(A1);
	mult_mat(A,A1,A2); //pour vérifier notre inverse, il suffit de vérifier que A*A1 = I.
	affiche_mat(A2);
	return 0;
}
