
#include<stdio.h>
#include<math.h>
#define ZERO 1.0E-20
#define true 1
#define false 0

double absval(double);
void INPUT(int *, double [][10], int *);
void OUTPUT(int, double [][10]);

main()
{
   double A[10][10], V[10], U[10], Z[10];
   double S,Q,RSQ,PROD;
   int N,I,J,K,KK,L,OK;

   INPUT(&OK, A, &N);
   if (OK) {
      /* STEP 1 */
      for (K=1; K<=N-2; K++) {
         Q = 0.0;
         KK = K + 1;
         /* STEP 2 */
         for (I=KK; I<=N; I++) Q = Q + A[I-1][K-1] * A[I-1][K-1];
         /* STEP 3 */
         /* S is used in place of alpha.  */
         if (absval(A[K][K-1]) <= ZERO) 
            S = sqrt(Q);
         else 
            S = A[K][K-1] / absval(A[K][K-1]) * sqrt(Q);
         /* STEP 4 */
         RSQ = (S + A[K][K-1]) * S;
         /* STEP 5 */
         V[K-1] = 0.0;
         V[K] = A[K][K-1]+S;
         for (J=K+2; J<=N; J++) V[J-1] = A[J-1][K-1];
         /* STEP 6 */
         for (J=K; J<=N; J++) {
            U[J-1] = 0.0;
            for (I=KK; I<=N; I++) U[J-1] = U[J-1] + A[J-1][I-1]*V[I-1];
            U[J-1] = U[J-1] / RSQ;
         }  
         /* STEP 7 */
         PROD = 0.0;
         for (I=K+1; I<=N; I++) PROD = PROD + V[I-1]*U[I-1];
         /* STEP 8 */
         for (J=K; J<=N; J++) Z[J-1] = U[J-1] - 0.5*PROD*V[J-1]/RSQ;
         /* STEP 9 */
         for (L=K+1; L<=N-1; L++) {
            /* STEP 10 */
            for (J=L+1; J<=N; J++) {
               A[J-1][L-1] = A[J-1][L-1]-V[L-1]*Z[J-1]-V[J-1]*Z[L-1];
               A[L-1][J-1] = A[J-1][L-1];
            }  
            /* STEP 11 */
            A[L-1][L-1] = A[L-1][L-1] - 2.0*V[L-1]*Z[L-1];
         }  
         /* STEP 12 */
         A[N-1][N-1] = A[N-1][N-1]-2.0*V[N-1]*Z[N-1];
         /* STEP 13 */
         for (J=K+2; J<=N; J++) {
            A[K-1][J-1] = 0.0;
            A[J-1][K-1] = 0.0;
         }  
         /* STEP 14 */
         A[K][K-1] = A[K][K-1]-V[K]*Z[K-1];
         A[K-1][K] = A[K][K-1];
      }  
      /* STEP 15 */
      OUTPUT(N, A);
   }
   return 0;
}

void INPUT(int *OK, double A[][10], int *N)
{
   int I, J, FLAG;
   char AA;
   char NAME[30];
   FILE *INP; 

   printf("Metodo de householder .\n");
   *OK = false;
   printf("la matriz simetrica sera un arreglo que se mostrara \n");
   printf("en el orden:\n");
   printf("              A(1,1), A(1,2), A(1,3), ..., A(1,n),\n");
   printf("                      A(2,2), A(2,3), ..., A(2,n),\n");
   printf("                              A(3,3), ..., A(3,n),\n");
   printf("                                      ..., A(n,n)\n\n");
   printf("coloque las entradas que desee en cada una de las lineas");
   printf("entradas con\n");
   printf("recuerde que debe de haber al menos un espacio en blanco.\n\n\n");
   printf("Usted ha creado la entrada? - continuar Y or N.\n");
   scanf("%c",&AA);
   if ((AA == 'Y') || (AA == 'y')) {
      printf("Ingrese el nombre del archivo en el formulario - drive:name.ext\n");
      printf("por ejemplo:   A:DATA.DTA\n");
      scanf("%s", NAME);
      INP = fopen(NAME, "r");
      *OK = false;
      while (!(*OK)) {
         printf("ingrese la dimension n.\n");
         scanf("%d", N);
         if (*N > 1) {
            for (I=1; I<=*N; I++) 
               for (J=I; J<=*N; J++) {
                  fscanf(INP, "%lf", &A[I-1][J-1]);
                  A[J-1][I-1] = A[I-1][J-1];
               }  
            fclose(INP);   
            *OK = true;
         }  
         else printf("recuerde que la dimension debe de ser mayor que  1.\n");
      }  
   }  
   else printf("El programa finalizará para que se pueda crear el archivo de entrada..\n");
}

void OUTPUT(int N, double A[][10])
{
   int FLAG,I,J;
   char NAME[30];
   FILE *OUP;

   printf("Seleccione el metodo de salida de resultado como mejor le convenga:\n");
   printf("1. Salida a pantalla\n");
   printf("2. \n");
   printf("Por favor ingresa  1 or 2.\n");
   scanf("%d", &FLAG);
   if (FLAG == 2) {
      printf("ingresa el nombre del archivo - drive:name.ext\n");
      printf("por ejemplo   A:OUTPUT.DTA\n");
      scanf("%s", NAME);
      OUP = fopen(NAME, "w");
   }
   else OUP = stdout;
   fprintf(OUP, "metodo de householder\n\n");
   fprintf(OUP, "la matriz tridiagonal sera - output by rows\n\n");
   for (I=1; I<=N; I++) { 
       for (J=1; J<=N; J++) fprintf(OUP, " %11.8f", A[I-1][J-1]);
       fprintf(OUP, "\n\n");
   }
   fclose(OUP);
}   

/* Absolute Value Function */
double absval(double val)
{
   if (val >= 0) return val;
   else return -val;
}
