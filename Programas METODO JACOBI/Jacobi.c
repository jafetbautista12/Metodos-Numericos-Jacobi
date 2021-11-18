#include <stdio.h>
#include <stdlib.h>

int dim;
float norma(float vector1[],float vector2[]);
float suma_jacobi(float Matriz[], float vector[], int componente);

int main(){
    int i,j,iteraciones=0;
    float error,epsilon;
    printf("\n METODO DE JACOBI DE RESOLUCION DE SISTEMAS Ax=b \n");

    printf("Dimension de la matriz A: ");
    scanf("%d",&dim);
    float A[dim][dim],b[dim],x[dim],x_prev[dim],aux[dim];

    printf("\n Elementos de la matriz A: \n");
    for(i=0;i<dim;i++) for(j=0;j<dim;j++){
        printf("A(%d,%d)=",i,j); scanf("%f",&A[i][j]);
    }

    printf("\n Elementos del vector b: \n");
    for(i=0;i<dim;i++){
        printf("b(%d)=",i); scanf("%f",&b[i]);
    }

    printf("\n Error de parada: \n");
    printf("E=",i); scanf("%f",&epsilon);
    error=epsilon+1;

    //cominezo algoritmo de Jacobi
    //Error se mide como la norma del vector diferenceia entre la iteracion i e i+1
    printf("\n Valor inicail de la iteracion: \n");
    for(i=0;i<dim;i++){
        printf("x0(%d)=",i); scanf("%f",&x_prev[i]);
    }
    while (error>epsilon){
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++) aux[j]=A[i][j];
            x[i]=(1/A[i][i])*(b[i]-suma_jacobi(aux,x_prev,i));
        }
        error=norma(x,x_prev);

        printf("\n\n Iteracion %d: \n",iteraciones);
        for(i=0;i<dim;i++){
            x_prev[i]=x[i];
            printf("X(%d)=%f \n",i,x[i]);
        }

        iteraciones++;
        if (iteraciones==10) error=epsilon-1;
    }

    printf("Solucion del sistema\n");
    printf("Numero de iteraciones: %d \n", iteraciones);
    for(i=0;i<dim;i++){
        printf("x(%d)=%f\n",i,x[i]);
    }
    return 1;
}

float norma(float vector1[],float vector2[]){
    float aux=0;
    int i;
    for(i=0;i<dim;i++){
        aux=aux+(vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
    }
    return aux;
}

float suma_jacobi(float Matriz[], float vector[], int componente)
{
    float aux=0;
    int i;
    for(i=0;i<dim;i++){
        if (componente!=i){
            aux=aux+Matriz[i]*vector[i];
        }
    }
    return aux;
}
