/*  ncuerpos calcula la trayectoria de n cuerpos debido a la atracción gravitacional entre ellos. En base a su masa y sus posiciones y velocidades inciales.
  Compilar: 
    $ g++ ncuerpos2.cpp -o ncuerpos
  Ejecutar:
    $ ./ncuerpos
  Estudiantes:
    Juan Andrés Guarín Rojas, 2201870.
    Juan Diego Figueroa Hernandez, 2200815.
    Gabriela Sánchez Ariza, 2200816.
  Escuela de física, Facultad de Ciencias. UIS. 2022.

*/

// Importación de librerias y definiciones.
#include<iostream>
#include<math.h>
#include <stdio.h>
#include "omp.h"

#define G (6.67E-11)  // Constante de la ley de gravitación universal.
using namespace std;

// Definición de funciones.
double* aceleracion(int i, int j, double x[], double y[], double m[]){

    // Función que halla la aceleración gravitacional del cuerpo j-ésimo debido a la fuerza ejercida por el cuerpo i-ésimo.
    // i,j son números enteros que indican el subíndice de cada cuerpo.
    // x,y son vectores de tamaño n, que contienen las posiciones de los n-cuerpos.
    // m es un vector de tamaño n, con los datos de las masas.
  
    double rij, ax, ay;
    double* a=new double[2];
    
    rij = pow(x[i]-x[j],2) + pow(y[i]-y[j],2);   //Distancia al cuadrado entre i y j.
    ax = -G*m[i]*(x[j]-x[i])/pow(rij,1.5);        //Se calculan las aceleraciones en "x" y "y".
    ay = -G*m[i]*(y[j]-y[i])/pow(rij,1.5);
    a[0] = ax;
    a[1] = ay;
    return a; // Aceleraciones en "x" y "y".
}

double* Fk(int k, double x[], double y[], double vx[], double vy[], int n, double m[]){

      // Esta función evalúa la superposición de fuerzas para hallar la aceleración neta. Hace las veces de f(t, y(t)).
    // k indica el subíndice del cuerpo al que se le halla la superposición de fuerzas.
    // x,y,vx,vy son vectores de tamaño n con las posiciones y velocidades de los n cuerpos.
    // n es la cantidad de cuerpos.
    // m es un vector de tamaño n con los datos de las masas.
  
  double akx=0, aky=0; // Inicializa aceleraciones en "x" y "y" del cuerpo k-ésimo en cero.
  double* dat=new double[4]; 
  double* aki=new double[2];
  for(int i=0; i<n; i++){ // Itera sobre los n-cuerpos.
    if (i != k){ // Evita calcular sobre el propio cuerpo de estudio.
    aki = aceleracion(i,k,x,y,m);  // Halla la aceleración sobre i debido al cuerpo k.
    akx += aki[0]; 
    aky += aki[1];
    }
  }
  dat[0]=vx[k];
  dat[1]=vy[k];
  dat[2]=akx;
  dat[3]=aky;
  return dat; // Velocidades y aceleraciones del cuerpo k-ésimo
}

double** siguiente_valor(double x[], double y[], double vx[], double vy[], int n, double h, double m[]){

    // Esta función halla los siguientes valores de posición y velocidad ejecutando el algoritmo de Runge-Kutta Cuarto Orden.
    // x,y,vx,vy son vectores de tamaño n con las posiciones y velocidades de los n Cuerpos.
    // n es la cantidad de cuerpos.
    // h es el paso de tiempo.
    // m es un vector de tamaño n con los datos de las masas.

  double k1[n][4];  // Se calculan los valores de k1
  double* dat=new double[4];
  
  for(int i=0; i<n; i++){
    dat = Fk(i,x,y,vx,vy,n,m);
    k1[i][0] = h*dat[0];
    k1[i][1] = h*dat[1];
    k1[i][2] = h*dat[2];
    k1[i][3] = h*dat[3];
  }
  
  double x1[n], y1[n], vx1[n], vy1[n]; //Listas que contendrán los valores asociados al término k1 del método RK4
  
  for(int i=0; i<n; i++){
    x1[i]=(x[i]+k1[i][0]/2);
    y1[i]=(y[i]+k1[i][1]/2);
    vx1[i]=(vx[i]+k1[i][2]/2);
    vy1[i]=(vy[i]+k1[i][3]/2);
  }
  
  double k2[n][4];  // Se calculan los valores de k2
  for (int i=0; i<n; i++){
    dat  = Fk(i,x1,y1,vx1,vy1,n,m);
    k2[i][0] = h*dat[0];
    k2[i][1] = h*dat[1];
    k2[i][2] = h*dat[2];
    k2[i][3] = h*dat[3];
  }
  double x2[n], y2[n], vx2[n], vy2[n]; //Listas que contendrán los valores asociados al término k2 del método RK4
  for (int i=0; i<n; i++){
    x2[i] = (x[i]+k2[i][0]/2);
    y2[i] = (y[i]+k2[i][1]/2);
    vx2[i] = (vx[i]+k2[i][2]/2);
    vy2[i] = (vy[i]+k2[i][3]/2);
  }
  
  double k3[n][4]; // Se calculan los valores de k3
  for(int i=0; i<n; i++){
    dat = Fk(i,x2,y2,vx2,vy2,n,m);
    k3[i][0] = h*dat[0];
    k3[i][1] = h*dat[1];
    k3[i][2] = h*dat[2];
    k3[i][3] = h*dat[3];
  }
  double x3[n], y3[n], vx3[n], vy3[n]; //Listas que contendrán los valores asociados al término k3 del método RK4
  for (int i=0; i<n; i++){
    x3[i] = (x[i] + k3[i][0]);
    y3[i] = (y[i] + k3[i][1]);
    vx3[i] = (vx[i] + k3[i][2]);
    vy3[i] = (vy[i] + k3[i][3]);
  }

  double k4[n][4];  // Se calculan los valores de k4
  for(int i=0; i<n; i++){
    dat = Fk(i,x3,y3,vx3,vy3,n,m);
    k4[i][0] = h*dat[0];
    k4[i][1] = h*dat[1];
    k4[i][2] = h*dat[2];
    k4[i][3] = h*dat[3];
  }
  double xf[n], yf[n], vxf[n], vyf[n];  // Listas que contendrán los valores finales del intervalo
  for(int i=0; i<n; i++){
    xf[i] =  x[i] + (k1[i][0]+2*k2[i][0]+2*k3[i][0]+k4[i][0])/6;
    yf[i] =  y[i] + (k1[i][1]+2*k2[i][1]+2*k3[i][1]+k4[i][1])/6;
    vxf[i] = vx[i] + (k1[i][2]+2*k2[i][2]+2*k3[i][2]+k4[i][2])/6;
    vyf[i] = vy[i] + (k1[i][3]+2*k2[i][3]+2*k3[i][3]+k4[i][3])/6;
  }
  
  double** array2D=0;      // Se guardan los resultados
  array2D = new double*[4];
  for (int i=0; i<4; i++)
    {
    array2D[i] = new double[n];
    }
  
  
  for(int i=0; i<n; i++){
    array2D[0][i] = xf[i];
    array2D[1][i] = yf[i];
    array2D[2][i] = vxf[i];
    array2D[3][i] = vyf[i];
  }
  return array2D;   // Se devuelven los resultados
}

// Funciones para hallar el maximo y minimo de un vector
double maximo(double x[]){

  int n = sizeof(x)/sizeof(x[0]);   // Tamaño del vector x.
  double max=x[0];
  for(int i=1; i<n; i++){
    if (x[i]>max){ // Si x[i] es más grande que el máximo, entonces x[i] es el nuevo máximo.
      max=x[i];
    }
  }
  return max;
}

double minimo(double x[]){
  double min=x[0];
  int n = sizeof(x)/sizeof(x[0]);   // Tamaño del vector x
  for(int i=1; i<n; i++){
    if (x[i]<min){ // Si x[i] es más pequeño que el mínimo, entonces x[i] es el nuevo mínimo.
      min=x[i];
    }
  }
  return min;
}


// Definicion función principal.
int main(){
   
  
  
   // Valores inciales de posición y velocidad en listas de números reales.
   double x[4] = {0,-400E9, 150E9,206E9};
   double y[4] = {0,160E9,0,0};
   double vx[4] = {0,4000,0,0};
   double vy[4] = {-2000,0,30000,24000};

   // Lista con masas de los cuerpos (n en total) con números reales.
   double m[4] = {2E30,2E14,6E24,6E23};
   
   // Parametros y constantes.
   int N=60000;                    // Cantidad de pasos. Aumentar para añadir tiempo de animación.
   double h=10000;                  // Paso de tiempo. Modificar según se necesite.
   int n = sizeof(x)/sizeof(x[0]);   // Tamaño de x. Indica la cantidad de cuerpos.
   
   // Listas que contendrán los datos de posición y velocidad
   double X[N][n], Y[N][n], Vx[N][n], Vy[N][n];
   
   // Inicializando las listas de posición y velocidad
   for(int i=0; i<n; i++){
     X[0][i] = x[i];
     Y[0][i] = y[i];
     Vx[0][i] = vx[i];
     Vy[0][i] = vy[i];
   }
   
   
   //Se hallan los datos de posición, velocidad y aceleración para todos los N intervalos de este código.
   double** sig;
    
   for(int t=0; t<N-1; t++){
     sig = siguiente_valor(X[t],Y[t],Vx[t],Vy[t],n,h,m);
     
     for(int i=0; i<n; i++){
       X[t+1][i]=sig[0][i];
       Y[t+1][i]=sig[1][i];
       Vx[t+1][i]=sig[2][i];
       Vy[t+1][i]=sig[3][i];
       
     }
   }
   
   
   bool centrado=0;
   
   return 0;
}


