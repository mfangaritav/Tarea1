#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
  C adaptation from C++ code written by Richard J. Gonsalves.
*/ 

double L = 4.0;                    // length of shock tube
double gama = 1.4;               // ratio of specific heats
int N = 5000;                     // number of grid points

double CFL = 1;                // Courant-Friedrichs-Lewy number
//double nu = 0.0;                 // artificial viscosity coefficient

//double vol = NULL;                     // for Roe solver


/*void allocate() {
  int j;

  U = malloc(N * sizeof(double *));
  newU = malloc(N * sizeof(double *));
  F = malloc(N * sizeof(double *));
  vol = malloc(N * sizeof(double *));
  for (j = 0; j < N; j++) {
    U[j] = malloc(3 * sizeof(double));
    newU[j] = malloc(3 * sizeof(double));
    F[j] = malloc(3 * sizeof(double));
  }
}*/
double cMax(double U[N][3]) {
    double uMax = 0;
    double rho, u, p, c;
    int i;
    for (i = 0; i < N; i++) {
      if (U[i][0] == 0)
	continue;

      rho = U[i][0];
      u = U[i][1] / rho;
      p = (U[i][2] - rho * u * u / 2) * (gama - 1);
      c = sqrt(gama * fabs(p) / rho);

      if (uMax < (c + fabs(u)))
	uMax = c + fabs(u);
    }    
    return uMax;
}

void initialize(double U[N][3]) {
  int j;
  double rho,p,u,e;

  //allocate();
  //h = 1.0 * L / (N - 1);
  for (j = 0; j < N; j++) {
    rho = 1; p = 1; u = 0;
    if (j > N / 2){
      rho = 0.125;
      p = 0.1;
    }
    e = p / (gama - 1) + rho * u * u / 2;
    U[j][0] = rho;
    U[j][1] = rho * u;
    U[j][2] = e;
    //vol[j] = 1;
  }
  //tau = CFL * h / cMax(U);
  //step = 0;
}
void obtenerUsol(double usol[N][3], double U[N][3]){
  int j;
  for(j=0;j<N;j++){
    usol[j][0] = U[j][0];
    usol[j][1] = U[j][1] / U[j][0];
    usol[j][2] = (U[j][2] - U[j][1] * U[j][1] / U[j][0] / 2)* (gama - 1.0);
  }
}

void multMatrices(double M1[3][3], double M2[3][3],double R[3][3]){
  int i,j;
  //printf("este es el resultado\n");
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      R[i][j]=M1[i][0]*M2[0][j]+M1[i][1]*M2[1][j]+M1[i][2]*M2[2][j];
      //printf("%f  ",R[i][j]);
      //printf("M1 es %f en %d %d\n",M1[i][j], i,j);
      //printf("M2 es %f en %d %d\n",M2[i][j], i,j);
    }
    //printf("\n");
  }
}
void matrizVector(double A[3][3],double v[3], double r[3]){
  int i;
  for(i=0;i<3;i++){
    r[i]=A[i][0]*v[0]+A[i][1]*v[1]+A[i][2]*v[2];
  }

}

int main(void)
{ 
    double h;                        // lattice spacing
    double tau;                      // time step
    printf("empece\n");

    double tMax=1.0;
    char filename[80];
    sprintf(filename,"UpwindGodunov");
    double U[N][3];
    initialize(U);
    double gamma=1.4;
    double t = 0.0;
    int step = 0;
    int plot = 0;
    int i;
    FILE *out;
    char filename_tmp[1024];
    int j;
    double rho_avg = 0.0, u_avg = 0.0, e_avg = 0.0, P_avg = 0.0;
    double rho, u, e, P;
    h = 1.0 * L / (N - 1);
    tau = CFL * h / cMax(U);
    sprintf(filename_tmp, "%s_step_%d.dat", filename, plot);
    if(!(out = fopen(filename_tmp, "w"))){
      fprintf(stderr, "problem opening file %s\n", filename);
      exit(1);
    }
    // write solution in plot files and print       
    for (j = 0; j < N; j++) {
      rho = U[j][0];
      u = U[j][1] / U[j][0];
      e = U[j][2];
      P = (U[j][2] - U[j][1] * U[j][1] / U[j][0] / 2)* (gama - 1.0);
        
      rho_avg += rho;
      u_avg += u;
      e_avg += e;
      P_avg += P;
      fprintf(out, "%d\t%f\t%f\t%f\t%f\n", j, rho, u, e, P);
    }

    fclose(out);    
    
    double usol[N][3];
    obtenerUsol(usol,U);
    plot=1;
    while (t < tMax+0.1) {
      
      tau = CFL * h / cMax(U);
      double htot[N];
      for(i=0;i<N;i++){
        htot[i]=(gamma/(gamma-1))*(usol[i][2]/usol[i][0])+0.5*pow(usol[i][1],2);
      }
      double Phi[N-1][3];
      for(j=0;j<N-1;j++){
        double r=sqrt(usol[j+1][0]/usol[j][0]);
        double rm=r*usol[j][0];
        double um=(r*usol[j+1][1]+usol[j][1])/(r+1);
        double hm=(r*htot[j+1]+htot[j])/(r+1);
        double am=sqrt((gamma-1)*(hm-0.5*um*um));
        double alfa1=(gamma-1)*um*um/(2*am*am);
        double alfa2=(gamma-1)/(am*am);
        double Udif[3];
        for(i=0;i<3;i++){
         Udif[i]=U[j+1][i]-U[j][i];
        }
        double Pinv[3][3];
        Pinv[0][0]=0.5*(alfa1+um/am);
        Pinv[0][1]=-0.5*(alfa2*um+1/am);
        Pinv[0][2]=alfa2/2;
        Pinv[1][0]=1-alfa1;
        Pinv[1][1]=alfa2*um;
        Pinv[1][2]=-alfa2;
        Pinv[2][0]=0.5*(alfa1-um/am);
        Pinv[2][1]=-0.5*(alfa2*um-1/am);
        Pinv[2][2]=alfa2/2;
        double P[3][3];
        P[0][0]=1.0;
        P[0][1]=1.0;
        P[0][2]=1.0;
        P[1][0]=um-am;
        P[1][1]=um;
        P[1][2]=um+am;
        P[2][0]=hm-am*um;
        P[2][1]=0.5*um*am;
        P[2][2]=hm+am*um;
        double L[3][3];
        L[0][0]=fabs(um-am);
        L[0][1]=0;
        L[0][2]=0;
        L[1][0]=0;
        L[1][1]=fabs(um);
        L[1][2]=0;
        L[2][0]=0;
        L[2][1]=0;
        L[2][2]=fabs(um+am);
        double R[3][3];
        multMatrices(P,L,R);
        double R1[3][3];
        multMatrices(R,Pinv,R1);
        double Phi1[3];
        matrizVector(R1,Udif,Phi1);
        for(i=0;i<3;i++){
           Phi[j][i]=Phi1[i];
        }
      }
      double F[N][3];
      for(j=0;j<N;j++){
        F[j][0]=U[j][1];
        double uloc1=U[j][1]/U[j][0];
        double rhou2=U[j][1]*uloc1;
        double press1=(gamma-1)*(U[j][2]-0.5*rhou2);
        F[j][1]=rhou2+press1;
        F[j][2]=(U[j][2]+press1)*uloc1;
      }
      for(j=0;j<N-1;j++){
        for(i=0;i<3;i++){
          Phi[j][i]*=-0.5;
          Phi[j][i]+=0.5*(F[j][i]+F[j+1][i]);
        }
      }
      for(j=1;j<N-1;j++){
        for(i=0;i<3;i++){
          U[j][i]=U[j][i]-tau/h*(Phi[j][i]-Phi[j-1][i]);
          //printf("%f  ",w[i][j]);
        }
      }
      if(t>=plot*0.2){
         sprintf(filename_tmp, "%s_step_%d.dat", filename, plot);
         out = fopen(filename_tmp, "w");
         plot++;
         for(i=0;i<N;i++){
           usol[i][0]=U[i][0];
           //printf("%f  ",usol[0][i]);
           usol[i][1]=U[i][1]/U[i][0];
           //printf("%f  ",usol[1][i]);
           usol[i][2]=(gamma-1)*(U[i][2]-0.5*U[i][1]*usol[i][1]);
           //printf("%f  \n",usol[2][i]);
           fprintf(out, "%d\t%f\t%f\t%f\t%f\n", i, usol[i][0], usol[i][1], U[i][2], usol[i][2]);
        }
      }      
      t += tau;
      obtenerUsol(usol,U);
      
      
    }
    
}


