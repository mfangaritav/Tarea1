#include <stdio.h>
#include <math.h>


void RK4_step(double delta_t, double e, double t, double *x, double *v,double *x3,double *v3);
void leapfrog_step(double delta_t,double e, double t, double *x, double *v,double *x3, double *v3);
double deriv_x(double t, double x, double v);
double deriv_v(double t, double x, double v,double e);
double deriv_v3(double t, double x3,double x, double v, double e);
FILE *file1;
FILE *file2;


int main(){
  double x_RK4, x_LF, x3_RK4,x3_LF;
  double v_RK4, v_LF, v3_RK4,v3_LF;
  double t;
  double x3_CI;
  
  double T=2800;
  double delta_t=6E-3;
  int n_step = (int)(T/delta_t);
  double x_step = 0.0, e=1;
  double v_step = 0.0;
  int i,j;

  file1 = fopen("salida3LF.dat", "w");
  file2 = fopen("salida3RK4.dat", "w"); 

 
  t=0.0;
  x_RK4=0.425;
  v_RK4=0.0;
  x3_RK4=0.0;
  v3_RK4=0.0;
  x3_LF=0.0;
  v_LF=0.0;
  x_LF=0.46;
  v3_LF=0.0;
  x3_CI=-0.2;
  for(j=0;j<100;j++){
    x3_LF=x3_CI;
    v_LF=0.0;
    x_LF=0.46;
    v3_LF=0.0;
    x3_RK4=x3_CI;
    v_RK4=0.0;
    x_RK4=0.425;
    v3_RK4=0.0;
    for(i=0;i<n_step;i++){ 
      if(fabs(v_LF)<0.0005)   
          fprintf(file1,"%f %.15e %.15e %.15e %.15e \n", t, x_LF, v_LF, x3_LF, v3_LF);
      if(fabs(v_RK4)<0.0005)   
          fprintf(file2,"%f %.15e %.15e %.15e %.15e \n", t, x_LF, v_LF, x3_RK4, v3_RK4);
      RK4_step(delta_t,e, t, &x_RK4, &v_RK4, &x3_RK4, &v3_RK4);
      leapfrog_step(delta_t,e, t, &x_LF, &v_LF,&x3_LF,&v3_LF);
      t += delta_t;
    }
    x3_CI+=0.03;
  }
  fclose(file1);
  fclose(file2);
  return 0;
}

double deriv_x(double t, double x, double v){
  return v;
}
double deriv_v(double t, double x, double v, double e){
  return -2.0*x/pow(4*x*x+e,1.5);
}
double deriv_x3(double t, double x3, double v3){
  return v3;
}
double deriv_v3(double t, double x3,double x, double v, double e){
  return (x-x3)/pow(pow((x-x3),2)+pow(e,2)/4,1.5)-(x+x3)/pow(pow((x+x3),2)+pow(e,2)/4,1.5);
}

void leapfrog_step(double delta_t,double e, double t, double *x, double *v,double *x3,double *v3){
  double x_in,x3_in;
  double v_in,v3_in;
  x_in = *x;
  v_in = *v;
  x3_in=*x3;
  v3_in=*v3;

  x_in += 1/(2*(2-pow(2,1/3)))* v_in * delta_t;
  v_in += 1/(2-pow(2,1/3)) * deriv_v(t,x_in,v_in,e) * delta_t;
  x_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * v_in * delta_t;
  v_in += -pow(2,1/3)/(2-pow(2,1/3)) * deriv_v(t, x_in, v_in,e) * delta_t;
  x_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * v_in * delta_t;
  v_in += 1/(2-pow(2,1/3)) * deriv_v(t, x_in, v_in,e) * delta_t;
  x_in += 1/(2*(2-pow(2,1/3))) * v_in * delta_t;
  
  x3_in += 1/(2*(2-pow(2,1/3))) * v3_in * delta_t;
  v3_in += 1/(2-pow(2,1/3)) * deriv_v3(t,x3_in,x_in,v_in,e) * delta_t;
  x3_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * v3_in * delta_t;
  v3_in += -pow(2,1/3)/(2-pow(2,1/3)) * deriv_v3(t, x3_in,x_in, v_in,e) * delta_t;
  x3_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * v3_in * delta_t;
  v3_in += 1/(2-pow(2,1/3)) * deriv_v3(t,x3_in,x_in,v_in,e) * delta_t;
  x3_in += 1/(2*(2-pow(2,1/3))) * v3_in * delta_t;

  
  /*v_in += 1/(2*(2-pow(2,1/3)))* deriv_v(t,x_in,v_in,e) * delta_t;
  x_in += 1/(2-pow(2,1/3))  * v_in * delta_t;
  v_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * deriv_v(t, x_in, v_in,e) * delta_t;
  x_in += -pow(2,1/3)/(2-pow(2,1/3)) * v_in * delta_t;
  v_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * deriv_v(t, x_in, v_in,e) * delta_t;
  x_in += 1/(2-pow(2,1/3)) * v_in * delta_t;
  v_in += 1/(2*(2-pow(2,1/3)))* deriv_v(t, x_in, v_in,e) * delta_t;

  v3_in += 1/(2*(2-pow(2,1/3))) * deriv_v3(t,x3_in,x_in,v_in,e) * delta_t;
  x3_in += 1/(2-pow(2,1/3)) * v3_in * delta_t;
  v3_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * deriv_v3(t, x3_in,x_in, v_in,e) * delta_t;
  x3_in += -pow(2,1/3)/(2-pow(2,1/3)) * v3_in * delta_t;
  v3_in += (1-pow(2,1/3))/(2*(2-pow(2,1/3))) * deriv_v3(t,x3_in,x_in,v_in,e) * delta_t;
  x3_in += 1/(2-pow(2,1/3)) * v3_in * delta_t;
  v3_in += 1/(2*(2-pow(2,1/3))) * deriv_v3(t,x3_in,x_in,v_in,e) * delta_t;*/

  *x = x_in;
  *v = v_in;
  *x3=x3_in;
  *v3=v3_in;
}

void RK4_step(double delta_t, double e, double t, double *x, double *v,double *x3,double *v3){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double m1, m2, m3, m4;
  double n1, n2, n3, n4;
  double x_in;
  double v_in;
  double x3_in;
  double v3_in;
  x_in = *x;
  v_in = *v;
  x3_in = *x3;
  v3_in = *v3;

  k1 = deriv_x(t,x_in,v_in);  
  l1 = deriv_v(t,x_in,v_in,e);  

  k2 = deriv_x(t + delta_t*0.5, x_in + k1*delta_t*0.5, v_in + l1*delta_t*0.5);
  l2 = deriv_v(t + delta_t*0.5, x_in + k1*delta_t*0.5, v_in + l1*delta_t*0.5,e) ;

  k3 = deriv_x(t + delta_t*0.5, x_in + k2*delta_t*0.5, v_in + l2*delta_t*0.5);
  l3 = deriv_v(t + delta_t*0.5, x_in + k2*delta_t*0.5, v_in + l2*delta_t*0.5,e);

  k4 = deriv_x(t + delta_t, x_in + k3*delta_t, v_in + l3*delta_t);
  l4 = deriv_v(t + delta_t, x_in + k3*delta_t, v_in + l3*delta_t,e);

  x_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  v_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*delta_t;
  
  m1 = deriv_x3(t,x3_in,v3_in);  
  n1 = deriv_v3(t,x3_in,x_in,v_in,e);  

  m2 = deriv_x3(t + delta_t*0.5, x3_in + m1*delta_t*0.5, v3_in + n1*delta_t*0.5);
  n2 = deriv_v3(t + delta_t*0.5, x3_in + m1*delta_t*0.5, x_in + k1*delta_t*0.5, v_in + l1*delta_t*0.5,e) ;

  m3 = deriv_x3(t + delta_t*0.5, x3_in + m2*delta_t*0.5, v3_in + n2*delta_t*0.5);
  n3 = deriv_v3(t + delta_t*0.5, x3_in + m2*delta_t*0.5, x_in + k2*delta_t*0.5, v_in + l2*delta_t*0.5,e);

  m4 = deriv_x3(t + delta_t, x3_in + m3*delta_t, v3_in + n3*delta_t);
  n4 = deriv_v3(t + delta_t, x3_in + m3*delta_t, x_in + k3*delta_t, v_in + l3*delta_t,e);

  x3_in += (m1/6.0 + m2/3.0 + m3/3.0 + m4/6.0)*delta_t;
  v3_in += (n1/6.0 + n2/3.0 + n3/3.0 + n4/6.0)*delta_t;  

  *x = x_in;
  *v = v_in;
  *x3 = x3_in;
  *v3 = v3_in;
}
