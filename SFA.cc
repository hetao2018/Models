#include <cstdio>

// define the parameter of fundamental electric field
#define w0 0.057
//#define E0 0.08
#define A0 1.4 //A0=E0/w0 //1.4
#define nc0 24

// define the parameter of harmonic electric field
#define w1 0.114
//#define E1 0.046 //conversion efficiency ~30%
#define A1 0.4 //0.4
#define nc1 48

// define the relative relationship of fundamental&second electric field
// #define Delta 50.
// #define Phi 1.

//global parameter
#define nt 28001
#define t0 0.
#define t1 20000
#define tc 10000
//const double dt = (t1-t0)/(nt-1) //1fs=42a.u. 0.6a.u.~0.014fs

#define Ip 0.58

#include <complex>
#define complex std::complex<double>


//Creation_Time_Grid
void Creation_Time_Grid (double* time){

  double dt = (t1-t0)/(nt-1);
  //time axis
  for (int it = 0; it < nt; it ++) {
    time[it] = t0 + it * dt;
  }

}
//Creation_Time_Grid

//Creation_Potential
void Creation_Potential (double* time, double* AfP, double* AsP, double* AfS, double* AsS, double Delta){

  //# fundamental vector potential
  double duration0 = 2. * M_PI / w0 * nc0;
  //define the potential of fundamental electric field
  for (int it = 0; it < nt; it ++) {
     auto phi = w0 * (time[it] - tc + duration0 / 2.);
    if (phi >= 0 && phi < w0 * duration0) {
       auto e = sin(phi / (2.*nc0));
      AfP[it] = A0 * e * e * sin(phi);
      AfS[it] = A0 * e * e * cos(phi);
    } else {
      AfP[it] = 0.;
      AfS[it] = 0.;
    }
  }
  //define the potential of fundamental electric field

  //# second harmonic vector potential
  double duration1 = 2. * M_PI / w1 * nc1;
  //define the potential of fundamental electric field
  for (int it = 0; it < nt; it ++) {
      auto phi = w1 * (time[it] - tc + duration1 / 2.);
    if (phi >= 0 && phi < w1 * duration1) {
       auto e = sin(phi / (2.*nc1));
      AsP[it] = A1 * e * e * sin(phi);
      AsS[it] = A1 * e * e * cos(phi);
    } else {
      AsP[it] = 0.;
      AsS[it] = 0.;
    }
  }
  //define the potential of fundamental electric field

}
//Creation_Potential

//Creation_Electric_Field
void Creation_Electric_Field (double* time, double* EfP, double* EsP, double* EfS, double* EsS,double Delta){

  //# fundamental electric field
  double duration0 = 2. * M_PI / w0 * nc0;
  //define fundamental electric field analytically
    for (int it = 0; it < nt; it ++){
      auto phi = w0 * (time[it] - tc + duration0 / 2.);
      if(phi >= 0 && phi < w0 * duration0){
        EfP[it] = -1 * A0 * w0 / (2*nc0) * sin(phi / nc0) * sin(phi)
                     - A0 * w0 * sin(phi / (2*nc0))
                     * sin(phi / (2*nc0)) * cos(phi);
        EfS[it] = -1 * A0 * w0 / (2*nc0) * sin(phi / nc0) * cos(phi)
                     + A0 * w0 * sin(phi / (2*nc0))
                     * sin(phi / (2*nc0)) * sin(phi);
      }else{
        EfP[it] = 0.;
        EfS[it] = 0.;
      }
    }
    //define fundamental electric field analytically

    //# second electric field
    double duration1 = 2. * M_PI / w1 * nc1;
    //define second electric field analytically
    for (int it = 0; it < nt; it ++){
       auto phi = w1 * (time[it] - tc + duration1 / 2.);
      if(phi >= 0 && phi < w1 * duration1){
        EsP[it] = -1 * A1 * w1 / (2*nc1) * sin(phi / nc1) * sin(phi)
                     - A1 * w1 * sin(phi / (2*nc1))
                     * sin(phi / (2*nc1)) * cos(phi);
        EsS[it] = -1 * A1 * w1 / (2*nc1) * sin(phi / nc1) * cos(phi)
                     + A1 * w1 * sin(phi / (2*nc1))
                     * sin(phi / (2*nc1)) * sin(phi);
      }else{
        EsP[it] = 0.;
        EsS[it] = 0.;
      }
    }
    //define second electric field analytically

}
//Creation_Electric_Field

//define the integrated Potential A to void artifical slope
// void Creation_Alpha(double* time, double* Alphaf, double* Alphas, double Delta){
//   //# fundamental Alphaf
//   double duration0 = 2. * M_PI / w0 * nc0;
//   //define fundamental Alphaf analytically
//   for (int it = 0; it < nt; it ++){
//     if(time[it] >= -duration0 / 2. && time[it] < duration0 / 2.){
//       Alphaf[it] =  A0 * w0 / (2. * w0 * (-1. + nc0 * nc0))
//       * (-1. + cos(w0 * (time[it]-duration0/2.)) * (1. - nc0*nc0 + nc0*nc0*cos(w0 * (time[it]-duration0/2.)/nc0)) +
//       nc0 * sin(w0 * (time[it]-duration0/2.)) * sin(w0 * (time[it]-duration0/2.)/nc0));
//     }else{
//       Alphaf[it] = 0.;
//     }
//   }
//
//   //second harmonic Alphas
//   double duration1 = 2. * M_PI / w1 * nc1;
//   //define second harmonic Alphas analytically
//   for (int it = 0; it < nt; it ++){
//     if(time[it] >= -duration1 / 2. && time[it] < duration1 / 2.){
//       Alphas[it] = A1 * w1 / (2. * w1 * (-1. + nc1 * nc1))
//       * (-1. + cos(w1 * (time[it]-Delta-duration1/2.)) * (1. - nc1*nc1 + nc1*nc1*cos(w1 * (time[it]-Delta-duration1/2.)/nc1)) +
//       nc1 * sin(w1 * (time[it]-Delta-duration1/2.)) * sin(w1 * (time[it]-Delta-duration1/2.)/nc1));
//     }else{
//       Alphas[it] = 0.;
//     }
//   }
//
// }
//define the integrated Potential A to void artifical slope

//define Trajectory function
void Trajectory (double* time, double* AfP, double* AsP, double* AfS, double* AsS, double* EfP, double* EsP,double* EfS, 
                 double* EsS, double Delta){
  Creation_Time_Grid(time);
  Creation_Potential(time, AfP, AsP, AfS, AsS, Delta);
  Creation_Electric_Field(time, EfP, EsP, EfS, EsS, Delta);
}
//define Trajectory function

//coordination transform
void Transform (double* LabX, double* LabY, double* PolfP, double* PolsP, double* PolfS, double* PolsS){
  for(int it = 0; it < nt; it ++){
    LabX[it] = PolfP[it]+PolsP[it];
    LabY[it] = PolfS[it]-PolsS[it];
  }
}

// prepare pulse & Ain
void Prepare_Ain (double* intA1x, double* intA1y, double* intA2x, double* intA2y,double* Ax, double* Ay){
     double dt = (t1-t0)/(nt-1); 
   for (int it = 0; it < nt; it++){
    if (it == 0){
     intA1x[it] = 0.;
     intA1y[it] = 0.;
     intA2x[it] = 0.;
     intA2y[it] = 0.;
   }else {
        intA1x[it] = intA1x[it-1] + Ax[it-1] * dt;
        intA1y[it] = intA1y[it-1] + Ay[it-1] * dt;
        intA2x[it] = intA2x[it-1] + Ax[it-1] * Ax[it-1] * dt;
        intA2y[it] = intA2y[it-1] + Ay[it-1] * Ay[it-1] * dt;
      }  
   }
  }

// void action_S (double psx, double psy, int itau, int it, int itp) {
//     return (0.5 * (psx * psx + psy * psy + Ip)) * itau * dt
//             + psx * (intA1x[it] - intA1x[itp]) + psy * (intA1y[it] - intA1y[itp])
//             + 0.5 * (intA2x[it] - intA2x[itp]) + 0.5 * (intA2y[it] - intA2y[itp]) ;
// }



//coordination transform
//定义计算S方向矩阵元函数
inline complex dipolex(double px, double py, double Ax, double Ay){
  double dnom = (px + Ax) * (px + Ax) + (py + Ay)*(py + Ay) + 2 * Ip;
  return -complex(0., 1.) * 11.3137 * pow(2.*Ip, 1.25) * (px + Ax)/(M_PI * dnom * dnom * dnom); 
}
//定义计算P方向矩阵元函数
inline complex dipoley(double px, double py, double Ax, double Ay){
  double dnom = (px + Ax) * (px + Ax) + (py + Ay) * (py + Ay) + 2 * Ip;
  return -complex(0., 1.) * 11.3137 * pow(2.*Ip, 1.25) * (py + Ay)/(M_PI * dnom * dnom * dnom);
}

//calculate high-harmonic generation
void HHG (double* Ax, double* Ay, double* Ex, double* Ey, double* Dx, double* Dy,
           double* intA1x, double* intA1y, double* intA2x, double* intA2y){
  complex HHGx[nt];
  complex HHGy[nt];
  const double epsilon = 0.1;
  double dt = (t1-t0)/(nt-1);
  for (int it = 0; it < nt; it ++){
    HHGx[it] = 0.;
    HHGy[it] = 0.;
    for (int itau = 1; itau <= it; itau ++){
          int itp = it - itau; //ionization instant
      double tau = itau * dt;
      //double psx = (Alphax[itp] - Alphax[it]) / (itau * dt);
      //double psy = (Alphay[itp] - Alphay[it]) / (itau * dt);
      double psx = - (intA1x[it] - intA1x[itp]) / tau;
      double psy = - (intA1y[it] - intA1y[itp]) / tau;
      double S = (0.5 * (psx * psx + psy * psy) + Ip) * tau
            + psx * (intA1x[it] - intA1x[itp]) + psy * (intA1y[it] - intA1y[itp])
            + 0.5 * (intA2x[it] - intA2x[itp]) + 0.5 * (intA2y[it] - intA2y[itp]) ;
      // S += (0.5*((psx + Ax[itp]) * (psx + Ax[itp]) + (psy + Ay[itp]) * (psy + Ay[itp])) + Ip) * (dt);
      if (it % 1000 == 0 && itau % 1000 == 0){printf("psx&psy: %f \n ", S);}
      HHGx[it] += conj(dipolex(psx, psy, Ax[it], Ay[it])) * (dipolex(psx,psy,Ax[itp],Ay[itp])*Ex[itp] + dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp])
      * exp(-complex(0., 1.) * S) * pow(M_PI / (epsilon + complex(0., 1.) * tau / 2.), 1.5);

      HHGy[it] += conj(dipoley(psx, psy, Ax[it], Ay[it])) * (dipolex(psx,psy,Ax[itp],Ax[itp])*Ex[itp] + dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp])
      * exp(-complex(0., 1.) * S) * pow(M_PI / (epsilon + complex(0., 1.) * tau / 2.), 1.5);
    }
    HHGx[it] *= -complex(0., 1.) * dt;
    HHGy[it] *= -complex(0., 1.) * dt;
    Dx[it] = HHGx[it].real();
    Dy[it] = HHGy[it].real();
    if (it % 1000 ==0) printf("HHGprogress: %f\n", 1.0 * it /nt);
  }

}
//calculate high-harmonic generation

int main(int argc, char const *argv[]){
  double* time = new double[nt];

  double* AfP = new double[nt];
  double* AsP = new double[nt];
  double* AfS = new double[nt];
  double* AsS = new double[nt];

  double* EfP = new double[nt];
  double* EsP = new double[nt];

  double* EfS = new double[nt];
  double* EsS = new double[nt];
    
  double* intA1x = new double[nt];
  double* intA1y = new double[nt];
  double* intA2x = new double[nt];
  double* intA2y = new double[nt];

  // double* Alphaf = new double[nt];
  // double* Alphas = new double[nt];

  //create the labatory coordinate
  double* Ax = new double[nt];
  double* Ay = new double[nt];

  double* Ex = new double[nt];
  double* Ey = new double[nt];

  // double* Alphax = new double[nt];
  // double* Alphay = new double[nt];

  double Delta; //the relative time between the fundamental and second harmonic pulse
  //double Phi;  //the angle between the fundamental and second harmonic wave

  //create dipole
  double* Dx = new double[nt];
  double* Dy = new double[nt];

  char filename[256];
  FILE* file0 = fopen("parameter.dat","w");

  for (int i = 0; i < 1; i ++){

  Delta = 41.3 * 1.35 / 1 * i;

  //Phi = 0.5 * M_PI / 40 * j;

  fprintf(file0, "%le \n", Delta);

  Trajectory(time, AfP, AsP, AfS, AsS, EfP, EsP, EfS, EsS, Delta);
  //Creation_Alpha(time, Alphaf, Alphas, Delta);
  Transform(Ax, Ay, AfP, AsP, AfS, AsS);
  Transform(Ex, Ey, EfP, EsP, EfS, EsS);
  Prepare_Ain(intA1x, intA1y, intA2x, intA2y, Ax, Ay);
  //Transform(Alphax, Alphay, Alphaf, Alphas, Phi);
  HHG(Ax, Ay, Ex, Ey, Dx, Dy, intA1x, intA1y, intA2x, intA2y);

  sprintf(filename, "HHG_%d.dat", i);
  FILE* file1 = fopen(filename, "w");
  for (int it = 0; it < nt; it ++){
    fprintf(file1, "%le %le %le %le %le %le %le %le %le \n",
               time[it], Ex[it], Ey[it], Ax[it], Ay[it], Dx[it], Dy[it], intA1x[it], intA1y[it]);
      }
    fclose(file1);
    printf("i = %d\n", i);
  }

  fclose(file0);

  // Delta = 0.;
  // Phi = M_PI / 2. * 1. * 1. / 2.;
  //
  // Trajectory(time, Af, As, Ef, Es, Delta);
  // //Creation_Alpha(time, Alphaf, Alphas, Delta);
  // Transform(Ax, Ay, Af, As, Phi);
  // Transform(Ex, Ey, Ef, Es, Phi);
  // //Transform(Alphax, Alphay, Alphaf, Alphas, Phi);
  // HHG(Ax, Ay, Ex, Ey, Dx, Dy);
  //
  // FILE* file = fopen("data2.dat", "w");
  // for (int it = 0; it < nt; it ++)
  //   fprintf(file, "%le %le %le %le %le %le %le\n",
  //           time[it], Ex[it], Ey[it], Ax[it], Ay[it], Dx[it], Dy[it]);
  // fclose(file);

  delete[] Dy;
  delete[] Dx;
  // delete[] Alphay;
  // delete[] Alphax;
  delete[] Ey;
  delete[] Ex;
  delete[] Ay;
  delete[] Ax;
  delete[] EsP;
  delete[] EfP;
  delete[] AsP;
  delete[] AfP;
  delete[] EsS;
  delete[] EfS;
  delete[] AsS;
  delete[] AfS;
  delete[] time;
  delete[] intA1x;
  delete[] intA1y;
  delete[] intA2x;
  delete[] intA2y;

  return 0;
}
