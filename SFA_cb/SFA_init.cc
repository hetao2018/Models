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
#define dt 0.6 //1fs=42a.u. 0.6a.u.~0.014fs
#define nt 5001

#define Ip 0.5

#include <complex>
#define complex std::complex<double>


//Creation_Time_Grid
void Creation_Time_Grid (double* time){

  double t0 = - dt * (nt - 1) / 2;
  //time axis
  for (int it = 0; it < nt; it ++) {
    time[it] = t0 + it * dt;
  }

}
//Creation_Time_Grid

//Creation_Potential
void Creation_Potential (double* time, double* Af, double* As, double Delta){

  //# fundamental vector potential
  double duration0 = 2. * M_PI / w0 * nc0;
  //define the potential of fundamental electric field
  for (int it = 0; it < nt; it ++) {
    if (time[it] >= -duration0 / 2. && time[it] < duration0 / 2.) {
      Af[it] = A0 * sin(w0 * (time[it]-duration0/2.) / (2.* nc0)) *
       sin(w0 * (time[it]-duration0/2.) / (2.* nc0)) * sin(w0 * (time[it]-duration0/2.));
    } else {
      Af[it] = 0.;
    }
  }
  //define the potential of fundamental electric field

  //# second harmonic vector potential
  double duration1 = 2. * M_PI / w1 * nc1;
  //define the potential of fundamental electric field
  for (int it = 0; it < nt; it ++) {
    if (time[it] >= -duration1 / 2. + Delta && time[it] < duration1 / 2. + Delta) {
      As[it] = A1 * sin(w1 * (time[it]-Delta-duration1/2.) / (2.* nc1)) *
       sin(w1 * (time[it]-Delta-duration1/2.) / (2.* nc1)) * sin(w1 * (time[it]-Delta-duration1/2.));
    } else {
      As[it] = 0.;
    }
  }
  //define the potential of fundamental electric field

}
//Creation_Potential

//Creation_Electric_Field
void Creation_Electric_Field (double* time, double* Ef, double* Es, double Delta){

  //# fundamental electric field
  double duration0 = 2. * M_PI / w0 * nc0;
  //define fundamental electric field analytically
    for (int it = 0; it < nt; it ++){
      if(time[it] >= -duration0 / 2. && time[it] < duration0 / 2. ){
        Ef[it] = -0.5 * A0 * w0 * sin(w0 * (time[it]-duration0/2.) / nc0) / nc0 * sin(w0 * (time[it]-duration0/2.))
                -A0 * sin(w0 * (time[it]-duration0/2.) / (2.* nc0)) * sin(w0 * (time[it]-duration0/2.) / (2.* nc0))
                * w0 * cos(w0 * (time[it]-duration0/2.));
      }else{
        Ef[it] = 0.;
      }
    }
    //define fundamental electric field analytically

    //# second electric field
    double duration1 = 2. * M_PI / w1 * nc1;
    //define second electric field analytically
    for (int it = 0; it < nt; it ++){
      if(time[it] >= -duration1 / 2. + Delta && time[it] < duration1 / 2. + Delta){
        Es[it] = -0.5 * A1 * w1 * sin(w1 * (time[it]-Delta-duration1/2.) / nc1) / nc1 * sin(w1 * (time[it]-Delta-duration1 /2.))
                -A1 * sin(w1 * (time[it]-Delta-duration1 /2.) / (2.* nc1)) * sin(w1 * (time[it]-Delta-duration1 /2.) / (2.* nc1))
                * w1 * cos(w1 * (time[it]-Delta-duration1 /2.));
      }else{
        Es[it] = 0.;
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
void Trajectory (double* time, double* Af, double* As, double* Ef, double* Es, double Delta){
  Creation_Time_Grid(time);
  Creation_Potential(time, Af, As, Delta);
  Creation_Electric_Field(time, Ef, Es, Delta);
}
//define Trajectory function

//coordination transform
void Transform (double* LabX, double* LabY, double* Polf, double* Pols, double Phi){
  for(int it = 0; it < nt; it ++){
    LabX[it] = Polf[it] * sin(Phi);
    LabY[it] = Polf[it] * cos(Phi) + Pols[it];
  }
}
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
void HHG (double* Ax, double* Ay, double* Ex, double* Ey, double* Dx, double* Dy){
  complex HHGx[nt];
  complex HHGy[nt];
  for (int it = 0; it < nt; it ++){
    HHGx[it] = 0.;
    HHGy[it] = 0.;
    complex S = 0.;
    double Aintx = 0.;
    double Ainty = 0.;
    for (int itau = 1; itau <= it; itau ++){
      int itp = it - itau; //ionization instant
      //double psx = (Alphax[itp] - Alphax[it]) / (itau * dt);
      //double psy = (Alphay[itp] - Alphay[it]) / (itau * dt);
      Aintx += Ax[itp] * dt;
      Ainty += Ay[itp] * dt;
      double psx = -Aintx / (itau * dt);
      double psy = -Ainty / (itau * dt);
      S += (0.5*((psx + Ax[itp]) * (psx + Ax[itp]) + (psy + Ay[itp]) * (psy + Ay[itp])) + Ip) * (dt);

      HHGx[it] += conj(dipolex(psx, psy, Ax[it], Ay[it])) * (dipolex(psx,psy,Ax[itp],Ay[itp])*Ex[itp] + dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp])
      * exp(-complex(0., 1.) * S) * pow(2.*M_PI/(itau * dt)*(1./complex(0.,1.)), 1.5);

      HHGy[it] += conj(dipoley(psx, psy, Ax[it], Ay[it])) * (dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp] + dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp])
      * exp(-complex(0., 1.) * S) * pow(2.*M_PI/(itau * dt)*(1./complex(0.,1.)), 1.5);
    }
    HHGx[it] *= -complex(0., 1.) * dt;
    HHGy[it] *= -complex(0., 1.) * dt;
    Dx[it] = HHGx[it].real();
    Dy[it] = HHGy[it].real();
    if (it % 100 ==0) printf("HHGprogress: %f\n", 1.0 * it /nt);
  }

}
//calculate high-harmonic generation

int main(int argc, char const *argv[]){

  double* time = new double[nt];

  double* Af = new double[nt];
  double* As = new double[nt];

  double* Ef = new double[nt];
  double* Es = new double[nt];

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
  double Phi;  //the angle between the fundamental and second harmonic wave

  //create dipole
  double* Dx = new double[nt];
  double* Dy = new double[nt];

  char filename[256];
  FILE* file0 = fopen("parameter.dat","w");

  for (int i = 0; i < 40; i ++){
    for (int j =0; j < 40; j ++){

      Delta = 41.3 * 1.35 / 40 * i;
      Phi = 0.5 * M_PI / 40 * j;

      fprintf(file0, "%le %le\n", Delta, Phi);

      Trajectory(time, Af, As, Ef, Es, Delta);
      //Creation_Alpha(time, Alphaf, Alphas, Delta);
      Transform(Ax, Ay, Af, As, Phi);
      Transform(Ex, Ey, Ef, Es, Phi);
      //Transform(Alphax, Alphay, Alphaf, Alphas, Phi);
      HHG(Ax, Ay, Ex, Ey, Dx, Dy);

      sprintf(filename, "HHG_%d_%d.dat", i, j);
      FILE* file1 = fopen(filename, "w");
      for (int it = 0; it < nt; it ++){
        fprintf(file1, "%le %le %le %le %le %le %le\n",
                   time[it], Ex[it], Ey[it], Ax[it], Ay[it], Dx[it], Dy[it]);
      }
      fclose(file1);
    }
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
  // delete[] Ay;
  // delete[] Ax;
  delete[] Es;
  delete[] Ef;
  delete[] As;
  delete[] Af;
  delete[] time;

  return 0;
}
