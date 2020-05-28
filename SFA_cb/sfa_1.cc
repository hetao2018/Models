#include <cstdio>
#include <math.h>

#define w 0.057
#define E1 0.08
#define E2 0.046
#define nc 15
#define t0 -900
#define t1 900
#define nt 1801
#define dt 1.   // dt = (t1-t0)/(nt-1)


#define Ip 0.58

#include <complex>
#define complex std::complex<double>

void Creation_Time_Grid (double* time){

  // double dt = (t1-t0)/(nt-1);
  //time axis
  for (int it = 0; it < nt; it ++) {
    time[it] = t0 + it * dt;
  }
}

void Electric_Field(double* time, double* EfP, double* EfS, double* EsP, double* EsS, double Delta){
    double duration = 2. * M_PI / w * nc;
        for (int it = 0; it < nt; it ++){
      if(time[it] >= -duration/2 && time[it] <= duration/2){
        EfP[it] = sin(M_PI/4) * E1 * sin(w*(time[it]-duration/2)/(2*nc)) * sin(w*(time[it]-duration/2)/(2*nc)) * cos(w*time[it]);
        EfS[it] = cos(M_PI/4) * E1 * sin(w*(time[it]-duration/2)/(2*nc)) * sin(w*(time[it]-duration/2)/(2*nc)) * sin(w*time[it]);
      }else{
        EfP[it] = 0.;
        EfS[it] = 0.;
      }
    }
    for (int it = 0; it < nt; it ++){
      if(time[it] >= -duration/2+Delta*41.34 && time[it] <= duration/2+Delta*41.34){
        EsP[it] = sin(M_PI/4) * E2 * sin(w*(time[it]-duration/2-Delta*41.34)/(2.*nc))* sin(w*(time[it]-duration/2-Delta*41.34)/(2.*nc)) * cos(2*w*(time[it]-Delta*41.34));
        EsS[it] = cos(M_PI/4) * E2 * sin(w*(time[it]-duration/2-Delta*41.34)/(2.*nc))* sin(w*(time[it]-duration/2-Delta*41.34)/(2.*nc)) * sin(2*w*(time[it]-Delta*41.34));
      }else{
        EsP[it] = 0.;
        EsS[it] = 0.;
      }
    }
    
}

void Potential(double* Ep, double* Es, double* Ap, double* As , double* EfP , double* EsP , double* EfS , double* EsS){
     for (int it = 0; it < nt; it ++){
      Ep[it] = EfP[it] + EsP[it];
      Es[it] = EfS[it] + EsS[it];             // co or counter rotating
     }
     Ap[0] = Ep[0];
     As[0] = Es[0];
     for (int it = 0; it < nt-1; it++){
     Ap[it+1] = Ap[it] + 0.5 * (Ep[it] + Ep[it+1]) * dt ;
     As[it+1] = As[it] + 0.5 * (Es[it] + Es[it+1]) * dt ;
     }

}

void Prepare_Ain (double* intA1x, double* intA1y, double* intA2x, double* intA2y,double* Ap, double* As){
    // double dt = (t1-t0)/(nt-1); 
   for (int it = 0; it < nt; it++){
    if (it == 0){
     intA1x[it] = 0.;
     intA1y[it] = 0.;
     intA2x[it] = 0.;
     intA2y[it] = 0.;
   }else {
        intA1x[it] = intA1x[it-1] + Ap[it-1] * dt;
        intA1y[it] = intA1y[it-1] + As[it-1] * dt;
        intA2x[it] = intA2x[it-1] + Ap[it-1] * Ap[it-1] * dt;
        intA2y[it] = intA2y[it-1] + As[it-1] * As[it-1] * dt;
      }  
   }
  }
  
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

void HHG (double* Ax, double* Ay, double* Ex, double* Ey, double* Dx, double* Dy,
           double* intA1x, double* intA1y, double* intA2x, double* intA2y){
  complex HHGx[nt];
  complex HHGy[nt];
  const double epsilon = 0.1;
  // double dt = (t1-t0)/(nt-1);
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
      // if (it % 1000 == 0 && itau % 1000 == 0){printf("psx&psy: %f \n ", S);}
      HHGx[it] += conj(dipolex(psx, psy, Ax[it], Ay[it])) * (dipolex(psx,psy,Ax[itp],Ay[itp])*Ex[itp] + dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp])
      * exp(-complex(0., 1.) * S) * pow(2.*M_PI / (complex(0., 1.) * tau), 1.5);

      HHGy[it] += conj(dipoley(psx, psy, Ax[it], Ay[it])) * (dipolex(psx,psy,Ax[itp],Ay[itp])*Ex[itp] + dipoley(psx,psy,Ax[itp],Ay[itp])*Ey[itp])
      * exp(-complex(0., 1.) * S) * pow(2.*M_PI / (complex(0., 1.) * tau), 1.5);
    }
    HHGx[it] *= -complex(0., 1.) * dt;
    HHGy[it] *= -complex(0., 1.) * dt;
    Dx[it] = HHGx[it].real();
    Dy[it] = HHGy[it].real();
    if (it % 100 ==0) printf("HHGprogress: %f\n", 1.0 * it /nt);
  }

}

int main(int argc, char const *argv[]){
  double* time = new double[nt];

  double* Ap = new double[nt];
  double* As = new double[nt];
  double* Ep = new double[nt];
  double* Es = new double[nt];

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
  double* Dx = new double[nt];
  double* Dy = new double[nt];
  double Delta;
  double k;

  char filename[256];
  FILE* file0 = fopen("parameter.dat","w");

  for (int i = 0; i < 1; i ++){
      for (int j = 0; j < 1; j ++){

          Delta = 41.34 * 1.35 / 40 * i;
          k = j * 0.1 + 0.001;
          //A1 = sqrt(k)*A0/2;
          //Phi = 0.5 * M_PI / 40 * j;

          fprintf(file0, "%le %le \n", Delta, k);
          
          Creation_Time_Grid (time);
          Electric_Field(time, EfP, EfS, EsP, EsS, Delta);
          Potential(Ep, Es, Ap, As, EfP, EsP, EfS, EsS);

          //Trajectory(time, AfP, AsP, AfS, AsS, EfP, EsP, EfS, EsS, Delta, A1);
          //Creation_Alpha(time, Alphaf, Alphas, Delta);
          //Transform(Ax, Ay, AfP, AsP, AfS, AsS);
          //Transform(Ex, Ey, EfP, EsP, EfS, EsS);
          Prepare_Ain(intA1x, intA1y, intA2x, intA2y, Ap, As);
          //Transform(Alphax, Alphay, Alphaf, Alphas, Phi);
          HHG(Ap, As, Ep, Es, Dx, Dy, intA1x, intA1y, intA2x, intA2y);

          sprintf(filename, "HHG_co_%d_%d.dat", i,j);
          FILE* file1 = fopen(filename, "w");
          for (int it = 0; it < nt; it ++){
            fprintf(file1, "%le %le %le %le %le %le %le %le %le \n",
                    time[it], Ep[it], Es[it], Ap[it], As[it], Dx[it], Dy[it], intA1x[it], intA1y[it]);
              }
            fclose(file1);
      }
    printf("i = %d\n", i);
  }

  fclose(file0);


  delete[] Dy;
  delete[] Dx;
  // delete[] Alphay;
  // delete[] Alphax;
  delete[] Ep;
  delete[] Es;
  delete[] Ap;
  delete[] As;
  delete[] EsP;
  delete[] EfP;
  //delete[] AsP;
  //delete[] AfP;
  delete[] EsS;
  delete[] EfS;
  //delete[] AsS;
  //delete[] AfS;
  delete[] time;
  delete[] intA1x;
  delete[] intA1y;
  delete[] intA2x;
  delete[] intA2y;

  return 0;
}











