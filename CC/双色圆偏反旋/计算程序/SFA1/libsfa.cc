#include <cstdio>
#include <complex>

#define complex std::complex<double>
#define I complex(0., 1.)
const double Ip = 0.5724;

  // 定义计算矩阵元函数
inline complex dipole_H_1s_S (double k_S, double k_P) {
  double dnom = (k_S*k_S + k_P*k_P + 2*Ip);
  // pow(2., 3.5) = 11.3137
  return - I * 11.3137 * pow(2.*Ip, 1.25) * k_S / (M_PI * dnom * dnom * dnom);
}

inline complex dipole_H_1s_P (double k_S, double k_P) {
  double dnom = (k_S*k_S + k_P*k_P + 2*Ip);
  // pow(2., 3.5) = 11.3137
  return - I * 11.3137 * pow(2.*Ip, 1.25) * k_P / (M_PI * dnom * dnom * dnom);
}

struct SFA_1d {
  long nt;
  double dt;
  double *t_grid, *ES, *AS, *alphaS, *EP, *AP, *alphaP;

  SFA_1d (long nt, double *t_grid, double *ES, double *AS, double *alphaS, double *EP, double *AP, double *alphaP):
  nt(nt), t_grid(t_grid), ES(ES), AS(AS), alphaS(alphaS), EP(EP), AP(AP), alphaP(alphaP) {
    dt = t_grid[1] - t_grid[0];
  }

  //S方向稳相点计算
  inline double k_sp_S (long it1, long it2) {
    return - (alphaS[it1] - alphaS[it2]) / ((it1-it2) * dt);   //saddle point
  }                   

  //P方向稳相点计算
  inline double k_sp_P (long it1, long it2) {
    return - (alphaP[it1] - alphaP[it2]) / ((it1-it2) * dt);   //saddle point	
  }


  // 2nd-order derivative of the C-C dipole
  void calc_dipole_CC_d2 (double* dipole_d_S_d2, double* dipole_d_P_d2) {

    for (long it = 0; it < nt; it ++) {
      complex d = 0.;
      complex term2S = 0.;
      complex term2P = 0.;
      printf("%ld of %ld\n", it, nt);
      for (long it1 = 0; it1 < it; it1 ++) {
        complex term1 = 0.;
        for (long itau = 1; itau <it - it1; itau ++) {
          double ks_S = k_sp_S (it1+itau,it1);
          double ks_P = k_sp_P (it1+itau,it1);
          double S = 0.;
          for (long it2 = it1; it2 <it1 + itau; it2 ++) {
            double k_S = ks_S + AS[it2];
            double k_P = ks_P + AP[it2];
            S += (0.5 * (k_S*k_S + k_P*k_P) + Ip) * dt;
          }	  
          term1 += pow (2.*M_PI / (itau * dt * I), 1.5) * exp (I * S)
  			*(conj (dipole_H_1s_S (ks_S+AS[it1], ks_P+AP[it1])) * ES[it1]+conj (dipole_H_1s_P (ks_S+AS[it1], ks_P+AP[it1]))*EP[it1])
  			*(dipole_H_1s_S (ks_S+AS[it1+itau], ks_P+AP[it1+itau])*ES[it1+itau]+dipole_H_1s_P (ks_S+AS[it1+itau], ks_P+AP[it1+itau])*EP[it1+itau])*dt;
        }
  	d += term1*dt;
  	double ks_S = k_sp_S(it, it1);
  	double ks_P = k_sp_P(it, it1);
  	double S = 0.;
  	for (long it2 = it1; it2 < it; it2++){
  	  double k_S = ks_S + AS[it2];
  	  double k_P = ks_P + AP[it2];
  	  S += (0.5*(k_S*k_S+k_P*k_P)+Ip)*dt;
  	}
  	term2S += pow (2.*M_PI/((it-it1)*dt*I), 1.5)*exp(I*S)
  	  *(conj (dipole_H_1s_S(ks_S+AS[it1], ks_P+AP[it1]))*ES[it1] + conj (dipole_H_1s_P(ks_S+AS[it1], ks_P+AP[it1]))*EP[it1])
  	  *(dipole_H_1s_S(ks_S+AS[it], ks_P+AP[it])*ES[it] + dipole_H_1s_P (ks_S+AS[it], ks_P+AP[it])*EP[it])*(ks_S+AS[it])*dt;
	
  	term2P += pow (2.*M_PI/((it-it1)*dt*I), 1.5)*exp(I*S)
  	  *(conj (dipole_H_1s_S(ks_S+AS[it1], ks_P+AP[it1]))*ES[it1] + conj(dipole_H_1s_P(ks_S+AS[it1], ks_P+AP[it1]))*EP[it1])
  	  *(dipole_H_1s_S(ks_S+AS[it], ks_P+AP[it])*ES[it] + dipole_H_1s_P(ks_S+AS[it], ks_P+AP[it])*EP[it])*(ks_P+AP[it])*dt;
      }
      complex d_S1 = d*ES[it];
      complex d_P1 = d*EP[it];
      dipole_d_S_d2[it] = 2.*real (d_S1) - 2.*real(term2S);
      dipole_d_P_d2[it] = 2.*real (d_P1) - 2.*real(term2P);
    }
  }

   // void calc_wt (double* dipole_d_d2) {
   //    FILE* file = fopen ("result.dat", "w");
   //    for (long it = 0; it < nt; it++) {
   //      printf("%ld of %ld\n", it, nt);
   //      for (long it1 = 0; it1 < it; it1 ++) {
   //        complex term1 = 0.;
   //        for (long itau = 1; itau < it - it1; itau ++) {
   //          double ks_S = k_sp_S (it1+itau,it1);
   //          double ks_P = k_sp_P (it1+itau,it1);
   //          double S = 0.;
   //          for (long it2 = it1; it2 < it1 + itau; it2 ++) {
   //            double k_S = ks_S + AS[it2];
   //            double k_P = ks_P + AP[it2];
   //            S += (0.5 * (k_S*k_S + k_P*k_P) + Ip) * dt;
   //          }  
   //           term1 += pow (2.*M_PI / (itau * dt * I), 1.5) * exp (I * S)
   // 	       * (conj (dipole_H_1s_S (ks_S+AS[it1], ks_P+AP[it1])) *ES[it1]+conj (dipole_H_1s_P (ks_S+AS[it1], ks_P+AP[it1]))*EP[it1])
   // 	       * (dipole_H_1s_S (ks_S+AS[it1+itau], ks_P+AP[it1+itau])*ES[it1+itau]+dipole_H_1s_P (ks_S+AS[it1+itau], ks_P+AP[it1+itau])*EP[it1+itau])*dt;
   //        }
   // 	  dipole_d_d2[it1] = 2.*real(term1);
   // 	  fprintf (file, "%le ", dipole_d_d2[it1]);
   //      }
   //      fprintf (file, "\n");
   //    }
   //     fclose (file);    
   // }

    void calc_Nt (double* dipole_d_d2) {
    for (long it = 0; it < nt; it ++) {
      complex d = 0.;
      printf("%ld of %ld\n", it, nt);
      for (long it1 = 0; it1 < it; it1 ++) {
        complex term1 = 0.;
        for (long itau = 1; itau <it - it1; itau ++) {
          double ks_S = k_sp_S (it1+itau,it1);
          double ks_P = k_sp_P (it1+itau,it1);
          double S = 0.;
          for (long it2 = it1; it2 <it1 + itau; it2 ++) {
            double k_S = ks_S + AS[it2];
            double k_P = ks_P + AP[it2];
            S += (0.5 * (k_S*k_S + k_P*k_P) + Ip) * dt;
          }	  
          term1 += pow (2.*M_PI / (itau * dt * I), 1.5) * exp (I * S)
  			*(conj (dipole_H_1s_S (ks_S+AS[it1], ks_P+AP[it1])) * ES[it1]+conj (dipole_H_1s_P (ks_S+AS[it1], ks_P+AP[it1]))*EP[it1])
  			*(dipole_H_1s_S (ks_S+AS[it1+itau], ks_P+AP[it1+itau])*ES[it1+itau]+dipole_H_1s_P (ks_S+AS[it1+itau], ks_P+AP[it1+itau])*EP[it1+itau])*dt;
        }
  	d += term1*dt;
      }
      complex nt = 2.*d;
      dipole_d_d2[it] = real (nt) ;
    }
  }


  
};

extern "C" void calc_dipole_CC_d2 (long nt, double* t_grid,
                                double* ES, double* AS, double* alphaS,
                                double* EP, double* AP, double* alphaP,
                                double* dipole_d_S_d2, double* dipole_d_P_d2) {
  SFA_1d calc_core (nt, t_grid, ES, AS, alphaS, EP, AP, alphaP);
  calc_core.calc_dipole_CC_d2 (dipole_d_S_d2, dipole_d_P_d2);
}


// extern "C" void calc_wt (long nt, double* t_grid,
//                          double* ES, double* AS, double* alphaS,
//                          double* EP, double* AP, double* alphaP,
//                          double* dipole_d_d2) {
//   SFA_1d calc_core (nt, t_grid, ES, AS, alphaS, EP, AP, alphaP);
//   calc_core.calc_wt (dipole_d_d2);
// } 

// extern "C" void calc_Nt (long nt, double* t_grid,
//                          double* ES, double* AS, double* alphaS,
//                          double* EP, double* AP, double* alphaP,
//                          double* dipole_d_d2) {
//   SFA_1d calc_core (nt, t_grid, ES, AS, alphaS, EP, AP, alphaP);
//   calc_core.calc_Nt (dipole_d_d2);
// } 

