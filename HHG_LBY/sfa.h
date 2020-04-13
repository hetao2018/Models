#pragma once
#include <complex>
#include <array>
#include "parameters.h"
using namespace parameters;
using Array = std::array<double, nt>;
using ArrayC = std::array<complex, nt>;

class SFA_HHG {
 public:
  Array time, Ef, Af, intA1, intA2;
  ArrayC d;

  void prepare_pulse (Pulse& pulse) {
    for (int it = 0; it < nt; it ++) {
      time[it] = t0 + it * dt;
      Ef[it] = pulse.E (time[it]);
      Af[it] = pulse.A (time[it]);
      if (it == 0) {
        // Af[it] = 0.;
        intA1[it] = 0.;
        intA2[it] = 0.;
      } else {
        // Af[it] = Af[it-1] - Ef[it-1] * dt;
        intA1[it] = intA1[it-1] + Af[it-1] * dt;
        intA2[it] = intA2[it-1] + Af[it-1] * Af[it-1] * dt;
      }
    }
  }

  virtual double action_S (double ps, int itau, int it, int itp) {
    double tau = itau * dt;
    return (0.5 * ps * ps + Ip) * tau
      + ps * (intA1[it] - intA1[itp])
      + 0.5 * (intA2[it] - intA2[itp]);
  }

  void calc_hhg_dx () {
    const double epsilon = 0.1;
    for (int it = 0; it < nt; it ++) {
      d[it] = 0.;
      // itau loop starts with 1, otherwise tau=0 renders ps singular
      for (int itau = 1; itau <= it; itau ++) {
        double tau = itau * dt;
        int itp = it - itau; // ionization time
        double ps = - (intA1[it] - intA1[itp]) / tau;
        double S = action_S (ps, itau, it, itp);
        if (it % 1000 == 0 && itau % 1000 == 0){printf("dipoley: %f \n ", S);}
        d[it] += conj (dipole (ps + Af[it])) * dipole (ps + Af[itp]) * Ef[itp]
          * exp (- I * S) * pow (M_PI / (epsilon + I * tau / 2.), 1.5);
      }
      d[it] *= I * (-dt);
      if (it % 1000 == 0)
        std::cout << "progress: " << 1.0 * it / nt << std::endl;
    }
  }

};
