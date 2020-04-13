#pragma once

namespace parameters {

#define complex std::complex<double>
const complex I = complex (0., 1.);

const double Ip = 0.58;
const double w0 = 0.057;
const double E0 = 0.08;
const int    nt = 28001;// 10001;
const double t0 = 0.;
  const double t1 = 20000;//4400*2*2;//2*M_PI/0.0008; // 2200;

// for THz
const double w1 = 0.0005;
const double E1 = 0.0001;
const double phi = M_PI/2.;
const double tc1 = M_PI/w1;
const double nc1 = 1.;

// for gaussian
const double tc = 10000; // 1000;
const double nc = 24;
const double fwhm = 744.12*2*2*2; // 744.12;
const double sigma = fwhm / (2. * sqrt (2. * log (2.)));
const double dt = (t1 - t0) / (nt - 1.);
const double A0 = E0 / w0;
const double Up = E0 * E0 / (4. * w0 * w0);

inline complex dipole (double p) {
  // return I * p;
  double dnom = (p*p + 2*Ip);
  return - I * 11.3137 * pow (2.*Ip, 1.25) * p
    / (M_PI * dnom * dnom * dnom); // pow(2., 3.5)
  // double alpha = 2.*Ip;
  // return I * pow (1./ (M_PI*alpha), 3./4.) * p / alpha * exp (- p*p / (2.*alpha));
}


}
