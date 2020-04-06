#pragma once

#include <cmath>
#include <complex>


#define complex std::complex<double>

class Pulse_mono {
  double w;   // frequency
  double E0;  // field amplitude
  double A0;
 public:
  Pulse_mono (double w, double E0) :
  w (w), E0 (E0) { A0 = E0 / w; }
  inline double A (double t) {
    return A0 * cos (w * t);
  }
  inline double E (double t) {
    return E0 * sin (w * t);
  }
};

class Pulse_sin2 {
  double w;   // frequency
  double E0;  // field amplitude
  double tc;  // center of pulse
  double nc;  // number of cycles
  double A0;
  double duration;
 public:
  Pulse_sin2 (double w, double E0, double tc, double nc) :
  w (w), E0 (E0), tc (tc), nc (nc) {
    A0 = E0 / w;
    duration = 2. * M_PI * nc / w;
  }
  double A (double t) {
    double result = 0.;
    auto phi = w * (t - tc + duration / 2.);
    if (phi >= 0 && phi < w * duration) {
      auto e = sin (phi / (2.*nc));
      result = A0 * e * e * sin (phi);
    }
    return result;
  }
  double E (double t) {
    double result = 0.;
    auto phi = w * (t - tc + duration / 2.);
    if (phi >= 0 && phi < w * duration) {
      auto phi_m = (1. - 1./nc) * phi;
      auto phi_p = (1. + 1./nc) * phi;
      result = - A0 / (4.*(nc*nc-1)*w)
        * (2.*(1.-nc*nc) * cos (phi)
           + nc * ((nc + 1.) * cos (phi_m) + (nc - 1.) * cos (phi_p)));
    }
    return result;
  }
};

class Pulse_gaussian {
  double w;   // frequency
  double E0;  // field amplitude
  double tc;  // center of pulse
  double fwhm;// FWHM
  double A0;
  double sigma;
 public:
  Pulse_gaussian (double w, double E0, double tc, double fwhm) :
    w (w), E0 (E0), tc (tc), fwhm (fwhm) {
    A0 = E0 / w;
    sigma = fwhm / (2. * sqrt (2. * log (2.)));
  }
  inline double A (double t) {
    auto tt = t - tc;
    double res = A0 * exp (- tt * tt / (2. * sigma * sigma)) * sin (w * tt);
    if (std::isnan (res)) res = 0.;
    return res;
  }
  inline double E (double t) {
    auto tt = t - tc;
    // double res = A0 * exp (-0.5 * sigma * sigma * w * w) * sqrt (0.5*M_PI) * sigma
    //   * imag (cerf (M_SQRT1_2 * complex (tt/sigma, sigma * w)));
    double sigma2 = sigma * sigma;
    double res = A0 / sigma2 * exp (- tt * tt / (2. * sigma2))
      * (tt * sin (w * tt) - w * sigma2 * cos (w * tt));
    if (std::isnan (res)) res = 0.;
    return res;

  }
};
