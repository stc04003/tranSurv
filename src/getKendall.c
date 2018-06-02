#include <R.h>
#include <Rmath.h>
#include <math.h>

void uCondKendall(double *x, double *y, int *n, double *out) {
  int i, j;  
  double M = 0.0;
  double tmp = 0.0;
  double tau = 0.0;
  for (int i = 0; i < (*n - 1); i++) {
    for (int j = (i + 1); j < *n; j++) {
      tmp = (y[i] - y[j]) * (x[i] - x[j]);
      tau += (tmp > 0) - (tmp < 0);
      M += 1;
    }
  }
  out[0] = tau / M;
  out;
}

void condKendall(double *t1, double *t2, double *d, int *n, double *weights, 
		    int *meth, double *out) {
  int i, j, k;
  double *bb = Calloc(*n * (*n - 1), double);
  double Uc; 
  double Um; 
  double tmp;
  double wgt;
  double sgn;
  double gamma;
  double v = 0.0;
  double v1 = 0.0;
  double v2 = 0.0;
  double M = 0.0;
  for (i = 0; i < (*n - 1); i++) {
    for (j = i + 1; j < *n; j++) {
      if (fmax(t1[i], t1[j]) <= fmin(t2[i], t2[j]) &&
	  d[i] * (t2[i] <= t2[j]) + d[j] * (t2[j] <= t2[i]) > 0 &&
	  weights[i] * weights[j] * weights[i + *n] * weights[j + *n] > 0) {
	tmp = (t1[i] - t1[j]) * (t2[i] - t2[j]);
	sgn = (tmp > 0) - (tmp < 0);
	if (*meth == 1) wgt = 1;
	if (*meth == 2) wgt = fmax(weights[j + *n], weights[i + *n]) *
			  fmax(weights[j + *n], weights[i + *n]) / (weights[j] * weights[i]);
	if (*meth == 3) {
	  wgt = weights[j + *n] * weights[i + *n] / (weights[j] * weights[i]);
	  Uc += d[i] * d[j] * sgn / wgt;
	  Um += d[i] * d[j] / wgt;
	  bb[i * (*n - 1) + j - 1] = bb[j * (*n - 1) + i] = d[i] * d[j] * sgn / wgt;
	} else {
	  Uc += sgn / wgt;
	  Um += 1 / wgt;
	  bb[i * (*n - 1) + j - 1] = bb[j * (*n - 1) + i] = sgn / wgt;
	}
	wgt = 0.0;
      }
    }
  }
  out[0] = Uc / Um;
  for (i = 0; i < *n; i++) {
    for (j = 0; j < (*n - 1); j++) {
      v1 += bb[i * (*n - 1) + j];
      v2 += bb[i * (*n - 1) + j] * bb[i * (*n - 1) + j];
    }
    v += (v1 * v1 - v2) / *n;
    v1 = 0.0; 
    v2 = 0.0;
  }
  out[1] = v * *n * (*n - 1) / (Um * Um * (*n - 2));
  Free(bb);
  out;
}

void pmcc(double *t1, double *t2, int *n, double *out) {
  int i, j, k;
  double STT = 0.0;
  double SXT = 0.0;
  double SXX = 0.0;
  double tmp = 0.0; 
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      if (fmax(t1[i], t1[j]) <= fmin(t2[i], t2[j])) {
	STT += (t1[i] - t1[j]) * (t1[i] - t1[j]);
	SXT += (t1[i] - t1[j]) * (t2[i] - t2[j]);
	SXX += (t2[i] - t2[j]) * (t2[i] - t2[j]);
      }
    }
  }
  out[0] = SXT / sqrt(STT * SXX); 
  for (i = 0; i < *n; i++) {
    tmp = 0.0;
    for (j = 0; j < (*n - 1); j++) {
      if (fmax(t1[i], t1[j]) <= fmin(t2[i], t2[j])) {
	tmp += (t2[i] - t2[j]) * (t2[i] - t2[j]) / SXX + (t1[i] - t1[j]) * (t1[i] - t1[j]) / STT -
	  2 * (t1[i] - t1[j]) * (t2[i] - t2[j]) / SXT;
      }
    }
    out[1] += tmp * tmp;     
  }
  out[1] = out[1] * out[0] * out[0]; 
  out;
}

