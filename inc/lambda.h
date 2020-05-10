#ifndef _LAMBDA_H_
#define _LAMBDA_H_

#ifdef __cplusplus
extern "C" {
#endif

int lambda(int n, int m, const double *a, const double *Q, double *F, double *s);
int lambda_search(int n, int m, const double *a, const double *Q, double *F, double *s);

#ifdef __cplusplus
}
#endif

#endif