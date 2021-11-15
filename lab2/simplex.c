#include <stdio.h>
#include <stdlib.h>

#define M 20
#define N 20

struct simplex_t{
  int m;
  int n;
  double** a;
  double *b;
  double *c;
  double *x;
  double y;
  int *var;

  /* double a[M][N+1]; */
  /* double b[M]; */
  /* double c[N]; */
  /* double x[N+1]; */
  /* double y; */
  /* int var[N+M+1]; */
};

int init(struct simplex_t * s, int m, int n, double ** a, double *b, double *c, double *x, double y, int *var){
  int i,k;
  s->n = n;
  s->m = m;
  s->a = a;
  s->b = b;
  s->c = c;
  s->x = x;
  s->y = y;
  s->var = var;
  if(s->var == NULL){
    s->var = calloc(m+n+1, sizeof(int));
    for(i=0;i<n+m+1;i++){
      s->var[i] = i;
    }
  }
  for(k =0, i=1; i<m; i++){
    if(b[i]<b[k]){
      k=i;
    }
  }
  return k;
}

int main(){
  struct simplex_t *s  = malloc(sizeof (struct simplex_t));
  return 0;
}
