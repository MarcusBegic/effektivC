#include <stdio.h>
#include <stdlib.h>

#define M 20
#define N 20
#define eps 10e-6
#define INF 10e12

double** make_matrix(int m, int n)
{
  double **     a;
  int           i;

  a = calloc(m, sizeof(double*));
  for(i=0; i<m; i+=1){
    a[i] = calloc(n, sizeof(double));
  }
  return a;
}

struct simplex_t{
  int m;
  int n;
  double** a;
  double *b;
  double *c;
  double *x;
  double y;
  int *var;
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
    for(i=0;i<n+m;i++){
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

int select_nonbasic(struct simplex_t *s){
  int i;
  for(i=0; i<s->n; i++){
    if(s->c[i] > eps){
      return i;
    }
  }
  return -1;
}

void prepare(struct simplex_t * s, int k){
  printf("%s", "reaches prepare");
  int m = s->m;
  int n = s->n;
  int i;
  for(i = m+n; i<n;i++){
    s->var[i] = s->var[i-1];
  }
  s->var[n] =m+n;
  n++;
  for(i=0; i<m; i++){
    s->a[i][n-1] = -1;
  }
  s->x = calloc(m+n, sizeof(double));
  s->c = calloc(n, sizeof(double));
  s->c[n-1] = -1;
  s->n = n;
}

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h); 

void pivot(struct simplex_t * s, int row, int col);


int initial(struct simplex_t * s, int m, int n, double ** a, double *b, double *c, double *x, double y, int *var){
  int i,j,k;
  double w;
  k = init(s,m,n,a,b,c,x,y,var);
  if(b[k]>=0) return 1;
  prepare(s,k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);
  for(i=0; i<m+n; i++){
    if(s->var[i] == m+n-1){
      if(s->x[i] > eps){
        free(s->x);
        free(s->c);
        return 0;
      }
      else{
        break;
      }
    }
  }
  if(i >= n){
    for(j=k=0; k<n; k++){
      if(abs(s->a[i-n][k]) > abs(s->a[i-n][j])){
        j=k;
      }
      pivot(s,i-n,j);
    }
  }
  if(i<n-1){
    k= s->var[i];  s->var[i] = s->var[n-1]; s->var[n-1] =k;
    for(k=0; k<m; k++){
      w = s->a[k][n-1];
      s->a[k][n-1] = s->a[k][i];
      s->a[k][i] =w;
    }
  }
  free(s->c);
  s->c=c;
  s->y=y;
  for(k=n-1;k<n+m+1;k++){
    s->var[k] = s->var[k+1];
  }
  n=s->n=s->n-1;
  double t[100]={0.0};

  for(k=0;k<n;k++){
    for(j=0; j<n; j++){
      if(k == s->var[j]){
        t[j] = t[j]+s->c[k];
        goto next_k;
      }
    }
    for(j=0; j<m; j++){
      if(s->var[n+j]==k){
        break;
      }
    }
    s->y = s->y+s->c[k]*s->b[j];
    for(i=0; i<n; i++){
      t[i] = t[i] - s->c[k] * s->a[j][i];
    }
    next_k:
  }
  for(i=0; i<n; i++){
    s->c[i] = t[i];
  }
  /* free(t); */
  free(s->x);
  return 1;
}

int glob=0;
void pivot(struct simplex_t * s, int row, int col) {
  double** a = s->a;
  double* b = s->b;
  double* c = s->c;

  int m = s->m;
  int n = s->n;
  int i,j,temp;

  temp = s->var[col];
  s->var[col] = s->var[n+row];
  s->var[n+row] = temp;
  s->y = s->y + c[col] * b[row] / a[row][col];
  for (i = 0; i<n; i++) {
    if(i != col)
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];
  }

  c[col] = -c[col] / a[row][col];
  for (i = 0; i < m; i++) {
    if (i != row) 
      b[i] = b[i] - a[i][col] * b[row] / a[row][col];
  }

  for (i = 0; i < m; i++) {
    if(i != row) {
      for(j = 0; j < n; j++) {
        if (j != col) {
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        }
      }
    }
  }
  glob+=1;

  for(i=0; i<m; i++){
    if(i!=row){
      a[i][col] = -a[i][col] / a[row][col];
    }
  }

  for(i=0; i<n; i++) {
    if (i != col) {
      a[row][i] = a[row][i] / a[row][col];
    }
  }
  glob+=1;
  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
}


double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
  struct simplex_t *s  = malloc(sizeof (struct simplex_t));
  int i, row, col;

  initial(s, m, n, a, b, c, x, y, var);

  while(1) {
    col = select_nonbasic(s);
    if (col < 0) {
      break;
    }
    row =-1;

    for (i = 0; i < m; i++) {
      if (a[i][col] > eps && (row < 0 || b[i]/a[i][col] < b[row] / a[row][col])) {
        row = i;
      }
    }

    if (row < 0) {
      return INF;
    }

    pivot(s, row, col);    
  }
  for(i=0; i<m; i++){
    free(s->a[i]);
  }
  free(s->a);
  free(s->c);

  if (h==0) {
    for(i = 0; i < n; i++) {
      if (s->var[i] < n) {
        x[s->var[i]] = 0;
      }
    }

    for(i = 0; i < m; i++) {
      if (s->var[n+i] < n) {
        x[s->var[n+i]] = s->b[i];
      }
    }

   free(s->var);

  } else {
    for(i = 0; i < n; i++)
      x[i] = 0;

    for(i = n; i<n+m; i++)
      x[i] = s->b[i-n];  
  }
  free(s->b);
  free(s->x);
  int temp = s->y;
  free(s);

  return temp;
  
}


double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}


int main(){


  int m, n;
  scanf("%d%d", &m, &n);

  double* c = calloc(n, sizeof(double));    
  for(int i=0; i<n; i+=1){
    scanf("%lf", &c[i]);
  }

  double ** a;

  a = calloc(m, sizeof(double*));  
  for (int j = 0; j<m; j++) {
    a[j] = calloc(n, sizeof(double));
  }
  for(int row=0; row<m; row+=1){
    for(int col=0; col<n; col+=1){
      scanf("%lf", &a[row][col]);
    }
  }

  double* b = calloc(m, sizeof(double));

  for(int i=0; i<m; i+=1){
    scanf("%lf", &b[i]);
  }

  double * x = calloc(n+1, sizeof(double));

  double y = simplex(m, n, a, b, c, x, 0.0);
  printf("%lf", y);

  return 0;
}
