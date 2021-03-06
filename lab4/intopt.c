#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define M 20
#define N 20
#define eps 0.000001

typedef struct node_t node_t;
typedef struct n_list n_list;

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

struct node_t {
  node_t* next;
  int m;
  int n;
  int h;
  int k;
  double xh;
  double ak;
  double bk;
  double* min;
  double* max;
  double** a;
  double* b;
  double* x;
  double* c;
  double z;
};

struct n_list {
  n_list* next;
  n_list* prev;
  node_t* data;
};

void pivot(struct simplex_t * s, int row, int col);
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);
void succ(node_t* p, n_list** h, int m, int n,double** a, double* b, double* c, int k, double ak, double bk, double *zp, double* x);

node_t* initial_node(int m, int n, double **a, double *b, double* c){

  node_t *p = malloc(sizeof(node_t));
  int i, y;

  p->a = calloc(m+1, sizeof(double*));
  for (i = 0; i < m + 1; i++) {
      p->a[i] = calloc(n+1, sizeof(double));
  }
  p->b = calloc(m+1, sizeof(double));
  p->c = calloc(n+1, sizeof(double));
  p->x = calloc(n+1, sizeof(double));
  p->min = calloc(n, sizeof(double));
  p->max = calloc(n, sizeof(double));
  p->m = m;
  p->n = n;
  p->next = NULL;

  for(i=0; i<m; i++){
    for(y=0; y<n; y++){
      p->a[i][y] = a[i][y];
    }
  }

  for(i=0; i<m; i++){
    p->b[i] = b[i];
  }

  for(i=0; i<n; i++){
    p->c[i] = c[i];
  }

  for(i=0; i< n; i++){
    p->min[i] = -INFINITY;
    p->max[i]= INFINITY;
  }

  return p;
}

node_t* extend (node_t* p,int m, int n, double** a, double* b, double* c, int k, double ak, double bk) {
  node_t *q = malloc(sizeof(node_t));
  int i;
  int j;
  q->k = k;
  q->ak = ak;
  q->bk = bk;
  if(ak > 0 && p->max[k] < INFINITY){
    q->m = p->m;
  }else if(ak < 0 && p->min[k] > 0){
    q->m = p->m;
  }else{
    q->m = p->m+1;
  }
  q->n = p->n;
  q->h = -1;

  q->a = calloc(q->m+1, sizeof(double*));
  for(i=0; i<q->m+1;i++){
    q->a[i] = calloc(q->n+1, sizeof(double));
  }

  q->b = calloc(q->m + 1, sizeof(double));
  q->c = calloc(q->n + 1, sizeof(double));
  q->x = calloc(q->n + 1, sizeof(double));
  q->min = calloc(n, sizeof(double));
  q->max = calloc(n, sizeof(double));


  
  // copy over p->min, p->max to q!
  for(i=0; i<q->n; i++) {
    q->min[i] = p->min[i];
    q->max[i] = p->max[i];
  }

  for(i=0; i<n+1;i++){
    q->c[i] = c[i];
  }

  // copy m first rows of parameter a to q->a
  //
  for(i=0; i<m+1; i++) {
    for (j=0; j<n+1; j++) {
      q->a[i][j] = a[i][j];
    }
  }

  for (i=0; i < q->m; i++) {
    q->b[i] = b[i];
  }

  if(ak > 0) {
    if (q->max[k] == INFINITY || bk < q->max[k])
      q->max[k] = bk;
  } else if(q->max[k] == -INFINITY || -bk > q->min[k]) {
    q->min[k] = -bk;
  }

  for(i=m, j = 0; j < n; j++) {
    if (q->min[j] > -INFINITY) {
      q->a[i][j] = -1;
      q->b[i] = -q->min[j];
      i++;
    }

    if (q->max[j] < INFINITY) {
      q->a[i][j] = 1;
      q->b[i] = q->max[j];
      i++;
    }

  }

  return q;
}

int is_integer(double* xp) {
  double x = *xp;
  double r = round(x);
  if (fabs(r - x) < eps) {
    *xp = r;
    return 1;
  }
  return 0;
}

int integer(struct node_t* p) {
  int i;
  for(i = 0; i < p->n; i++) {
    if (!is_integer(&p->x[i])) {
      return 0;
    }
  }
  return 1;
}

void bound(struct node_t* p, n_list** h, double* zp, double* x) {
  if (p->z > *zp) {
    *zp = p->z;
    for(int i=0; i<p->n;i++) {
      x[i] = p->x[i];
    }

    n_list* current = *h;
    while (current != NULL) {
        n_list* next = current->next;
        if (current->data->z < p->z) {
            if (current->prev != NULL) {
                current->prev->next = current->next;
            } else {
                *h = current->next;
            }
            if (current->next != NULL) {
                current->next->prev = current->prev;
            }
            free(current->data->min);
            free(current->data->max);
            free(current->data);
            free(current);
        }
        current = next;
    }

  }
}

int branch(node_t *q, double z) {

  if (q->z < z) {
    return 0;
  }

  double min, max;
  int h;
  for (h=0; h < q->n; h++) {
    if (!is_integer(&q->x[h])) {
      if(q->min[h] == -INFINITY) {
        min = 0;
      } else {
        min = q->min[h];
      }
      max = q->max[h];
      if (floor(q->x[h]) < min || ceil(q->x[h]) > max) {
        continue;
      }
      q->h = h;
      q->xh = q->x[h];

      for(int i = 0; i < q->m+1; i++) {
        free(q->a[i]);
      }

      free(q->a);
      free(q->b);
      free(q->c);
      free(q->x);
      return 1;
    }
  }

  return 0;
}

void succ(node_t* p, n_list** h, int m, int n,double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x) {

  struct node_t* q = extend(p, m, n, a, b, c, k, ak, bk);
  if (q==NULL) {
    return;
  }

  q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
  if (isfinite(q->z)) {

    if (integer(q)) {
      bound(q, h, zp, x);

    } else if (branch(q, *zp)) {
        n_list* el = calloc(1, sizeof(n_list));
        el->data = q;
        el->next = *h;

        if (*h != NULL) {
            (*h)->prev = el;
        }

        *h = el;
        return;
    }
  }

  free(q->min);
  free(q->max);
  /* free(q->x); */
  /* free(q->b); */
  /* free(q->c); */
  /* for(int i=0; i < q->m+1; i++) { */
  /*   free(q->a[i]); */
  /* } */
  /* free(q->a); */
  free(q);
}


double intopt(int m, int n, double** a, double* b, double* c, double* x) {

  node_t* p = initial_node(m,n,a,b,c);
  n_list* h = calloc(1, sizeof(n_list));
  h->data = p;

  double z = -INFINITY;
  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0.0);  

  if (integer(p) || !isfinite(p->z)) {
    z = p->z;
    if (integer(p)) {
      for (int i=0; i<p->n; i++) {
        x[i] = p->x[i];
      }
    }

    for(int i=0; i<p->m + 1;i++){
      free(p->a[i]);
    }

    free(p->a);
    free(p->b);
    free(p->c);
    free(p->x);
    free(p->min);
    free(p->max);
    free(p);
    free(h);
    return z;
  }
  
  int br = branch(p, z);
  while(h != NULL) {
    p = h->data;
    n_list* old = h;
    if (h->next != NULL) {
        h->next->prev = NULL;
    }
    h = h->next;
    free(old);

    succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
    succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
    free(p->min);
    free(p->max);
    free(p);
  }

  if(z <= -INFINITY) {
    return NAN;
  }

  return z;
}



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

  int m = s->m;
  int n = s->n;
  int i;
  for(i = m+n; i>n;i--){
    s->var[i] = s->var[i-1];
  }
  s->var[n] =m+n;
  n++;
  for(i=0; i<m; i++){
    s->a[i][n-1] = -1;
  }
  s->x = calloc(m+n, sizeof(double));
  s->c = calloc(n+m, sizeof(double));
  s->c[n-1] = -1;
  s->n = n;
  pivot(s,k,n-1);
}

int initial(struct simplex_t * s, int m, int n, double ** a, double *b, double *c, double *x, double y, int *var){
  int i,j,k;
  double w;
  k = init(s,m,n,a,b,c,x,y,var);  
  if(b[k]>=0){
    return 1;
  }

  prepare(s,k);  
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);
  for(i=0; i<m+n; i++){
    if(s->var[i] == m+n-1){
      if(fabs(s->x[i]) > eps){
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
    for(j=0,k=0; k<n; k++){
      if(fabs(s->a[i-n][k]) > fabs(s->a[i-n][j])){
        j=k;
      }
      pivot(s,i-n,j);
      i=j;
    }
  }

  if(i<n-1){
    k= s->var[i];
    s->var[i] = s->var[n-1];
    s->var[n-1] =k;

    for(k=0; k<m; k++){      
      w = s->a[k][n-1];
      s->a[k][n-1] = s->a[k][i];
      s->a[k][i] =w;
    }
  }
  free(s->c);
  s->c=c;
  s->y=y;
  for(k=n-1;k<n+m-1;k++){
    s->var[k] = s->var[k+1];
  }
  n = s->n = s->n-1;
  double *t = calloc(n, sizeof(double));

  for(k=0;k<n;k++){
    for(j=0; j<n; j++){
      if(k == s->var[j]){
        t[j] = t[j] + s->c[k];
        goto next_k;
      }
    }
    for(j=0; j<m; j++){
      if(s->var[n+j]==k){
        break;
      }
    }
    s->y = s->y + s->c[k] * s->b[j];
    
    for(i=0; i<n; i++){
      t[i] = t[i] - s->c[k] * s->a[j][i];
    }
    next_k:;
  }

  for(i=0; i<n; i++){
    s->c[i] = t[i];
  }
  free(t);
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
 struct simplex_t s;
  int i, row, col;
  
  initial(&s, m, n, a, b, c, x, y, var);

  
  while(1) {
    col = select_nonbasic(&s);
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
      free(s.var);
      s.var = NULL;
      return INFINITY;
    }

    pivot(&s, row, col);
  }

  if (h==0) {
    for(i = 0; i < n; i++) {
      if (s.var[i] < n) {
        x[s.var[i]] = 0;
      }
    }

    for(i = 0; i < m; i++) {
      if (s.var[n+i] < n) {
        x[s.var[n+i]] = s.b[i];
      }
    }

   free(s.var);

  } else {
    for(i = 0; i < n; i++)
      x[i] = 0;

    for(i = n; i<n+m; i++)
      x[i] = s.b[i-n];
  }
   return s.y;

}


double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

