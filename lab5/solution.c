#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include<assert.h>
 
 
#define EPS 0.000001
#define UNBOUNDED 2147483647
 
typedef struct simplex_t simplex_t;
typedef struct node_t node_t;
typedef struct list_n list_n;
 
struct list_n {
    node_t* content;
    list_n* next;
    list_n* prev;
};
struct simplex_t{
    int n;
    int m;
    int* var;
    double** a;
    double* b;
    double* c;
    double* x;
    double y;
 
};
struct node_t { 
    node_t* next;
    int m; 
    int n;
    int k;
    int h;
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
    int id;
};
 
 
double** make_matrix(int m, int n); 
 
int is_integer(double* xp);
 
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);
 
int integer(node_t* p);
  
 
node_t* initial_node(int m, int n, double** a, double* b, double* c) {
    node_t* p = malloc(sizeof(node_t));
    int i, y;
 
    p->a = calloc(m + 1, sizeof(double*));
    for (i = 0; i < m + 1; i += 1) {
        p->a[i] = calloc(n + 1, sizeof(double));
    } 
 
    p->b = calloc(m+1, sizeof(double));
    p->c = calloc(n+1, sizeof(double));
    p->x = calloc(n+1, sizeof(double));
    p->min = calloc(n, sizeof(double));
    p->max = calloc(n, sizeof(double));
    p->m = m;
    p->n = n;
    p->next = NULL;
 
    for (i = 0; i < m; i += 1) {
        for (y = 0; y < n; y += 1) {
            p->a[i][y] = a[i][y];
        }
    }
 
    for (i = 0; i < m; i += 1) {
        p->b[i] = b[i];
    }
 
    for (i = 0; i < n; i += 1) {
        p->c[i] = c[i];
    }
 
 
    for (i = 0; i < n; i++){
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY; 
    }
 
 
    return p;
}
  
node_t* extend (node_t* p,int m, int n, double** a, double* b, double* c, int k, double ak, double bk) {
 
    assert(n>0);
    assert(k<n);
    assert(p != NULL);
 
    node_t* q;
    q = malloc(sizeof(node_t));
 
    int i,j;
 
    q->k = k;
    q->ak = ak;
    q->bk = bk;
 
    if(ak > 0 && p->max[k] < INFINITY){
        q->m = p->m;
    }else if(ak < 0 && p->min[k] > 0){
        q->m = p->m;
    }else{
        q->m = p->m + 1;
    }
 
    q->n = p->n;
    q->h = -1;
 
    q->a = calloc(q->m + 1, sizeof(double*));
    for (i = 0; i < q->m + 1; i += 1) {
        q->a[i] = calloc(q->n + 1, sizeof(double));
    }
 
    q->b = calloc(q->m+1, sizeof(double));
    q->c = calloc(q->n+1, sizeof(double));
    q->x = calloc(q->n+1, sizeof(double));
    q->min = calloc(n, sizeof(double));
    q->max = calloc(n, sizeof(double));
 
    memcpy(q->max, p->max, n * sizeof(double));
    memcpy(q->min, p->min, n * sizeof(double));
    for(i = 0; i < m+1; i++){
        memcpy(q->a[i],a[i], (n+1) * sizeof(double));
    }
 
    memcpy(q->b, b, m * sizeof(double));
    memcpy(q->c, c, (n+1) * sizeof(double));
 
    if(ak > 0){
        if(q->max[k] == INFINITY || bk < q->max[k]){
            q->max[k] = bk;    
        }
    }else if(q->max[k] == -INFINITY || -bk > q->min[k]){
            q->min[k] = -bk;
    }
 
    for(i = m, j = 0; j < n; j++){
        if(q->min[j] > -INFINITY){
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i += 1;        
        }
        if(q->max[j] < INFINITY){
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i +=1;        
        }
    } 
 
    return q;
}
 
 
void bound(node_t* p,list_n** h, double* zp, double* x){
    if(p->z > *zp){
        *zp = p->z; //TODD: maybe select best x
        for (int i = 0; i < p->n; i += 1) {
            x[i] = p->x[i];
        }
 
 
        list_n* current = *h;
        while (current != NULL) {
            list_n* next = current->next;
            if (current->content->z < p->z) {
                if (current->prev != NULL) {
                    current->prev->next = current->next;
                } else {
                    *h = current->next;
                }
                if (current->next != NULL) {
                    current->next->prev = current->prev;
                }
                free(current->content->min);
                free(current->content->max);
                free(current->content);
                free(current);
            }
            current = next;
        }
    }
}
 
 
 
 
int branch(node_t* q, double z){
    double min, max;
  
    if(q->z < z){
        return 0;
    }
 
    int h = 0;
    do {
        if (!is_integer(&q->x[h])) {
            if (q->min[h] == -INFINITY) {
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
            for (int i = 0; i < q->m + 1; i += 1) {
                free(q->a[i]);
            }
            free(q->a);
            free(q->b);
            free(q->c);
            free(q->x);
            return 1;
        }
        h++;
    } while (h < q->n);
    return 0;
}
 
void succ(node_t* p, list_n** h, int m, int n, double** a, double* b, double* c,int k, double ak, double bk, double* zp, double* x){
    node_t* q = extend(p,m,n,a,b,c,k,ak,bk);
    if (q == NULL){
        return;
    }
    q->z = simplex(q->m,q->n,q->a,q->b,q->c,q->x,0);
 
    if(isfinite(q->z)){
        if(integer(q)){
            bound(q,h,zp,x);
        }else if (branch(q,*zp)){
            list_n* new_list_element = calloc(1, sizeof(list_n));
            new_list_element->content = q;
            new_list_element->next = *h;
            if (*h != NULL) {
                (*h)->prev = new_list_element;
            }
            *h = new_list_element;
            return;
        }
    }
    for (int i = 0; i < q->m + 1; i += 1) {
        free(q->a[i]);
    }
    free(q->a);
    free(q->b);
    free(q->c);
    free(q->x);
    free(q->min);
    free(q->max);
    free(q);
}
 
int init(simplex_t* s, int m, int n,double** a, double* b, double* c,double* x, double y, int* var){
 
    int i,k;
 
    s->n = n;
    s->m = m;
    s->a = a;
    s->b = b;
    s->c = c;
    s->var = var;
    s->x = x;
    s->y = y;
 
    if (s->var == NULL){
        s->var = calloc(n + m + 1,sizeof(int));
        for (int i = 0; i < n + m; i++){
            s->var[i] = i;
        }
    }
    for (k = 0, i = 1; i < m; i++){
        if (b[i] < b[k]){
            k = i;
        }
    }
 
 
    return k;
}
 
void pivot(simplex_t* s, int row, int col)
{
    double** a = s->a;
    double* b = s->b;
    double* c = s->c;
    int m = s->m;
    int n = s->n;
    int i,j,t;
 
    t = s->var[col];
    s->var[col] = s->var[n + row];
    s->var[n + row] = t;
    double f = 1/a[row][col]; 
    s->y = s->y + c[col] * b[row]*f;
    for (int i = 0; i < n; i++)
    {
        if (i != col){
            c[i] = c[i] - c[col] * a[row][i]*f;
        }
    }
 
    c[col] = -c[col]*f;
 
    for (int i = 0; i < m; i++){
        if(i != row){
            b[i] = b[i] - a[i][col]*b[row]*f;
        }
    }
 
    for (i = 0; i < m; i++){
        if (i != row){
            for(j= 0; j<n; j++){
                if(j != col){
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] *f;
                }
            }
        }
    }
 
    for (i = 0; i < m; i++){
        if (i != row){
            a[i][col] = -a[i][col] *f;
        }
    }
 
    for (i = 0; i < n; i++){
        if (i != col){
            a[row][i] = a[row][i] *f;
        }   
     }
    b[row] = b[row] *f;
    a[row][col] = f;
}
 
void prepare(simplex_t* s,int k){
    //printf("I am preparing\n");
 
    int m = s->m;
    int n = s->n;
    int i;
 
    for (i =m+n; i>n;i--){
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m+n;
    n++;
 
    for (i = 0; i < m; i++){
        s->a[i][n-1] = -1;  
    }
 
    s->x = calloc(m+n, sizeof(double));
    s->c = calloc(m+n, sizeof(double));
 
    s->c[n-1] = -1;
    s->n = n;
    pivot(s,k,n-1);
}
 
double select_nonbasic(simplex_t* s)
{
    for(int i=0; i < s->n;i++)
    {
        if(s->c[i] > 0.0000001)
        {
            return i;
        }
    }
    return -1;
}
 
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);
 
 
int initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var){
 
    int i,j,k; 
    double w;
 
    k = init(s, m, n, a, b, c, x, y, var);
 
 
    if(b[k] >= 0){
        return 1;
    }
    // Rest of func inplemented in lab3
    prepare(s,k);
 
 
    n = s->n;
    s->y = xsimplex(m,n,s->a,s->b,s->c,s->x,0,s->var,1);
 
    for (i = 0; i < m+n; i++){
        if(s->var[i] == m + n - 1){
            if(fabs(s->x[i]) > EPS){
                // TODO: Is This right?
                free(s->x);
                free(s->c);
                return 0;
            }else{
                break;
            }
        }
    }
 
    if (i >= n){
        // xn+m is basic. find good nonbasic.
        for(j = 0, k = 0; k < n; k++){
            if(fabs(s->a[i-n][k]) > fabs(s->a[i-n][j])){
                j = k;
            }
        }
        pivot(s,i-n,j);
        i =j;
    }
 
    if(i < n-1){
        // xn+m is nonbasic and not last. swap columns i and n-1
        k = s->var[i];
        s->var[i] = s->var[n-1];
        s->var[n-1] = k;
        for (k = 0; k<m; k++){
            w = s->a[k][n-1];
            s->a[k][n-1] = s->a[k][i];
            s->a[k][i] = w;
        }
    }
        // xn+m is nonbasic and last. forget it.
        // TODO: This is baced on that the sevdo code is wrong. I put everything below in the else statment
        // but not sure 
 
        free(s->c);
        s->c = c;
        s->y = y;
 
        for (k = n-1; k < n+m-1; k++){
            s->var[k] = s->var[k+1];
        }
        n = s->n = s->n - 1;
 
 
        double* t;
        t = calloc(n, sizeof(double));
 
 
        for (k = 0; k < n; k++){
            for (j = 0; j < n; j++){
                if(k == s->var[j]){
                    t[j] = t[j] + s->c[k];
                    goto next_k; 
                }
            }
 
            for (j = 0; j < m; j++){
                if(s->var[n+j] == k){
                    // TODO: THIS BREAK MAKES NO sesne
                    break;
                }
            }
            s->y = s->y +s->c[k]*s->b[j];
 
            for (i = 0; i < n; i++){
                t[i] = t[i] - s->c[k]*s->a[j][i];
            }
        next_k:;
        }
        for (i = 0; i < n; i++){
            s->c[i] = t[i];
        }
        // TODO: are we going to free all variables that is writen as DELETE? YES
        free(t);
        free(s->x);
        return 1;
}
 
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h){
 
    simplex_t s; 
 
    int i,row,col;
 
 
    if(!initial(&s, m, n, a, b, c, x, y, var)){
        free(s.var);
        s.var = NULL;
        return NAN;
    }
 
    while ((col = select_nonbasic(&s)) >= 0 )
    {
        row =  -1; 
        for (i = 0; i < m; i++){
            if((a[i][col] > EPS) && ((row < 0) || ((b[i] / (a[i][col]) < (b[row]/ a[row][col] ))))) {
                row = i;
            }
        }
        if (row < 0){
            free(s.var);
            s.var = NULL;
            return INFINITY;
        }
        pivot(&s, row, col); 
 
    }
 
    if(h == 0){
        for(i = 0; i < n;i++){
            if(s.var[i] < n){
                x[s.var[i]] = 0;
            }
        }
        for(i = 0; i < m;i++){
            if(s.var[n+i] < n){
                x[s.var[n+ i]] = s.b[i];
            }
        }
        free(s.var);
        s.var = NULL;
    }else{
        for (i = 0; i < n; i++){
            x[i] = 0;
        }
        for (i = n; i < n+m; i++){
            x[i] = s.b[i-n];
        }
    }
    return s.y;
}
 
 
 
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y){
    return xsimplex(m,n,a,b,c,x,y,NULL,0);
}
 
 
double** make_matrix(int m, int n)
{
    double**       a;
    int            i;
    a = calloc(m, sizeof(double*));
    for (i = 0; i < m; i += 1){
            a[i] = calloc(n, sizeof(double)); 
}
    return a;
}
 
int is_integer(double* xp){
    double x = *xp;
    double r = round(x);
    if(fabs(x-r) < EPS){
        *xp = r;
        return 1;
    }else{
        return 0;
    }
}
 
int integer(node_t* p){
    int i;
    for(i = 0; i < p->n; i++){
        if(!is_integer(&p->x[i])){
            return 0;
        }
    }
    return 1;
}
 
double intopt(int m, int n, double** a, double* b, double* c, double* x){
 
    node_t* p = initial_node(m,n,a,b,c);
 
    // TODO: maybe change as in the working code
    list_n* h = calloc(1, sizeof(list_n));
    h->content = p;
 
    double z = -INFINITY;
    p->z = simplex(p->m,p->n,p->a,p->b,p->c,p->x,0);
 
    if(integer(p) || !isfinite(p->z)){
        z = p->z;
        if(integer(p)){
            for (int i = 0; i < p->n; i += 1) {
                x[i] = p->x[i];
            }
        }
 
 
 
        for (int i = 0; i < p->m + 1; i += 1) {
            free(p->a[i]);
        }
        free(p->a);
        free(p->b);
        free(p->c);
        free(p->x);
        free(p->min);
        free(p->max);
        free(p);
        free(h); // TODO: Not precent in "working code"
        return z;     // Same comment as above
    }
    branch(p,z);
 
    while(h != NULL){
       // printf("inside while Head: %d\n",h->size);
 
        // Pop
        p = h->content;
        list_n* old = h;
        if (h->next != NULL) {
            h->next->prev = NULL;
        }
        h = h->next;
        free(old);
 
        succ(p,&h,m,n,a,b,c,p->h,1,floor(p->xh),&z,x);
        succ(p,&h,m,n,a,b,c,p->h,0,-ceil(p->xh),&z,x);
 
        // Maybe free more
        free(p->min);
        free(p->max);
        free(p);
 
    }   
    //free(h);
    if(z == -INFINITY){
        return NAN;
    }else{
        return z;
    }
}


int main(int argc, char** argv){
    int     m;
    int     n;
    scanf("%d %d\n", &m, &n);
    printf("m = %d, n = %d\n",m,n);
    double** matrix;
    
    matrix = calloc(m+1, sizeof(double*));
    for (int i = 0; i < m +1; i++){
        matrix[i] = calloc(n+1, sizeof(double));
    }
    
    
    double* c_vector = calloc(n+1, sizeof(double));
    double* b_vector = calloc(m+1, sizeof(double));
    double* x_vector = calloc(n, sizeof(double));
    double y;
    int i;
    int j; 
    
    for(i = 0; i < n; i++){
        scanf("%lf", &c_vector[i]);
    }
    
    
    for(i = 0;i < m; i++){
        for(j = 0; j < n; j++){
            scanf("%lf", &matrix[i][j]);
        }
    }
    
    
    for(j = 0; j < m; j++){
        scanf("%lf", &b_vector[j]);
    }
    
    
    printf("max z =");
    
    for(j = 0; j < n; j++){
        printf("%10.3lf x%d ", c_vector[j],j);
    }
    
    printf("\n");
    
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            printf("%10.3lf  x%d", matrix[i][j], j); 
        }
        printf(" <= %10.3lf", b_vector[i]);
        printf("\n");
    }
    
    
    
    int simp = intopt(m, n, matrix, b_vector, c_vector, x_vector);
    
    printf("%d \n", simp); 
    
    free(c_vector);
    free(b_vector);
    free(x_vector);
    for(i = 0; i < m+1; i++){
        free(matrix[i]);
    }
    free(matrix);
} 
 
