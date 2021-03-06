#include <stdio.h>
#include <stdlib.h>

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

int main()
{
  int m, n;
  scanf("%d%d", &m, &n);

  double* c_vector = calloc(n, sizeof(double));    
  double ** matrix;
  matrix = calloc(m, sizeof(double*));  
  for (int j = 0; j<m; j++) {
    matrix[j] = calloc(n, sizeof(double));
  }
  //free(j);

  double* b_vector = calloc(m, sizeof(double));

  for(int i=0; i<n; i+=1){
    scanf("%lf", &c_vector[i]);
  }

  for(int row=0; row<m; row+=1){
    for(int col=0; col<n; col+=1){
      scanf("%lf", &matrix[row][col]);
    }
  }

  for(int i=0; i<m; i+=1){
    scanf("%lf", &b_vector[i]);
  }

  printf("%s", "max z = ");
  for(int i=0; i<n; i+=1){
      printf("%10.3lf ", c_vector[i]);
    if(i != n-1){
      printf("%s%d + ", "x", i);
    }else{
      printf("%s%d ", "x", i);
    }
  }
  free(c_vector);
  printf("\n");
  printf("\n");

  for(int row=0; row<m; row+=1){
    for(int col=0; col<n; col+=1){
      if(col != n-1){
        printf("%10.3lfx%d + ", matrix[row][col], col);
      }else{
        printf("%10.3lfx%d ", matrix[row][col], col);
      }
    }
    printf("\u2264 %10.3lf", b_vector[row]);
    printf("\n");
  }
  for(int i=0; i<m; i++){
    free(matrix[i]);
  }
  free(matrix);
  free(b_vector);  

  return 0;
}


