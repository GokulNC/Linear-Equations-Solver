/*******************************************************
*                                                      *
* File Description  : Standalone utility module for    *
*                     managing matrix operations       *
*                                                      *
* Written by        : Gokul N.C.                       *
*                     ( gokulnc@ymail.com )            *
*                     http://about.me/GokulNC          *
*                                                      *
*******************************************************/

#include <stdlib.h>

float** getNewMatrix(int rows, int columns) {
    float **mat = (float**) malloc(sizeof(float*)*rows);
    for(int i=0; i<rows; ++i)
        mat[i] = (float*) malloc(sizeof(float)*columns);
    return mat;
}

void getCofactorMatrix(float **mat, float **cofactorMatrix, int dim, int row, int column) {
    int fill_i=0, fill_j;
    for(int i=0; i<dim; ++i) {
        if(i==row) continue;
        for(int j=0, fill_j=0; j<dim; ++j) {
            if(j==column) continue;
            cofactorMatrix[fill_i][fill_j++] = mat[i][j];
        }
        fill_i++;
    }
}

float determinant(float **mat, int n) {
    if(n==1) return mat[0][0];
    float result = 0.0;
    char sign = 1;
    float **cofactorMatrix = getNewMatrix(n-1, n-1);
    for(int j=0; j<n; ++j) {
        getCofactorMatrix(mat, cofactorMatrix, n, 0, j);
        result += sign*mat[0][j]*determinant(cofactorMatrix, n-1);
        sign = -sign;
    }
    free(cofactorMatrix);
    return result;
}

void getAdjointMatrix(float **mat, float **adj, int n) {
    char sign;
    float **cofactorMatrix = getNewMatrix(n-1, n-1);
    for(int i=0; i<n; ++i) {
        sign = i&1 ? -1 : +1;
        for(int j=0; j<n; ++j) {
            getCofactorMatrix(mat, cofactorMatrix, n, i, j);
            adj[j][i] = sign * determinant(cofactorMatrix, n-1);
            sign = -sign;
        }
    }
    free(cofactorMatrix);
}

void elementWiseMultiply(float val, float **mat, int r, int c) {
    for(int i=0; i<r; ++i)
        for(int j=0; j<c; ++j)
            mat[i][j] *= val;
}

// Returns -1 if matrix cannot be inverted.
int getInverseMatrix(float **mat, float **inverse, int n) {
    getAdjointMatrix(mat, inverse, n);
    float det = determinant(mat, n);
    if(det==0.0) {
        // No inverse exists for given matrix.
        return -1;
    }
    elementWiseMultiply(1.0/det, inverse, n, n);
    return 0;
}

float** multiplyMatrices(float **A, float **B, int r1, int c1, int r2, int c2) {
    float sum;
    if(c1 != r2) {
        printf("ERROR: Cannot multiply the given matrices!\n");
        exit(-1);
    }
    float **product = getNewMatrix(r1, c2);
    for(int i=0; i<r1; ++i)
        for(int j=0; j<c2; ++j) {
            sum = 0.0;
            for(int k=0; k<c1; ++k)
                sum += A[i][k]*B[k][j];
            product[i][j] = sum;
        }
    return product;
}

void printMatrix(float **mat, int rows, int columns) {
    for(int i=0; i<rows; ++i) {
        for(int j=0; j<columns; ++j)
            printf("%0.2f\t", mat[i][j]);
        printf("\n");
    }
}

void printVector(float *vec, int n) {
    for(int i=0; i<n; ++i)
        printf("%0.3f\n", vec[i]);
    printf("\n");
}
