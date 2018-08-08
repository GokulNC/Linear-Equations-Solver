/*******************************************************
*                                                      *
* File Description  : Program for solving a set of     *
*                      multivariate linear equations   *
*                                                      *
* Written by        : Gokul N.C.                       *
*                     ( gokulnc@ymail.com )            *
*                     http://about.me/GokulNC          *
*                                                      *
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrixOps.c"
#define MAX_EQNS 50

typedef struct {
    char alpha;
    float val;
} coeff;

typedef struct {
    coeff *vars;
    float constant;
    int numVars;
} equation;

typedef struct {
    equation *eqns;
    int numEquations;
} familyOfEqns;

int countVariables(char *line) {
    int count = 0;
    for(int i=0; line[i]!='\0'; ++i)
        if(line[i]>='a' && line[i]<='z') ++count;
    return count;
}

//parses number from str[i] into val and return no. of chars processed
int parseNumber(char *str, int i, float *val, float defaultVal) {
    int sign = 1, floatLocation = -1, parseLen = 0;
    *val = 0;
    while(1) {
        if(str[i] >= '0' && str[i] <='9') {
            if(floatLocation>0) {
                *val += (str[i]-'a')/((float) pow(10.0, floatLocation));
                floatLocation += 1;
            } else {
                floatLocation = 0;
                *val = (*val)*10.0 + str[i]-'0';
            }
        }
        else if(str[i]=='+') sign = 1;
        else if(str[i]=='-') sign = -1;
        else if(str[i]=='.') floatLocation = 1;
        else {
            if(floatLocation==-1) {
                //if no number is given, assume defaultVal
                *val = sign*defaultVal;
            } else {
                *val = sign*(*val);
            }
            return parseLen;
        }
        ++parseLen; ++i;
    }
}

// parses and stores coefficient & corresponding variable
// name at str[i] into 'var' (a struct container for both)
int parseToken(char *str, int i, coeff *var) {
    if(str[i]=='=') return 0;
    int parseLen = parseNumber(str, i, &(var->val), 1.0);
    if(parseLen==-1) {
        printf("\nProblem parsing the equation at index %d\n", i);
    }
    i += parseLen;
    if(str[i]>='a' && str[i]<='z') {
        var->alpha = str[i];
        return parseLen+1;
    } else {
        return -1;
    }
}

void parseEquation(char *line, equation *eqn) {
    eqn->numVars = countVariables(line);
    eqn->vars = (coeff*) malloc(sizeof(coeff)*eqn->numVars);
    int index = 0, parseLen, varIndex = 0;
    while(line[index] != '\0') {
        parseLen = parseToken(line, index, eqn->vars+varIndex);
        if(parseLen>0) {
            index += parseLen;
            ++varIndex;
        } else if(parseLen==0) {
            index += 1;
            parseNumber(line, index, &(eqn->constant), 0.0);
            break;
        } else {
            printf("\nProblem parsing this line: %s\n", line);
            exit(-1);
        }
    }

    if(eqn->numVars != varIndex) {
        printf("\nError parsing variable coeffs\n");
        exit(-1);
    }
}

int verifyEquations(familyOfEqns *eqns) {
    if(eqns->numEquations==0) return -1;
    equation targetEqn = eqns->eqns[0];
    if(eqns->numEquations != targetEqn.numVars) return -2;
    for(int i=1; i< (eqns->numEquations); ++i) {
        if(eqns->eqns[i].numVars != targetEqn.numVars) {
            printf("Equation %d does not match with the previous equations.\n", i);
            return -3;
        }
        for(int j=0; j< targetEqn.numVars; ++j) {
            if(eqns->eqns[i].vars[j].alpha != targetEqn.vars[j].alpha) {
                printf("Equation %d does not have variable %c in proper sequence.\n", i, targetEqn.vars[j].alpha);
                return -4;
            }
	}
    }
    return 0;
}

familyOfEqns readEquations() {
    familyOfEqns eqnsContainer;
    eqnsContainer.numEquations = 0;
    eqnsContainer.eqns = (equation *) malloc(sizeof(equation) * MAX_EQNS);
    int i = 0;
    char line[101], end[]="END";
    while(i<MAX_EQNS) {
        scanf("%s[^\n]\r\n", line);
        if(strcmp(end, line)==0) break;
        parseEquation(line, eqnsContainer.eqns+i);
        ++i;
    }
    if(i>=MAX_EQNS) {
        printf("ERROR: Increase the MAX_EQNS limit to solve for more than %d variables.\n", MAX_EQNS);
        exit(-1);
    }
    eqnsContainer.numEquations = i;
    int verifyStatus = verifyEquations(&eqnsContainer);
    if(verifyStatus<0) {
        printf("ERROR: Unable to verify equation's integrity\n");
        exit(-1);
    }
    return eqnsContainer;
}

void getCoeffsMatrix(familyOfEqns *eqns, float **mat) {
    int n = eqns->numEquations;
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            mat[i][j] = eqns->eqns[i].vars[j].val;
}

void getConstantsVector(familyOfEqns *eqns, float **vec) {
    int n = eqns->numEquations;
    for(int i=0; i<n; ++i)
        vec[i][0] = eqns->eqns[i].constant;
}

equation* getResultPairs(float **vec, coeff *vars, int n) {
    equation *solution = (equation*) malloc(sizeof(equation));
    solution->vars = (coeff*) malloc(sizeof(coeff)*n);
    solution->numVars = n;
    for(int i=0; i<n; ++i) {
        solution->vars[i].alpha = vars[i].alpha;
        solution->vars[i].val = vec[i][0];
    }
    return solution;
}

equation* solveEquations(familyOfEqns *eqnSet) {
    int n = (*eqnSet).numEquations;
    float **coeffsMatrix = getNewMatrix(n, n);
    float **constantsVector = getNewMatrix(n, 1);
    getCoeffsMatrix(eqnSet, coeffsMatrix);
    getConstantsVector(eqnSet, constantsVector);
    float **inverse = getNewMatrix(n, n);
    int status = getInverseMatrix(coeffsMatrix, inverse, n);
    if(status<0) {
        printf("The given equations cannot be solved.\n");
        exit(-1);
    }
    float **result = multiplyMatrices(inverse, constantsVector, n, n, n, 1);
    equation *solution = getResultPairs(result, eqnSet->eqns[0].vars, eqnSet->eqns[0].numVars);
    free(result);
    free(coeffsMatrix);
    free(constantsVector);
    free(inverse);
    return solution;
}


void printEquationCoeffs(equation *s) {
    for(int i=0; i< (s->numVars); ++i) 
        printf("%c = %0.2f\n", s->vars[i].alpha, s->vars[i].val);
}

int main() {
    familyOfEqns eqnSet = readEquations();
    equation *solution = solveEquations(&eqnSet);
    printf("\nSolution:\n");
    printEquationCoeffs(solution);
    free(solution);
    //TODO: Release all dynamic allocations properly
    return 0;
}
