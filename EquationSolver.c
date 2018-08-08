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
#define MAX_EQNS 52

char errorCodes[][100] = {
    "Equations are verified.",
    "No equations found.",
    "Infinite no. of solutions exists.",
    "More no. of equations given than number of variables."
};

typedef struct {
    float coeffs[MAX_EQNS];
    char isValid[MAX_EQNS];
    float constant;
    int numVars;
} equation;

typedef struct {
    equation *eqns;
    int numEquations;
    char isVarPresent[MAX_EQNS];
    int varsEncountered;
    char *varsToSolve;
} familyOfEqns;

int isUpperCase(char c) {
    return (c >= 'A' && c <= 'Z');
}

int isLowerCase(char c) {
    return (c >= 'a' && c <= 'z');
}

int isChar(char c) {
    return (isLowerCase(c) || isUpperCase(c));
}

int isDigit(char c) {
    return (c >= '0' && c <='9');
}

int encodeVar(char c) {
    if(isLowerCase(c))
        return c-'a';
    if(isUpperCase(c))
        return c-'A'+26;
    return -1;
}

char decodeVar(int n) {
    if(n>=0 && n<26) return 'a'+n;
    if(n>=26 && n<52) return 'A'+n-26;
    return '\0';
}

int countVariables(char *line) {
    int count = 0;
    for(int i=0; line[i]!='\0'; ++i)
        if(isChar(line[i])) ++count;
    return count;
}

void storeCoeff(equation *eqn, char c, float val) {
    int location = encodeVar(c);
    eqn->coeffs[location] += val;
    eqn->isValid[location] = 1;
}

float getCoeff(equation *eqn, char c) {
    int location = encodeVar(c);
    return eqn->coeffs[location];
}

void markVariablePresence(familyOfEqns *eqnsContainer, char c) {
    int location = encodeVar(c);
    if(eqnsContainer->isVarPresent[location]) return;
    eqnsContainer->isVarPresent[location] = 1;
    eqnsContainer->varsToSolve[eqnsContainer->varsEncountered] = c;
    ++(eqnsContainer->varsEncountered);
}

//parses number from str[i] into val and return no. of chars processed
int parseNumber(char *str, int i, float *val, float defaultVal) {
    int sign = 1, floatLocation = -1, parseLen = 0;
    *val = 0;
    while(1) {
        if(isDigit(str[i])) {
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

// Parses and stores coefficient & corresponding variable
// name at str[i] into 'var' (a struct container for both)
int parseToken(char *str, int i, equation *eqn) {
    if(str[i]=='=') return 0;
    float value;
    int parseLen = parseNumber(str, i, &value, 1.0);
    if(parseLen==-1) {
        printf("\nProblem parsing the equation at index %d\n", i);
    }
    i += parseLen;
    if(isChar(str[i])) {
        storeCoeff(eqn, str[i], value);
        return parseLen+1;
    } else {
        return -1;
    }
}

// Checks if given character is present in string str
int isCharInString(char c, char *str) {
    for(int i=0; str[i] != '\0'; ++i) {
        if(c == str[i]) return 1;
    }
    return 0;
}

void removeJunkCharacters(char *str) {
    char removeChars[] = {' ', '\t', '\r', '\n'};
    int k=0;
    for(int i=0; str[i] != '\0'; ++i) {
        if(!isCharInString(str[i], removeChars))
            str[k++] = str[i];
    }
    str[k] = '\0';
}

// Checks if equation contains only valid characters
// Returns 0 if it contains illegal characters
int checkEquationValidity(char* eqn, familyOfEqns *eqnsContainer) {
    int wasEqualsEncountered = 0;
    char allowedArithmeticChars[] = {'+', '-', '.'};
    for(int i=0; eqn[i] != '\0'; ++i) {
        if((isChar(eqn[i]) && !isChar(eqn[i+1]) && !wasEqualsEncountered)) {
            markVariablePresence(eqnsContainer, eqn[i]);
            continue;
        } else if(isDigit(eqn[i])) {
            continue;
        } else if(isCharInString(eqn[i], allowedArithmeticChars) && !isCharInString(eqn[i+1], allowedArithmeticChars))
            continue;
        else if(eqn[i]=='=' && !wasEqualsEncountered) {
            wasEqualsEncountered = 1;
            continue;
        } else {
            return 0;
        }
    }
    return (wasEqualsEncountered);
}

void parseEquation(char *line, familyOfEqns *eqnsContainer, int i) {
    equation *eqn = eqnsContainer->eqns + i;
    removeJunkCharacters(line);
    int status = checkEquationValidity(line, eqnsContainer);
    if(!status) {
        printf("\nERROR: Equation in illegal format.\n");
        exit(-1);
    }
    eqn->numVars = countVariables(line);
    memset(eqn->coeffs, 0, MAX_EQNS*sizeof(eqn->coeffs[0]));
    memset(eqn->isValid, 0, MAX_EQNS*sizeof(eqn->isValid[0]));
    int index = 0, parseLen, varIndex = 0;
    while(line[index] != '\0') {
        parseLen = parseToken(line, index, eqn);
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

int verifyEqns(familyOfEqns *eqns) {
    if(eqns->numEquations==0) return 1;
    if(eqns->numEquations < eqns->varsEncountered) return 2;
    if(eqns->numEquations > eqns->varsEncountered) return 3;
    eqns->varsToSolve[eqns->varsEncountered] = '\0';
    return 0;
}

familyOfEqns readEquations() {
    familyOfEqns eqnsContainer;
    eqnsContainer.numEquations = 0;
    eqnsContainer.eqns = (equation *) malloc(sizeof(equation) * MAX_EQNS);
    memset(eqnsContainer.isVarPresent, 0, MAX_EQNS*sizeof(eqnsContainer.isVarPresent[0]));
    eqnsContainer.varsToSolve = malloc(MAX_EQNS);
    eqnsContainer.varsEncountered = 0;
    int i = 0, baseEquation=0;
    char line[101], end[]="END";
    while(i<MAX_EQNS) {
        scanf("%[^\n]", line);
        getchar();
        if(strcmp(end, line)==0) break;
        parseEquation(line, &eqnsContainer, i);
        if(eqnsContainer.eqns[i].numVars > eqnsContainer.eqns[baseEquation].numVars)
            baseEquation = i;
        ++i;
    }
    printf("\n");
    if(i>=MAX_EQNS) {
        printf("ERROR: Cannot solve for more than %d variables.\n", MAX_EQNS);
        exit(-1);
    }
    eqnsContainer.numEquations = i;
    int verifyStatus = verifyEqns(&eqnsContainer);
    if(verifyStatus) {
        printf("%s\n", errorCodes[verifyStatus]);
        printf("The given equations cannot be solved.\n");
        exit(-1);
    }
    return eqnsContainer;
}

void getCoeffsMatrix(familyOfEqns *eqns, float **mat) {
    int n = eqns->numEquations;
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            mat[i][j] = getCoeff(eqns->eqns+i, eqns->varsToSolve[j]);
}

void getConstantsVector(familyOfEqns *eqns, float **vec) {
    int n = eqns->numEquations;
    for(int i=0; i<n; ++i)
        vec[i][0] = eqns->eqns[i].constant;
}

float** solveEquations(familyOfEqns *eqnSet) {
    int n = (*eqnSet).numEquations;
    float **coeffsMatrix = getNewMatrix(n, n);
    float **constantsVector = getNewMatrix(n, 1);
    getCoeffsMatrix(eqnSet, coeffsMatrix);
    getConstantsVector(eqnSet, constantsVector);
    float **inverse = getNewMatrix(n, n);
    int status = getInverseMatrix(coeffsMatrix, inverse, n);
    if(status<0) {
        printf("The given equations have no solution.\n");
        exit(-1);
    }
    float **result = multiplyMatrices(inverse, constantsVector, n, n, n, 1);
    free(coeffsMatrix);
    free(constantsVector);
    free(inverse);
    return result;
}

void printSolution(float **solution, familyOfEqns *eqnSet) {
    for(int i=0; i<eqnSet->varsEncountered; ++i)
        printf("%c: %0.2f\n", eqnSet->varsToSolve[i], solution[i][0]);
}

int main() {
    familyOfEqns eqnSet = readEquations();
    float **solution = solveEquations(&eqnSet);
    printf("Solution:\n");
    printSolution(solution, &eqnSet);
    free(solution);
    //TODO: Release all dynamic allocations properly
    return 0;
}
