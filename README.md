# Linear Equations Solver

A simple program to solve the linear equations given via stdin.

## Steps to build

```
gcc EquationSolver.c -o EquationSolver.o -lm
```

## Sample I/O

```
$ ./EquationSolver.o 
x + y = 2
2x+3y = 5        
END

Solution:
x: 1.00
y: 1.00
```

## Instructions

- Equations have to be in format `ax+by+...=k` where `x,y,...` are the variables to solve.
- Each equation has to be separated by a new line.
- After entering all equations, type `END` to stop reading input.

## Features

- Can have any number of variables (MAX 52, since each variable can only be an alphabet including upper & lower cases.
- Variables in the equations can be in any order.

## How the equations are solved

The equations are constructed into matrices, to obtain a matrix equation of the form:
```
AX = C
```

where
- `A` is the coefficient matrix
- `X` is the vector of variables
- `C` is the vector of constants (from RHS of each equation)

The solution is found just by calculating
```
X = A⁻¹.C
```

where 
- `A⁻¹` is the inverse of the coefficient matrix

## Footnotes

This was done by me as an assignment while training at [MulticoreWare](https://multicorewareinc.com/).  
Feel free to contribute. 
