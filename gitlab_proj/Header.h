#pragma once
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>

void VectorPlusVector(double* a, double* b, double* res, int N)
{
    for (int i = 0; i < N; i++)
        res[i] = a[i] + b[i];
}

double ScalarMultiplication(double* a, double* b, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
        sum += a[i] * b[i];
    return sum;
}

void VectorXVector(double* a, double* b, double* res, int N)
{
    for (int i = 0; i < N; i++)
        res[i] = a[i] * b[i];
}

void CoeffXVector(double* vec, double k, double* res, int N)
{
    for (int i = 0; i < N; i++)
        res[i] = k * vec[i];
}

void MatrixXVector(int* ia, int* ja, double* di, double* al, double* au, double* vec, double* res, int N)
{
    for (int i = 0; i < N; i++)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        res[i] = di[i] * vec[i];
        for (int k = i0; k < i1; k++)
        {
            int j = ja[k];
            res[i] += al[k] * vec[j];
            res[j] += au[k] * vec[i];
        }
    }
}
void Gilbert(double** g, int N, double* vec)
{
    for (int i = 1; i <= N; i++)
    {
        vec[i - 1] = 0;
        for (int j = 1; j <= N; j++)
        {
            g[i - 1][j - 1] = 1 / double(i + j - 1);
            vec[i - 1] += g[i - 1][j - 1] * (j);
        }
    }
}

void profile(double** g, int* ia, int* ja, double* al, double* au, double* di, int N)
{
    ia[0] = 0;
    for (int i = 0; i < N; i++)
    {
        ia[i + 1] = ia[i] + i;
        for (int j = 0; j < i; j++)
        {
            al[ia[i] + j] = g[i][j];
            au[ia[i] + j] = al[ia[i] + j];
            ja[ia[i] + j] = j;
        }
        di[i] = g[i][i];
    }
}
