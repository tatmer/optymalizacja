#pragma once

#include"ode_solver.h"
#include<cmath>
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

// --- Laboratorium 1 ---
matrix ff1T(matrix, matrix = NAN, matrix = NAN); // Testowa funkcja celu 
matrix ff1R(matrix, matrix = NAN, matrix = NAN); // Rzeczywista funkcja celu (zbiorniki) 
matrix df1(double, matrix, matrix = NAN, matrix = NAN); // Równania różniczkowe (zbiorniki) 

// --- Labolatorium 2 ---

matrix ff_lab2_T(matrix, matrix = NAN, matrix = NAN);
// Zadanie 5b - problem rzeczywisty (robot)
matrix ff_lab2_R(matrix, matrix = NAN, matrix = NAN);
matrix df_lab2_R(double, matrix, matrix = NAN, matrix = NAN);

// --- Labolatorium 3 ---
// --- Laboratorium 3 ---
matrix ff3T(matrix, matrix = NAN, matrix = NAN); // Funkcja z ograniczeniami
matrix ff3R(matrix, matrix = NAN, matrix = NAN); // Rzeczywista funkcja celu (piłka)
matrix df3(double, matrix, matrix = NAN, matrix = NAN); // Równania ruchu (piłka)

// --- Labolatorium 4 ---
// --- Funkcja Testowa (Wielomian) ---
// Zwraca wartość funkcji celu
matrix func_test(matrix x, matrix ud1, matrix ud2);
// Zwraca wektor gradientu
matrix grad_test(matrix x, matrix ud1, matrix ud2);
// Zwraca macierz Hessego (dla metody Newtona)
matrix hess_test(matrix x, matrix ud1, matrix ud2);

// --- Problem Rzeczywisty (Klasyfikator) ---
// Funkcja kosztu logistycznego
matrix func_real(matrix theta, matrix X, matrix Y);
// Gradient funkcji kosztu
matrix grad_real(matrix theta, matrix X, matrix Y);