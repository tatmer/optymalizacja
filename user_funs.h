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

// --- Laboratorium 2 ---
matrix ff_lab2_T(matrix, matrix = NAN, matrix = NAN); // Funkcja celu dla laboratorium 2