#include "opt_alg.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm> // do min/max

using namespace std;

// --- Zmienne globalne ---
int g_calls_cnt = 0;
int H_calls_cnt = 0;
std::vector<matrix>* g_path_trace = nullptr; // Wskaźnik do zapisu trajektorii (Lab 4)

// --- Pomocnicze funkcje do metody Fibonacciego (Lab 1) ---
const int FIB_MAX = 100;
double F[FIB_MAX + 1];
bool fib_init = false;

void calculate_fib() {
    if (fib_init) return;
    F[0] = 0; F[1] = 1;
    for (int i = 2; i <= FIB_MAX; ++i) F[i] = F[i - 1] + F[i - 2];
    fib_init = true;
}

// --- Pomocnicza funkcja: Złoty podział (Line Search dla Lab 4) ---
double line_search_golden(matrix(*ff)(matrix, matrix, matrix), matrix x, matrix d, matrix ud1, matrix ud2) {
    double a = 0.0;
    double b = 1.0; 
    double epsilon = 1e-7;
    double alpha = (sqrt(5.0) - 1.0) / 2.0;

    double c = b - alpha * (b - a);
    double d_pt = a + alpha * (b - a);

    double fc = ff(x + c * d, ud1, ud2)(0);
    double fd = ff(x + d_pt * d, ud1, ud2)(0);
    solution::f_calls += 2; 

    while ((b - a) > epsilon) {
        if (fc < fd) {
            b = d_pt; d_pt = c; fd = fc;
            c = b - alpha * (b - a);
            fc = ff(x + c * d, ud1, ud2)(0);
            solution::f_calls++;
        } else {
            a = c; c = d_pt; fc = fd;
            d_pt = a + alpha * (b - a);
            fd = ff(x + d_pt * d, ud1, ud2)(0);
            solution::f_calls++;
        }
    }
    return (a + b) / 2.0;
}


// ===========================================================================
//                          IMPLEMENTACJE ALGORYTMÓW
// ===========================================================================

// --- Monte Carlo (Lab 0) ---
solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y(0) < epsilon) { Xopt.flag = 1; break; }
            if (solution::f_calls > Nmax) { Xopt.flag = 0; break; }
        }
        return Xopt;
    } catch (string ex_info) { throw ("solution MC(...):\n" + ex_info); }
}

// --- Metoda Ekspansji (Lab 1) ---
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
    try {
        double* p = new double[2];
        solution X0(x0), X1(x0 + d);
        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);

        if (X1.y(0) == X0.y(0)) { p[0] = x0; p[1] = x0 + d; return p; }
        if (X1.y(0) > X0.y(0)) {
            d = -d;
            X1.x = x0 + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y(0) >= X0.y(0)) { p[0] = x0 + d; p[1] = x0 - d; return p; }
        }

        solution X2;
        int i = 0;
        while (true) {
            i++;
            if (solution::f_calls > Nmax) break;
            X2.x = x0 + pow(alpha, i) * d;
            X2.fit_fun(ff, ud1, ud2);
            if (X2.y(0) >= X1.y(0)) break;
            X0 = X1; X1 = X2;
        }
        if (d > 0) { p[0] = X0.x(0); p[1] = X2.x(0); }
        else { p[0] = X2.x(0); p[1] = X0.x(0); }
        return p;
    } catch (string ex_info) { throw ("double* expansion(...):\n" + ex_info); }
}

// --- Metoda Fibonacciego (Lab 1) ---
solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
    try {
        calculate_fib();
        solution Xopt;
        double L = b - a;
        int n = 0;
        while (F[n] <= L / epsilon) n++;
        
        // Zabezpieczenie przed wyjściem poza tablicę
        if (n > FIB_MAX) n = FIB_MAX;

        double x1 = a + (double)F[n - 2] / F[n] * (b - a);
        double x2 = a + (double)F[n - 1] / F[n] * (b - a);
        
        matrix m1(1, 1, x1), m2(1, 1, x2);
        solution X1(m1), X2(m2);
        
        X1.fit_fun(ff, ud1, ud2);
        X2.fit_fun(ff, ud1, ud2);

        for (int k = 0; k < n - 3; ++k) {
            if (X1.y(0) < X2.y(0)) {
                b = x2;
                x2 = x1;
                X2 = X1;
                x1 = a + (double)F[n - k - 3] / F[n - k - 1] * (b - a);
                X1.x(0) = x1;
                X1.fit_fun(ff, ud1, ud2);
            } else {
                a = x1;
                x1 = x2;
                X1 = X2;
                x2 = a + (double)F[n - k - 2] / F[n - k - 1] * (b - a);
                X2.x(0) = x2;
                X2.fit_fun(ff, ud1, ud2);
            }
        }
        Xopt = (X1.y(0) < X2.y(0)) ? X1 : X2;
        Xopt.ud = matrix(2, 1);
        Xopt.ud(0) = a; Xopt.ud(1) = b;
        return Xopt;
    } catch (string ex_info) { throw ("solution fib(...):\n" + ex_info); }
}

// --- Metoda Lagrange'a (Lab 1) ---
solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
    try {
        solution Xopt;
        double x1 = a, x2 = b;
        double x3 = a + gamma * (b - a);
        
        matrix m1(1, 1, x1), m2(1, 1, x2), m3(1, 1, x3);
        solution X1(m1), X2(m2), X3(m3);
        
        X1.fit_fun(ff, ud1, ud2);
        X2.fit_fun(ff, ud1, ud2);
        X3.fit_fun(ff, ud1, ud2);

        while (true) {
            double L = X1.y(0) * (pow(X2.x(0), 2) - pow(X3.x(0), 2)) +
                       X2.y(0) * (pow(X3.x(0), 2) - pow(X1.x(0), 2)) +
                       X3.y(0) * (pow(X1.x(0), 2) - pow(X2.x(0), 2));
            double M = X1.y(0) * (X2.x(0) - X3.x(0)) +
                       X2.y(0) * (X3.x(0) - X1.x(0)) +
                       X3.y(0) * (X1.x(0) - X2.x(0));

            if (abs(M) <= 1e-10) { Xopt = X3; Xopt.flag = 0; break; }
            double x4 = 0.5 * L / M;
            
            if (x4 < a || x4 > b) { Xopt = X3; Xopt.flag = 0; break; }

            matrix m4(1, 1, x4);
            solution X4(m4);
            X4.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) { Xopt = X3; Xopt.flag = 0; break; }
            if (abs(X4.x(0) - X3.x(0)) < epsilon) { Xopt = X4; Xopt.flag = 1; break; }

            if (X4.y(0) < X3.y(0)) {
                if (X4.x(0) < X3.x(0)) {
                    X2 = X3; X3 = X4;
                } else {
                    X1 = X3; X3 = X4;
                }
            } else {
                if (X4.x(0) < X3.x(0)) {
                    X1 = X4;
                } else {
                    X2 = X4;
                }
            }
        }
        Xopt.ud = matrix(2, 1);
        Xopt.ud(0) = X1.x(0); Xopt.ud(1) = X2.x(0);
        return Xopt;
    } catch (string ex_info) { throw ("solution lag(...):\n" + ex_info); }
}

// --- Metoda Hooke'a-Jeevesa (Lab 2) ---
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
    int n = get_len(XB.x);
    for (int i = 0; i < n; ++i) {
        XB.x(i) += s;
        XB.fit_fun(ff, ud1, ud2);
        double y_plus = XB.y(0);
        XB.x(i) -= s;

        if (y_plus < XB.y(0)) {
            XB.x(i) += s;
        } else {
            XB.x(i) -= s;
            XB.fit_fun(ff, ud1, ud2);
            double y_minus = XB.y(0);
            XB.x(i) += s;
            
            if (y_minus < XB.y(0)) {
                XB.x(i) -= s;
                XB.fit_fun(ff, ud1, ud2);
            }
        }
    }
    return XB;
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try {
        bool record = (!isnan(ud1(0,0)) && ud1(0,0) == 1.0);
        ofstream Traj;
        if(record) Traj.open("trajektoria_HJ.csv");

        solution X(x0), XB(x0);
        X.fit_fun(ff, ud1, ud2);
        XB.fit_fun(ff, ud1, ud2);

        while (true) {
            if(record) Traj << XB.x(0) << ";" << XB.x(1) << endl;
            if (solution::f_calls > Nmax) { X.flag = 0; break; }

            X = HJ_trial(ff, XB, s, ud1, ud2);
            
            if (X.y(0) < XB.y(0)) {
                while (true) {
                    if(record) Traj << X.x(0) << ";" << X.x(1) << endl;
                    if (solution::f_calls > Nmax) break;
                    
                    solution XB_old = XB;
                    XB = X;
                    X.x = 2.0 * XB.x - XB_old.x;
                    X.fit_fun(ff, ud1, ud2);
                    X = HJ_trial(ff, X, s, ud1, ud2);
                    
                    if (X.y(0) >= XB.y(0)) break;
                }
                X = XB;
            } else {
                s *= alpha;
            }
            if (s < epsilon) { X.flag = 1; break; }
        }
        if(record) Traj.close();
        return X;
    } catch (string ex_info) { throw ("solution HJ(...):\n" + ex_info); }
}

// --- Metoda Rosenbrocka (Lab 2) ---
solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try {
        solution X(x0);
        X.fit_fun(ff, ud1, ud2);
        int n = get_len(x0);
        matrix d = ident_mat(n);
        matrix lambda(n, 1), p(n, 1);
        matrix s = s0;

        while (true) {
            if (solution::f_calls > Nmax) { X.flag = 0; break; }
            bool change = true;
            for(int i=0; i<n; ++i) {
                if (abs(s(i)) < epsilon) change = false;
            }
            if (!change) { X.flag=1; break; }

            for (int i = 0; i < n; ++i) {
                solution Xt(X.x + s(i) * d[i]);
                Xt.fit_fun(ff, ud1, ud2);
                if (Xt.y(0) < X.y(0)) {
                    X = Xt;
                    lambda(i) += s(i);
                    s(i) *= alpha;
                } else {
                    s(i) *= -beta;
                    p(i) += 1;
                }
            }
        }
        return X;
    } catch (string ex_info) { throw ("solution Rosen(...):\n" + ex_info); }
}

// --- Funkcja Kary (Lab 3) ---
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try {
        double alpha = dc;
        solution X(x0);
        double c_curr = c;
        
        while (true) {
            if (get_len(ud2) < 1) ud2 = matrix(1, 1);
            ud2(0) = c_curr;
            
            solution X_new = sym_NM(ff, X.x, 0.5, 1.0, 0.5, 2.0, 0.5, epsilon, Nmax, ud1, ud2);
            
            if (solution::f_calls > Nmax) { X = X_new; X.flag = 0; break; }
            if (norm(X.x - X_new.x) < epsilon) { X = X_new; X.flag = 1; break; }
            
            X = X_new;
            c_curr *= alpha;
        }
        return X;
    } catch (string ex_info) { throw ("solution pen(...):\n" + ex_info); }
}

// --- Metoda Sympleksu Neldera-Meada (Lab 3) ---
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try {
        int n = get_len(x0);
        vector<solution> p(n + 1);
        
        p[0].x = x0;
        p[0].fit_fun(ff, ud1, ud2);

        for (int i = 1; i <= n; ++i) {
            matrix ei(n, 1); ei(i - 1) = 1.0;
            p[i].x = p[0].x + s * ei;
            p[i].fit_fun(ff, ud1, ud2);
        }

        while (true) {
            if (solution::f_calls > Nmax) return p[0];

            sort(p.begin(), p.end(), [](const solution& a, const solution& b) {
                return a.y(0) < b.y(0);
            });
            
            int best = 0;
            int worst = n;
            
            if (norm(p[worst].x - p[best].x) < epsilon) { p[best].flag = 1; return p[best]; }

            matrix p_bar(n, 1);
            for (int i = 0; i < n; ++i) p_bar = p_bar + p[i].x;
            p_bar = p_bar / (double)n;

            solution p_r;
            p_r.x = p_bar + alpha * (p_bar - p[worst].x);
            p_r.fit_fun(ff, ud1, ud2);

            if (p_r.y(0) < p[best].y(0)) {
                solution p_e;
                p_e.x = p_bar + gamma * (p_r.x - p_bar);
                p_e.fit_fun(ff, ud1, ud2);
                if (p_e.y(0) < p_r.y(0)) p[worst] = p_e;
                else p[worst] = p_r;
            } else {
                if (p_r.y(0) < p[n - 1].y(0)) {
                    p[worst] = p_r;
                } else {
                    solution p_c;
                    p_c.x = p_bar + beta * (p[worst].x - p_bar);
                    p_c.fit_fun(ff, ud1, ud2);
                    
                    if (p_c.y(0) < p[worst].y(0)) {
                        p[worst] = p_c;
                    } else {
                        for (int i = 1; i <= n; ++i) {
                            p[i].x = p[best].x + delta * (p[i].x - p[best].x);
                            p[i].fit_fun(ff, ud1, ud2);
                        }
                    }
                }
            }
        }
    } catch (string ex_info) { throw ("solution sym_NM(...):\n" + ex_info); }
}

// ===========================================================================
//                          LAB 4: SD, CG, NEWTON
// ===========================================================================

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    solution Xopt;
    Xopt.x = x0;
    Xopt.fit_fun(ff, ud1, ud2);

    // Zapis punktu startowego
    if (g_path_trace) g_path_trace->push_back(Xopt.x);

    matrix d(x0);
    solution X1;

    while (solution::f_calls < Nmax)
    {
        matrix grad = gf(Xopt.x, ud1, ud2);
        g_calls_cnt++;

        if (norm(grad) < epsilon) break;
        d = -grad;

        double h;
        if (h0 < 0) h = line_search_golden(ff, Xopt.x, d, ud1, ud2);
        else h = h0;

        X1.x = Xopt.x + h * d;
        X1.fit_fun(ff, ud1, ud2);

        if (g_path_trace) g_path_trace->push_back(X1.x);

        if (norm(X1.x - Xopt.x) < epsilon) { Xopt = X1; break; }
        Xopt = X1;
    }
    return Xopt;
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    solution Xopt;
    Xopt.x = x0;
    Xopt.fit_fun(ff, ud1, ud2);
    
    // Zapis punktu startowego
    if (g_path_trace) g_path_trace->push_back(Xopt.x);

    matrix grad = gf(Xopt.x, ud1, ud2);
    g_calls_cnt++;

    matrix d = -grad;
    matrix grad_old = grad;
    solution X1;

    while (solution::f_calls < Nmax)
    {
        double h;
        if (h0 < 0) h = line_search_golden(ff, Xopt.x, d, ud1, ud2);
        else h = h0;

        X1.x = Xopt.x + h * d;
        X1.fit_fun(ff, ud1, ud2);

        if (g_path_trace) g_path_trace->push_back(X1.x);

        if (norm(X1.x - Xopt.x) < epsilon) { Xopt = X1; break; }
        Xopt = X1;

        matrix grad_new = gf(Xopt.x, ud1, ud2);
        g_calls_cnt++;

        if (norm(grad_new) < epsilon) break;

        double beta = pow(norm(grad_new), 2) / pow(norm(grad_old), 2);
        d = -grad_new + beta * d;
        grad_old = grad_new;
    }
    return Xopt;
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    solution Xopt;
    Xopt.x = x0;
    Xopt.fit_fun(ff, ud1, ud2);

    // Zapis punktu startowego
    if (g_path_trace) g_path_trace->push_back(Xopt.x);

    solution X1;

    while (solution::f_calls < Nmax)
    {
        matrix grad = gf(Xopt.x, ud1, ud2);
        g_calls_cnt++;
        if (norm(grad) < epsilon) break;

        matrix H = Hf(Xopt.x, ud1, ud2);
        H_calls_cnt++;

        matrix d = -inv(H) * grad;

        double h;
        if (h0 < 0) h = line_search_golden(ff, Xopt.x, d, ud1, ud2);
        else {
             h = 1.0; 
             if (h0 != 1.0) h = h0;
        }

        X1.x = Xopt.x + h * d;
        X1.fit_fun(ff, ud1, ud2);
        
        if (g_path_trace) g_path_trace->push_back(X1.x);

        if (norm(X1.x - Xopt.x) < epsilon) { Xopt = X1; break; }
        Xopt = X1;
    }
    return Xopt;
}

// --- Funkcje niewykorzystane ---
solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }