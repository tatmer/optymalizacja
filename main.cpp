/*********************************************
Kod stanowi uzupełnienie materiałów do ćwiczeń
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include"user_funs.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<random>
#include<cmath>
#include<iomanip>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- ZMIENNE EXTERN ---
// Pozwala używać liczników z opt_alg.cpp bez modyfikacji opt_alg.h
extern int g_calls_cnt;
extern int H_calls_cnt;
extern std::vector<matrix>* g_path_trace;
// ----------------------

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
    try
    {
        lab4(); // Uruchomienie Lab 4
    }
    catch (string EX_INFO)
    {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    return 0;
}

// Funkcja pomocnicza dla Lab 4
double calculate_accuracy(matrix theta, matrix X, matrix Y) {
    int m = 100;
    int correct = 0;
    for (int i = 0; i < m; ++i) {
        double z = 0.0;
        for (int j = 0; j < 3; ++j) z += theta(j) * X(j, i);
        double h = 1.0 / (1.0 + exp(-z));
        int prediction = (h >= 0.5) ? 1 : 0;
        if (prediction == (int)Y(0, i)) correct++;
    }
    return (double)correct / m * 100.0;
}

// --- POZOSTAŁE LABORATORIA ---
void lab0()
{
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1; a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    Nmax = 1000; epsilon = 1e-2; lb = 0; ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
    matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix(); Y[1].~matrix();
}

void lab1() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distr_x0(-100.0, 100.0);
    std::uniform_real_distribution<> distr_da_cm2(1.0, 100.0);

    double epsilon = 1e-4;
    int Nmax = 1000;

    cout << "5. a) Optymalizacja funkcji testowej..." << endl;
    vector<double> alphas = { 1.2, 1.5, 2.0 };
    ofstream file_test("wyniki_testowe.csv");
    file_test << "Alpha;x0;A_fib;B_fib;Xopt_fib;Yopt_fib;Iter_fib;A_lag;B_lag;Xopt_lag;Yopt_lag;Iter_lag" << endl;

    for (double alpha : alphas) {
        cout << "  -> Przetwarzanie dla alpha = " << alpha << endl;
        for (int i = 0; i < 100; ++i) {
            double x0 = distr_x0(gen);
            double d0 = 1.0; 
            double* interval = expansion(ff1T, x0, d0, alpha, Nmax);
            double a = interval[0];
            double b = interval[1];
            delete[] interval;

            solution::clear_calls();
            solution opt_fib = fib(ff1T, a, b, epsilon, NAN, NAN);
            int iter_fib = solution::f_calls;

            solution::clear_calls();
            solution opt_lag = lag(ff1T, a, b, epsilon, 0.5, Nmax, NAN, NAN);
            int iter_lag = solution::f_calls;

            file_test << alpha << ";" << x0 << ";"
                << opt_fib.ud(0) << ";" << opt_fib.ud(1) << ";" << opt_fib.x(0) << ";" << opt_fib.y(0) << ";" << iter_fib << ";"
                << opt_lag.ud(0) << ";" << opt_lag.ud(1) << ";" << opt_lag.x(0) << ";" << opt_lag.y(0) << ";" << iter_lag << endl;
        }
    }
    file_test.close();

    cout << "\n5. b) Optymalizacja problemu rzeczywistego (zbiorniki)..." << endl;
    matrix TB_target(1, 1, 50.0);
    
    double a_real = 1.0; 
    double b_real = 100.0;

    solution::clear_calls();
    solution opt_fib_real = fib(ff1R, a_real, b_real, epsilon, TB_target, NAN);
    cout << "  -> Fibonacciego D_A: " << opt_fib_real.x(0) << " cm^2, Blad: " << opt_fib_real.y(0) << ", Iter: " << solution::f_calls << endl;

    solution::clear_calls();
    solution opt_lag_real = lag(ff1R, a_real, b_real, epsilon, 0.5, Nmax, TB_target, NAN);
    cout << "  -> Lagrange'a D_A: " << opt_lag_real.x(0) << " cm^2, Blad: " << opt_lag_real.y(0) << ", Iter: " << solution::f_calls << endl;

    solution opt_best = (opt_fib_real.y(0) < opt_lag_real.y(0)) ? opt_fib_real : opt_lag_real;

    cout << "\n  -> Przeprowadzam symulacje dla optymalnego D_A = " << opt_best.x(0) << " cm^2" << endl;
    double DA_best_m2 = opt_best.x(0) / 10000.0;
    matrix Y0(3, 1);
    Y0(0) = 5.0; Y0(1) = 1.0; Y0(2) = 20.0;
    matrix ud2_solver(1, 1, DA_best_m2);

    matrix* S = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, NAN, ud2_solver);
    int N = get_len(S[0]);

    ofstream file_sim("symulacja_zbiorniki.csv");
    file_sim << "Czas;VA;VB;TB" << endl;
    double max_TB = 0.0;
    for (int i = 0; i < N; ++i) {
        file_sim << S[0](i) << ";" << S[1](0, i) << ";" << S[1](1, i) << ";" << S[1](2, i) << endl;
        if (S[1](2, i) > max_TB) max_TB = S[1](2, i);
    }
    file_sim.close();
    delete[] S;
    cout << "  -> Maksymalna osiagnieta temperatura TB: " << max_TB << " st. C" << endl;
}

void lab2()
{
    ofstream Sout("wyniki_lab2_partA.csv");
    if (!Sout.is_open()) { cerr << "Blad zapisu partA" << endl; return; }

    cout << "--- Zadanie 5a: Funkcja Testowa (Statystyka) ---" << endl;
    Sout << fixed << setprecision(6);
    cout << fixed << setprecision(6);

    int N_runs = 100;
    int Nmax = 20000;
    double epsilon = 1e-3;
    matrix lb(2, 1, -1.0);
    matrix ub(2, 1, 1.0);

    std::vector<matrix> start_points;
    for (int i = 0; i < N_runs; ++i)
    {
        matrix x0 = rand_mat(2, 1);
        x0(0) = (ub(0) - lb(0)) * x0(0) + lb(0); 
        x0(1) = (ub(1) - lb(1)) * x0(1) + lb(1); 
        start_points.push_back(x0);
    }

    double s_hj_all[] = { 1.0, 0.1, 0.01 };
    double alpha_hj = 0.5;
    matrix s0_rosen_all[] = { matrix(2, 1, 0.5), matrix(2, 1, 0.1), matrix(2, 1, 0.05) };
    double alpha_rosen = 2.0;
    double beta_rosen = 0.5;

    for (double s_hj : s_hj_all)
    {
        Sout << "\nMetoda Hooke'a-Jeevesa, s = " << s_hj << "\n";
        Sout << "x1_start;x2_start;x1_opt;x2_opt;y_opt;f_calls;flag\n";
        int global_found = 0;
        for (const auto& x0 : start_points)
        {
            solution::clear_calls();
            solution opt = HJ(ff_lab2_T, x0, s_hj, alpha_hj, epsilon, Nmax);
            Sout << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";" << opt.y << ";" << solution::f_calls << ";" << opt.flag << "\n";
            if (opt.flag == 1 && norm(opt.x) < 0.1 && opt.y < 0.01) global_found++;
        }
        cout << "HJ (s=" << s_hj << "): Globalne = " << global_found << "/" << N_runs << endl;
    }

    for (const auto& s0_rosen : s0_rosen_all)
    {
        double s_val = s0_rosen(0);
        Sout << "\nMetoda Rosenbrocka, s = " << s_val << "\n";
        Sout << "x1_start;x2_start;x1_opt;x2_opt;y_opt;f_calls;flag\n";
        int global_found = 0;
        for (const auto& x0 : start_points)
        {
            solution::clear_calls();
            solution opt = Rosen(ff_lab2_T, x0, s0_rosen, alpha_rosen, beta_rosen, epsilon, Nmax);
            Sout << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";" << opt.y << ";" << solution::f_calls << ";" << opt.flag << "\n";
            if (opt.flag == 1 && norm(opt.x) < 0.1 && opt.y < 0.01) global_found++;
        }
        cout << "Rosenbrock (s=" << s_val << "): Globalne = " << global_found << "/" << N_runs << endl;
    }
    Sout.close();

    ofstream SoutR("wyniki_lab2_partB.csv");
    if (!SoutR.is_open()) { cerr << "Blad zapisu partB" << endl; return; }

    cout << "\n--- Zadanie 5b: Problem Rzeczywisty ---" << endl;
    SoutR << fixed << setprecision(6);

    matrix x0_robot(2, 1); x0_robot(0) = 10.0; x0_robot(1) = 10.0;
    double s_robot = 0.5; 
    matrix s0_robot(2, 1, s_robot);

    solution::clear_calls();
    solution opt_HJ = HJ(ff_lab2_R, x0_robot, s_robot, alpha_hj, epsilon, Nmax);
    int hj_f_calls = solution::f_calls; 
    
    solution::clear_calls();
    solution opt_Ros = Rosen(ff_lab2_R, x0_robot, s0_robot, alpha_rosen, beta_rosen, epsilon, Nmax);
    
    SoutR << "Metoda;k1_start;k2_start;k1_opt;k2_opt;Q_min;f_calls\n";
    SoutR << "HJ;" << x0_robot(0) << ";" << x0_robot(1) << ";" << opt_HJ.x(0) << ";" << opt_HJ.x(1) << ";" << opt_HJ.y << ";" << hj_f_calls << "\n";
    SoutR << "Rosenbrock;" << x0_robot(0) << ";" << x0_robot(1) << ";" << opt_Ros.x(0) << ";" << opt_Ros.x(1) << ";" << opt_Ros.y << ";" << solution::f_calls << "\n";
    SoutR.close();

    solution best_opt = (opt_HJ.y < opt_Ros.y) ? opt_HJ : opt_Ros;
    matrix Y0(2, 1); 
    matrix* Y_sim = solve_ode(df_lab2_R, 0, 0.1, 100, Y0, best_opt.x);

    ofstream SimOut("symulacja_lab2.csv");
    if (SimOut.is_open())
    {
        SimOut << fixed << setprecision(6);
        int n_steps = get_len(Y_sim[0]);
        for (int i = 0; i < n_steps; ++i)
            SimOut << Y_sim[0](i, 0) << ";" << Y_sim[1](i, 0) << ";" << Y_sim[1](i, 1) << "\n";
        SimOut.close();
    }
    Y_sim[0].~matrix(); Y_sim[1].~matrix();

    matrix x0_traj = start_points[0];
    matrix flag_record(1, 1, 1.0);
    solution::clear_calls();
    HJ(ff_lab2_T, x0_traj, 0.5, alpha_hj, epsilon, Nmax, flag_record);
}

void lab3()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis_x(-10.0, 10.0);
    uniform_real_distribution<> dis_y(-10.0, 10.0);

    double a_vals[] = { 4.0, 4.4934, 5.0 };
    int N_runs = 100;
    int Nmax = 5000;
    double epsilon = 1e-3;

    ofstream fout_zew("wyniki_lab3_zew.csv");
    fout_zew << fixed << setprecision(5) << "a;x1_start;x2_start;x1_kon;x2_kon;y;L_wyw;Flaga;r\n";
    ofstream favg_zew("wyniki_lab3_srednie_zew.csv");
    favg_zew << fixed << setprecision(5) << "a;x1_sr;x2_sr;y_sr;L_wyw_sr;r_sr\n";

    ofstream fout_wew("wyniki_lab3_wew.csv");
    fout_wew << fixed << setprecision(5) << "a;x1_start;x2_start;x1_kon;x2_kon;y;L_wyw;Flaga;r\n";
    ofstream favg_wew("wyniki_lab3_srednie_wew.csv");
    favg_wew << fixed << setprecision(5) << "a;x1_sr;x2_sr;y_sr;L_wyw_sr;r_sr\n";

    cout << "--- Lab 3: Zadanie 5a (Testowa) ---" << endl;

    for (double a : a_vals)
    {
        cout << "Przetwarzanie a = " << a << "..." << endl;
        double zew_s_x1 = 0, zew_s_x2 = 0, zew_s_y = 0, zew_s_calls = 0, zew_s_r = 0; int zew_success = 0;
        double wew_s_x1 = 0, wew_s_x2 = 0, wew_s_y = 0, wew_s_calls = 0, wew_s_r = 0; int wew_success = 0;

        for (int i = 0; i < N_runs; ++i)
        {
            solution::clear_calls();
            matrix x0(2, 1);
            do {
                x0(0) = dis_x(gen);
                x0(1) = dis_y(gen);
            } while (x0(0) <= 1.001 || x0(1) <= 1.001 || norm(x0) >= a - 0.001);

            matrix x0_base = x0;
            matrix ud1(1, 1); ud1(0) = a;
            matrix clean_ud(2, 1); clean_ud(0) = 0; clean_ud(1) = 0; 

            solution::clear_calls();
            matrix ud2_zew(2, 1); ud2_zew(0) = 1.0; ud2_zew(1) = 1.0;
            solution opt_zew = pen(ff3T, x0_base, 1.0, 2.0, epsilon, Nmax, ud1, ud2_zew);

            double y_zew = ff3T(opt_zew.x, ud1, clean_ud)(0);
            double r_zew = norm(opt_zew.x);
            fout_zew << a << ";" << x0_base(0) << ";" << x0_base(1) << ";" << opt_zew.x(0) << ";" << opt_zew.x(1) << ";" << y_zew << ";" << solution::f_calls << ";" << opt_zew.flag << ";" << r_zew << "\n";

            if (opt_zew.flag == 1) { zew_s_x1 += opt_zew.x(0); zew_s_x2 += opt_zew.x(1); zew_s_y += y_zew; zew_s_calls += solution::f_calls; zew_s_r += r_zew; zew_success++; }

            solution::clear_calls();
            matrix ud2_wew(2, 1); ud2_wew(0) = 1.0; ud2_wew(1) = 2.0;
            solution opt_wew = pen(ff3T, x0_base, 1.0, 0.5, epsilon, Nmax, ud1, ud2_wew);

            double y_wew = ff3T(opt_wew.x, ud1, clean_ud)(0);
            double r_wew = norm(opt_wew.x);
            fout_wew << a << ";" << x0_base(0) << ";" << x0_base(1) << ";" << opt_wew.x(0) << ";" << opt_wew.x(1) << ";" << y_wew << ";" << solution::f_calls << ";" << opt_wew.flag << ";" << r_wew << "\n";

            if (opt_wew.flag == 1) { wew_s_x1 += opt_wew.x(0); wew_s_x2 += opt_wew.x(1); wew_s_y += y_wew; wew_s_calls += solution::f_calls; wew_s_r += r_wew; wew_success++; }
        }
        if (zew_success > 0) favg_zew << a << ";" << zew_s_x1 / zew_success << ";" << zew_s_x2 / zew_success << ";" << zew_s_y / zew_success << ";" << zew_s_calls / zew_success << ";" << zew_s_r / zew_success << "\n";
        if (wew_success > 0) favg_wew << a << ";" << wew_s_x1 / wew_success << ";" << wew_s_x2 / wew_success << ";" << wew_s_y / wew_success << ";" << wew_s_calls / wew_success << ";" << wew_s_r / wew_success << "\n";
    }
    fout_zew.close(); favg_zew.close();
    fout_wew.close(); favg_wew.close();

    cout << "\n--- Zadanie 5b: Problem rzeczywisty ---" << endl;
    matrix x0_real(2, 1); x0_real(0) = -5.0; x0_real(1) = 6.0;

    solution::clear_calls();
    matrix ud2_real(1, 1); ud2_real(0) = 1.0;
    solution opt_real = pen(ff3R, x0_real, 1.0, 2.0, 1e-3, 5000, NAN, ud2_real);

    matrix final_ud(1, 1); final_ud(0) = 0.0;
    double x_end_opt = -ff3R(opt_real.x, NAN, final_ud)(0);

    matrix Y0(4, 1); Y0(0) = 0.0; Y0(1) = opt_real.x(0); Y0(2) = 100.0; Y0(3) = 0.0;
    matrix omega_opt(1, 1); omega_opt(0) = opt_real.x(1);

    matrix* S = solve_ode(df3, 0, 0.01, 7.0, Y0, NAN, omega_opt);

    double x_at_50 = 0.0;
    int n = get_len(S[0]);
    for (int i = 0; i < n - 1; ++i)
    {
        double y1 = S[1](i, 2); double y2 = S[1](i + 1, 2);
        if (y1 >= 50.0 && y2 <= 50.0)
        {
            double x1 = S[1](i, 0); double x2 = S[1](i + 1, 0);
            x_at_50 = x1 + (x2 - x1) * (50.0 - y1) / (y2 - y1);
            break;
        }
    }

    ofstream f_tab3("wyniki_lab3_tabela_3.csv");
    f_tab3 << fixed << setprecision(5);
    f_tab3 << "v0x(0);omega(0);v0x_opt;omega_opt;x_end;x_dla_y_50;L_wyw\n";
    f_tab3 << x0_real(0) << ";" << x0_real(1) << ";" << opt_real.x(0) << ";" << opt_real.x(1) << ";" << x_end_opt << ";" << x_at_50 << ";" << solution::f_calls << "\n";
    f_tab3.close();

    ofstream f_sim("lab3_symulacja.csv");
    f_sim << fixed << setprecision(4);
    f_sim << "t;x;vx;y;vy\n";
    for (int i = 0; i < n; ++i)
    {
        if (S[1](i, 2) >= -0.5)
            f_sim << S[0](i) << ";" << S[1](i, 0) << ";" << S[1](i, 1) << ";" << S[1](i, 2) << ";" << S[1](i, 3) << "\n";
    }
    f_sim.close();
    delete[] S;
}

void lab4() {
    // 1. Wczytywanie danych
    matrix X(3, 100);
    matrix Y(1, 100);
    for (int i = 0; i < 100; ++i) X(0, i) = 1.0; 

    try {
        ifstream fileX("XData.txt");
        ifstream fileY("YData.txt");
        if (!fileX.is_open() || !fileY.is_open()) throw runtime_error("Brak plikow danych.");
        
        char dummy; double val;
        for (int row = 1; row <= 2; ++row) {
            for (int col = 0; col < 100; ++col) {
                fileX >> val; X(row, col) = val; fileX >> dummy;
            }
        }
        for (int col = 0; col < 100; ++col) {
            fileY >> val; Y(0, col) = val; fileY >> dummy;
        }
    } catch (...) { }

    // 2. Tabela 1 i 2 (Funkcja Testowa)
    ofstream f1("wyniki_lab4_tabela1.csv");
    ofstream f2("wyniki_lab4_tabela2.csv");
    
    f1 << "Dlugosc kroku,Lp.,x1(0),x2(0),Metoda najszybszego spadku,,,,,,Metoda gradientow sprzezonych,,,,,,Metoda Newtona,,,,,," << endl;
    f1 << ",,,,x1*,x2*,y*,f_calls,g_calls,MinGlob,x1*,x2*,y*,f_calls,g_calls,MinGlob,x1*,x2*,y*,f_calls,g_calls,H_calls,MinGlob" << endl;

    f2 << "Dlugosc kroku,Metoda najszybszego spadku,,,,,,Metoda gradientow sprzezonych,,,,,,Metoda Newtona,,,,,," << endl;
    f2 << ",x1*,x2*,y*,f_calls,g_calls,Liczba MinGlob,x1*,x2*,y*,f_calls,g_calls,Liczba MinGlob,x1*,x2*,y*,f_calls,g_calls,H_calls,Liczba MinGlob" << endl;

    random_device rd; mt19937 gen(rd());
    uniform_real_distribution<> dis(-2.0, 2.0);
    
    vector<matrix> start_points(100);
    for(int i=0; i<100; ++i) { start_points[i] = matrix(2, 1); start_points[i](0)=dis(gen); start_points[i](1)=dis(gen); }

    struct RunConfig { double h; string name; };
    vector<RunConfig> configs = { {0.05, "0.05"}, {0.25, "0.25"}, {-1.0, "M. zk."} };

    for (const auto& cfg : configs) {
        double sd_sum[5] = {0}; int sd_glob = 0; 
        double cg_sum[5] = {0}; int cg_glob = 0;
        double new_sum[6] = {0}; int new_glob = 0; 

        for (int i = 0; i < 100; ++i) {
            matrix x0 = start_points[i];
            
            if (i == 0) f1 << cfg.name; 
            f1 << "," << (i + 1) << "," << x0(0) << "," << x0(1) << ",";

            // --- SD ---
            solution::clear_calls(); g_calls_cnt = 0;
            solution sol_sd = SD(func_test, grad_test, x0, cfg.h, 1e-3, 10000);
            bool sd_is_glob = (abs(sol_sd.y(0)) < 1e-2); 
            f1 << sol_sd.x(0) << "," << sol_sd.x(1) << "," << sol_sd.y(0) << "," << sol_sd.f_calls << "," << g_calls_cnt << "," << (sd_is_glob ? "TAK" : "NIE") << ",";
            
            if (sd_is_glob) {
                sd_sum[0]+=sol_sd.x(0); sd_sum[1]+=sol_sd.x(1); sd_sum[2]+=sol_sd.y(0); 
                sd_sum[3]+=sol_sd.f_calls; sd_sum[4]+=g_calls_cnt; sd_glob++;
            }

            // --- CG ---
            solution::clear_calls(); g_calls_cnt = 0;
            solution sol_cg = CG(func_test, grad_test, x0, cfg.h, 1e-3, 10000);
            bool cg_is_glob = (abs(sol_cg.y(0)) < 1e-2);
            f1 << sol_cg.x(0) << "," << sol_cg.x(1) << "," << sol_cg.y(0) << "," << sol_cg.f_calls << "," << g_calls_cnt << "," << (cg_is_glob ? "TAK" : "NIE") << ",";

            if (cg_is_glob) {
                cg_sum[0]+=sol_cg.x(0); cg_sum[1]+=sol_cg.x(1); cg_sum[2]+=sol_cg.y(0); 
                cg_sum[3]+=sol_cg.f_calls; cg_sum[4]+=g_calls_cnt; cg_glob++;
            }

            // --- Newton ---
            solution::clear_calls(); g_calls_cnt = 0; H_calls_cnt = 0;
            solution sol_new = Newton(func_test, grad_test, hess_test, x0, cfg.h, 1e-6, 10000);
            bool new_is_glob = (abs(sol_new.y(0)) < 1e-2);
            f1 << sol_new.x(0) << "," << sol_new.x(1) << "," << sol_new.y(0) << "," << sol_new.f_calls << "," << g_calls_cnt << "," << H_calls_cnt << "," << (new_is_glob ? "TAK" : "NIE");

            if (new_is_glob) {
                new_sum[0]+=sol_new.x(0); new_sum[1]+=sol_new.x(1); new_sum[2]+=sol_new.y(0); 
                new_sum[3]+=sol_new.f_calls; new_sum[4]+=g_calls_cnt; new_sum[5]+=H_calls_cnt; new_glob++;
            }

            f1 << endl;
        }

        auto avg = [](double sum, int count) { return count > 0 ? sum/count : 0; };
        
        f2 << cfg.name << ",";
        f2 << avg(sd_sum[0], sd_glob) << "," << avg(sd_sum[1], sd_glob) << "," << avg(sd_sum[2], sd_glob) << "," << avg(sd_sum[3], sd_glob) << "," << avg(sd_sum[4], sd_glob) << "," << sd_glob << ",";
        f2 << avg(cg_sum[0], cg_glob) << "," << avg(cg_sum[1], cg_glob) << "," << avg(cg_sum[2], cg_glob) << "," << avg(cg_sum[3], cg_glob) << "," << avg(cg_sum[4], cg_glob) << "," << cg_glob << ",";
        f2 << avg(new_sum[0], new_glob) << "," << avg(new_sum[1], new_glob) << "," << avg(new_sum[2], new_glob) << "," << avg(new_sum[3], new_glob) << "," << avg(new_sum[4], new_glob) << "," << avg(new_sum[5], new_glob) << "," << new_glob << endl;
    }
    f1.close(); f2.close();

    // 2.5 Generowanie danych do wykresów (dla jednego punktu startowego)
    ofstream f_graphs("wyniki_lab4_wykresy.csv");
    f_graphs << "Nr iteracji,Metoda najszybszego spadku,,,,,,Metoda gradientow sprzezonych,,,,,,Metoda Newtona,,,,,," << endl;
    f_graphs << ",0.05,,0.25,,M. zk.,,0.05,,0.25,,M. zk.,,0.05,,0.25,,M. zk.," << endl;
    f_graphs << ",x1*,x2*,x1*,x2*,x1*,x2*,x1*,x2*,x1*,x2*,x1*,x2*,x1*,x2*,x1*,x2*,x1*,x2*" << endl;

    matrix x0_plot = start_points[0]; // Wybrany punkt startowy
    vector<matrix> plot_results[9]; // 9 wariantów: 3 metody * 3 kroki
    
    int plot_idx = 0;
    double steps_arr[] = {0.05, 0.25, -1.0};

    // SD
    for(double h : steps_arr) {
        g_path_trace = &plot_results[plot_idx];
        solution::clear_calls();
        SD(func_test, grad_test, x0_plot, h, 1e-6, 10000);
        g_path_trace = nullptr;
        plot_idx++;
    }
    // CG
    for(double h : steps_arr) {
        g_path_trace = &plot_results[plot_idx];
        solution::clear_calls();
        CG(func_test, grad_test, x0_plot, h, 1e-6, 10000);
        g_path_trace = nullptr;
        plot_idx++;
    }
    // Newton
    for(double h : steps_arr) {
        g_path_trace = &plot_results[plot_idx];
        solution::clear_calls();
        Newton(func_test, grad_test, hess_test, x0_plot, h, 1e-6, 10000);
        g_path_trace = nullptr;
        plot_idx++;
    }

    size_t max_iter = 0;
    for(int i=0; i<9; ++i) if(plot_results[i].size() > max_iter) max_iter = plot_results[i].size();

    for(size_t i=0; i<max_iter; ++i) {
        f_graphs << i;
        for(int j=0; j<9; ++j) {
            f_graphs << ",";
            if(i < plot_results[j].size()) {
                f_graphs << plot_results[j][i](0) << "," << plot_results[j][i](1);
            } else {
                f_graphs << ","; 
            }
        }
        f_graphs << endl;
    }
    f_graphs.close();

    // 3. Klasyfikator
    ofstream f3("wyniki_lab4_klasyfikator.csv");
    f3 << "Dlugosc kroku,Metoda gradientow sprzezonych,,,,," << endl;
    f3 << ",Theta0*,Theta1*,Theta2*,J(Theta*),P(Theta*),g_calls" << endl;

    vector<double> steps = { 0.01, 0.001, 0.0001 };
    matrix theta_start(3, 1); theta_start(0)=0; theta_start(1)=0; theta_start(2)=0;

    for (double step : steps) {
        solution::clear_calls(); g_calls_cnt = 0;
        solution sol = CG(func_real, grad_real, theta_start, step, 1e-6, 10000, X, Y);
        double acc = calculate_accuracy(sol.x, X, Y);
        
        f3 << step << "," << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << "," << sol.y(0) << "," << acc << "," << g_calls_cnt << endl;
    }
    f3.close();

    cout << "Wygenerowano pliki:\n- wyniki_lab4_tabela1.csv\n- wyniki_lab4_tabela2.csv\n- wyniki_lab4_wykresy.csv\n- wyniki_lab4_klasyfikator.csv" << endl;
}

void lab5() {}
void lab6() {}