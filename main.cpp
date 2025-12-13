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
#include<iomanip>
#include<fstream>
#include<vector>
#include<cmath> 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
		lab3(); // Uruchomienie wszystkich zadań z Lab 2
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	// (Kod lab0 bez zmian...)
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
	// Ustawienie generatora liczb losowych
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> distr_x0(-100.0, 100.0);
	std::uniform_real_distribution<> distr_da_cm2(1.0, 100.0);

	double epsilon = 1e-4; // Dokładność
	int Nmax = 1000;      // Maksymalna liczba wywołań funkcji celu (dla Lagrange'a)

	// ------------------------------------------------------------------------------------------------
	// 5. a) FUNKCJA TESTOWA - 100 OPTYMALIZACJI Z EKSPANSJĄ
	// ------------------------------------------------------------------------------------------------
	cout << "5. a) Optymalizacja funkcji testowej..." << endl;

	// Współczynniki ekspansji do przetestowania
	vector<double> alphas = { 1.2, 1.5, 2.0 };

	// Nagłówek pliku CSV dla Tabeli 1
	ofstream file_test("wyniki_testowe.csv");
	file_test << "Alpha;x0;A_fib;B_fib;Xopt_fib;Yopt_fib;Iter_fib;A_lag;B_lag;Xopt_lag;Yopt_lag;Iter_lag" << endl;

	for (double alpha : alphas) {
		cout << "  -> Przetwarzanie dla alpha = " << alpha << endl;
		for (int i = 0; i < 100; ++i) {

			// Losowy punkt startowy dla funkcji testowej
			double x0 = distr_x0(gen);

			// 1. Wstępna ekspansja
			double d0 = 1.0; // Początkowy krok
			double* interval = expansion(ff1T, x0, d0, alpha, Nmax);
			double a = interval[0];
			double b = interval[1];
			delete[] interval;

			solution::clear_calls();

			// 2. Metoda Fibonacciego (z zawężeniem)
			solution opt_fib = fib(ff1T, a, b, epsilon, NAN, NAN);
			int iter_fib = solution::f_calls;

			// 3. Metoda Lagrange'a (z zawężeniem)
			solution::clear_calls();
			solution opt_lag = lag(ff1T, a, b, epsilon, 0.5, Nmax, NAN, NAN);
			int iter_lag = solution::f_calls;

			// Zapis do CSV (Tabela 1)
			file_test << alpha << ";"
				<< x0 << ";"
				<< opt_fib.ud(0) << ";" << opt_fib.ud(1) << ";" << opt_fib.x(0) << ";" << opt_fib.y(0) << ";" << iter_fib << ";"
				<< opt_lag.ud(0) << ";" << opt_lag.ud(1) << ";" << opt_lag.x(0) << ";" << opt_lag.y(0) << ";" << iter_lag << endl;
		}
	}
	file_test.close();

	// ------------------------------------------------------------------------------------------------
	// 5. b) PROBLEM RZECZYWISTY - OPTYMALIZACJA ZBIORNIKÓW
	// ------------------------------------------------------------------------------------------------
	cout << "\n5. b) Optymalizacja problemu rzeczywistego (zbiorniki)..." << endl;

	// Cel: TB_max = 50 st. C
	matrix TB_target(1, 1, 50.0);

	// Losowy punkt startowy DA [cm^2]
	double DA0_cm2 = distr_da_cm2(gen);
	double a_real = 1.0; // [cm^2]
	double b_real = 100.0; // [cm^2]

	// 1. Optymalizacja Metodą Fibonacciego
	solution::clear_calls();
	solution opt_fib_real = fib(ff1R, a_real, b_real, epsilon, TB_target, NAN);
	cout << "  -> Fibonacciego D_A: " << opt_fib_real.x(0) << " cm^2, Blad: " << opt_fib_real.y(0) << ", Iter: " << solution::f_calls << endl;

	// 2. Optymalizacja Metodą Lagrange'a
	solution::clear_calls();
	solution opt_lag_real = lag(ff1R, a_real, b_real, epsilon, 0.5, Nmax, TB_target, NAN);
	cout << "  -> Lagrange'a D_A: " << opt_lag_real.x(0) << " cm^2, Blad: " << opt_lag_real.y(0) << ", Iter: " << solution::f_calls << endl;

	// Wybór najlepszego wyniku (z mniejszym błędem)
	solution opt_best = (opt_fib_real.y(0) < opt_lag_real.y(0)) ? opt_fib_real : opt_lag_real;

	// 3. Pełna symulacja dla optymalnego D_A
	cout << "\n  -> Przeprowadzam symulacje dla optymalnego D_A = " << opt_best.x(0) << " cm^2" << endl;

	// Wywołanie symulacji (jak w ff1R)
	double DA_best_m2 = opt_best.x(0) / 10000.0;
	matrix Y0(3, 1);
	Y0(0) = 5.0; Y0(1) = 1.0; Y0(2) = 20.0;
	matrix ud2_solver(1, 1, DA_best_m2);

	matrix* S = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, NAN, ud2_solver);
	int N = get_len(S[0]);

	// Zapis symulacji do pliku CSV (Arkusz "Symulacja")
	ofstream file_sim("symulacja_zbiorniki.csv");
	file_sim << "Czas;VA;VB;TB" << endl;

	double max_TB = 0.0;
	for (int i = 0; i < N; ++i) {
		file_sim << S[0](i) << ";" << S[1](0, i) << ";" << S[1](1, i) << ";" << S[1](2, i) << endl;
		if (S[1](2, i) > max_TB) {
			max_TB = S[1](2, i);
		}
	}
	file_sim.close();
	delete[] S;

	cout << "  -> Maksymalna osiagnieta temperatura TB: " << max_TB << " st. C" << endl;
	cout << "\nWyniki zapisano do plikow: wyniki_testowe.csv i symulacja_zbiorniki.csv" << endl;
}

void lab2()
{
	// ========================================================================
	// CZĘŚĆ 1: Zadanie 5a - Statystyka (Tabela 1 i 2)
	// ========================================================================

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

	// Generowanie 100 stałych punktów startowych dla rzetelnego porównania
	std::vector<matrix> start_points;
	for (int i = 0; i < N_runs; ++i)
	{
		matrix x0 = rand_mat(2, 1);
		x0(0) = (ub(0) - lb(0)) * x0(0) + lb(0); // Skalowanie do [-1, 1]
		x0(1) = (ub(1) - lb(1)) * x0(1) + lb(1); // Skalowanie do [-1, 1]
		start_points.push_back(x0);
	}

	double s_hj_all[] = { 1.0, 0.1, 0.01 };
	double alpha_hj = 0.5;

	matrix s0_rosen_all[] = { matrix(2, 1, 0.5), matrix(2, 1, 0.1), matrix(2, 1, 0.05) };
	double alpha_rosen = 2.0;
	double beta_rosen = 0.5;

	// --- Metoda Hooke'a-Jeevesa (Statystyka) ---
	for (double s_hj : s_hj_all)
	{
		Sout << "\nMetoda Hooke'a-Jeevesa, s = " << s_hj << "\n";
		Sout << "x1_start;x2_start;x1_opt;x2_opt;y_opt;f_calls;flag\n";
		int global_found = 0;
		for (const auto& x0 : start_points)
		{
			solution::clear_calls();
			solution opt = HJ(ff_lab2_T, x0, s_hj, alpha_hj, epsilon, Nmax);

			Sout << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";"
				<< opt.y << ";" << solution::f_calls << ";" << opt.flag << "\n";

			// Kryterium sukcesu: blisko (0,0) i wartość bliska 0
			if (opt.flag == 1 && norm(opt.x) < 0.1 && opt.y < 0.01) global_found++;
		}
		cout << "HJ (s=" << s_hj << "): Globalne = " << global_found << "/" << N_runs << endl;
	}

	// --- Metoda Rosenbrocka (Statystyka) ---
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

			Sout << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";"
				<< opt.y << ";" << solution::f_calls << ";" << opt.flag << "\n";

			if (opt.flag == 1 && norm(opt.x) < 0.1 && opt.y < 0.01) global_found++;
		}
		cout << "Rosenbrock (s=" << s_val << "): Globalne = " << global_found << "/" << N_runs << endl;
	}
	Sout.close();


	// ========================================================================
	// CZĘŚĆ 2: Zadanie 5b - Problem Rzeczywisty (Robot)
	// ========================================================================

	ofstream SoutR("wyniki_lab2_partB.csv");
	if (!SoutR.is_open()) { cerr << "Blad zapisu partB" << endl; return; }

	cout << "\n--- Zadanie 5b: Problem Rzeczywisty ---" << endl;
	SoutR << fixed << setprecision(6);

	// Punkt startowy k1=10, k2=10 (środek przedziału [0,20])
	matrix x0_robot(2, 1);
	x0_robot(0) = 10.0;
	x0_robot(1) = 10.0;

	double s_robot = 0.5; // Krok startowy dla robota
	matrix s0_robot(2, 1, s_robot);

	// --- Optymalizacja HJ ---
	solution::clear_calls();
	solution opt_HJ = HJ(ff_lab2_R, x0_robot, s_robot, alpha_hj, epsilon, Nmax);
	int hj_f_calls = solution::f_calls; // Zapamiętujemy liczbę wywołań
	cout << "HJ zakonczone. y = " << opt_HJ.y << endl;

	// --- Optymalizacja Rosenbrocka ---
	solution::clear_calls();
	solution opt_Ros = Rosen(ff_lab2_R, x0_robot, s0_robot, alpha_rosen, beta_rosen, epsilon, Nmax);
	cout << "Rosenbrock zakonczony. y = " << opt_Ros.y << endl;

	// Zapis do Tabeli 3
	SoutR << "Metoda;k1_start;k2_start;k1_opt;k2_opt;Q_min;f_calls\n";
	SoutR << "HJ;" << x0_robot(0) << ";" << x0_robot(1) << ";"
		<< opt_HJ.x(0) << ";" << opt_HJ.x(1) << ";" << opt_HJ.y << ";" << hj_f_calls << "\n";
	SoutR << "Rosenbrock;" << x0_robot(0) << ";" << x0_robot(1) << ";"
		<< opt_Ros.x(0) << ";" << opt_Ros.x(1) << ";" << opt_Ros.y << ";" << solution::f_calls << "\n";
	SoutR.close();

	// --- Symulacja dla najlepszego wyniku ---
	// Wybieramy metodę, która dała mniejszy błąd Q
	solution best_opt = (opt_HJ.y < opt_Ros.y) ? opt_HJ : opt_Ros;
	cout << "Generowanie symulacji dla lepszego wyniku (Metoda: "
		<< ((opt_HJ.y < opt_Ros.y) ? "HJ" : "Rosenbrock") << ")..." << endl;

	matrix Y0(2, 1); // Warunki początkowe
	matrix* Y_sim = solve_ode(df_lab2_R, 0, 0.1, 100, Y0, best_opt.x);

	ofstream SimOut("symulacja_lab2.csv");
	if (SimOut.is_open())
	{
		SimOut << fixed << setprecision(6);
		// Zapis: czas; pozycja; prędkość
		int n_steps = get_len(Y_sim[0]);
		for (int i = 0; i < n_steps; ++i)
		{
			SimOut << Y_sim[0](i, 0) << ";" << Y_sim[1](i, 0) << ";" << Y_sim[1](i, 1) << "\n";
		}
		SimOut.close();
		cout << "Dane symulacji zapisane w 'symulacja_lab2.csv'. (Wklej do arkusza Symulacja)" << endl;
	}
	Y_sim[0].~matrix(); Y_sim[1].~matrix();


	// ========================================================================
	// CZĘŚĆ 3: GENEROWANIE TRAJEKTORII (Dla arkusza "Wykres" - wykres poziomic)
	// ========================================================================
	cout << "\n--- Generowanie trajektorii dla wybranego przypadku (Zad 5a Wykres) ---" << endl;

	// Wybieramy pierwszy punkt startowy z naszej listy (dla powtarzalności)
	matrix x0_traj = start_points[0];

	// Ustawiamy flagę zapisu (ud1 = 1.0) - to uruchomi kod zapisu w opt_alg.cpp
	matrix flag_record(1, 1, 1.0);

	// 1. Trajektoria HJ (najczęściej wymagana do wykresu poziomic)
	solution::clear_calls();
	// Używamy przykładowego kroku 0.5
	HJ(ff_lab2_T, x0_traj, 0.5, alpha_hj, epsilon, Nmax, flag_record);
	cout << "Utworzono plik 'trajektoria_HJ.csv' z punktami posrednimi." << endl;
	cout << "Skopiuj zawartosc tego pliku na wykres poziomic w Excelu." << endl;
}

void lab3()
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis_x(-10.0, 10.0);
	uniform_real_distribution<> dis_y(-10.0, 10.0);

	// Parametry Zadania 5a
	double a_vals[] = { 4.0, 4.4934, 5.0 };
	int N_runs = 100;
	int Nmax = 5000;
	double epsilon = 1e-3;

	// Otwarcie plików dla Tabel 1 i 2
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
			// Losowanie punktu startowego dopuszczalnego
			do {
				x0(0) = dis_x(gen);
				x0(1) = dis_y(gen);
			} while (x0(0) <= 1.001 || x0(1) <= 1.001 || norm(x0) >= a - 0.001); // Ściśle wewnątrz

			matrix x0_base = x0;
			matrix ud1(1, 1); ud1(0) = a;
			matrix clean_ud(2, 1); clean_ud(0) = 0; clean_ud(1) = 0; // Do obliczenia czystej wartości f(x)

			// --- Metoda Zewnętrzna ---
			solution::clear_calls();
			// ud2 przechowuje: (0) współczynnik kary c, (1) typ kary (1-zew, 2-wew)
			matrix ud2_zew(2, 1); ud2_zew(0) = 1.0; ud2_zew(1) = 1.0;
			// pen: x0, c=1.0, dc=2.0 (alpha), epsilon, Nmax
			solution opt_zew = pen(ff3T, x0_base, 1.0, 2.0, epsilon, Nmax, ud1, ud2_zew);

			double y_zew = ff3T(opt_zew.x, ud1, clean_ud)(0);
			double r_zew = norm(opt_zew.x);
			fout_zew << a << ";" << x0_base(0) << ";" << x0_base(1) << ";" << opt_zew.x(0) << ";" << opt_zew.x(1) << ";" << y_zew << ";" << solution::f_calls << ";" << opt_zew.flag << ";" << r_zew << "\n";

			// Statystyka (tylko dla sukcesów flag=1 lub blisko optimum)
			if (opt_zew.flag == 1) { zew_s_x1 += opt_zew.x(0); zew_s_x2 += opt_zew.x(1); zew_s_y += y_zew; zew_s_calls += solution::f_calls; zew_s_r += r_zew; zew_success++; }

			// --- Metoda Wewnętrzna ---
			solution::clear_calls();
			matrix ud2_wew(2, 1); ud2_wew(0) = 1.0; ud2_wew(1) = 2.0;
			// pen: x0, c=1.0, dc=0.5 (alpha < 1 dla wew), epsilon, Nmax
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

	// ==========================================================
	// ZADANIE 5B: PROBLEM RZECZYWISTY
	// ==========================================================
	cout << "\n--- Zadanie 5b: Problem rzeczywisty ---" << endl;

	// 1. Weryfikacja modelu (v=5, w=10)
	{
		matrix x_ver(2, 1); x_ver(0) = 5.0; x_ver(1) = 10.0;
		matrix ud_null(1, 1); ud_null(0) = 0.0; // c=0 (bez kary)
		double ver = -ff3R(x_ver, NAN, ud_null)(0);
		cout << "[TEST] v=5, w=10 -> x_end = " << ver << " m (Oczekiwano: ~41.41)" << endl;
	}

	// 2. Optymalizacja
	// Punkt startowy
	matrix x0_real(2, 1);
	x0_real(0) = -5.0;
	x0_real(1) = 6.0;

	cout << "Start optymalizacji: v0x=" << x0_real(0) << ", omega=" << x0_real(1) << endl;

	solution::clear_calls();
	// Zewnętrzna funkcja kary: ud2(0) = c.
	matrix ud2_real(1, 1); ud2_real(0) = 1.0;

	// pen: c=1.0, dc=2.0 (alpha)
	solution opt_real = pen(ff3R, x0_real, 1.0, 2.0, 1e-3, 5000, NAN, ud2_real);

	// Wynik optymalny bez kary
	matrix final_ud(1, 1); final_ud(0) = 0.0;
	double x_end_opt = -ff3R(opt_real.x, NAN, final_ud)(0);

	// 3. Symulacja dla wykresu i tabeli
	matrix Y0(4, 1);
	Y0(0) = 0.0; Y0(1) = opt_real.x(0); Y0(2) = 100.0; Y0(3) = 0.0;
	matrix omega_opt(1, 1); omega_opt(0) = opt_real.x(1);

	matrix* S = solve_ode(df3, 0, 0.01, 7.0, Y0, NAN, omega_opt);

	double x_at_50 = 0.0;
	int n = get_len(S[0]);
	for (int i = 0; i < n - 1; ++i)
	{
		double y1 = S[1](i, 2); double y2 = S[1](i + 1, 2);
		if (y1 >= 50.0 && y2 <= 50.0) // Przejście przez y=50
		{
			double x1 = S[1](i, 0); double x2 = S[1](i + 1, 0);
			x_at_50 = x1 + (x2 - x1) * (50.0 - y1) / (y2 - y1);
			break;
		}
	}

	// Tabela 3
	ofstream f_tab3("wyniki_lab3_tabela_3.csv");
	f_tab3 << fixed << setprecision(5);
	f_tab3 << "v0x(0);omega(0);v0x_opt;omega_opt;x_end;x_dla_y_50;L_wyw\n";
	f_tab3 << x0_real(0) << ";" << x0_real(1) << ";"
		<< opt_real.x(0) << ";" << opt_real.x(1) << ";"
		<< x_end_opt << ";" << x_at_50 << ";"
		<< solution::f_calls << "\n";
	f_tab3.close();

	// Wykres
	ofstream f_sim("lab3_symulacja.csv");
	f_sim << fixed << setprecision(4);
	f_sim << "t;x;vx;y;vy\n";
	for (int i = 0; i < n; ++i)
	{
		if (S[1](i, 2) >= -0.5) // Zapisz do momentu uderzenia w ziemię
			f_sim << S[0](i) << ";" << S[1](i, 0) << ";" << S[1](i, 1) << ";" << S[1](i, 2) << ";" << S[1](i, 3) << "\n";
	}
	f_sim.close();
	delete[] S;

	cout << "Zakonczono. Wszystkie pliki wygenerowane." << endl;
}

void lab4() {}
void lab5() {}
void lab6() {}