/*********************************************
Kod stanowi uzupe nienie materia  w do  wicze
w ramach przedmiotu metody optymalizacji.
Kod udost pniony na licencji CC BY-SA 3.0
Autor: dr in .  ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

matrix rand_range(int n, double min_val, double max_val)
{
	// Używamy rand_mat() z biblioteki i skalujemy
	matrix R = rand_mat(n, 1); // Generuje macierz Nx1 z [0, 1]
	return R * (max_val - min_val) + min_val;
}

int main()
{
	try
	{
		//lab0();
		lab2();
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
	//Funkcja testowa
	double epsilon = 1e-2;									// dok adno  
	int Nmax = 10000;										// maksymalna liczba wywo a  funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g rne ograniczenie
		a(2, 1);											// dok adne rozwi zanie optymalne
	solution opt;											// rozwi zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik w

	//Wahadlo
	Nmax = 1000;											// dok adno  
	epsilon = 1e-2;											// maksymalna liczba wywo a  funkcji celu
	lb = 0, ub = 5;											// dolne oraz g rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si y dzia aj cy na wahad o oraz czas dzia ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi zujemy r wnanie r niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie  do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie 
	Y[0].~matrix();											// usuwamy z pami ci rozwi zanie RR
	Y[1].~matrix();
}

void lab1()
{

	// Ustawienia losowania
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-100.0, 100.0); 

	// --- Zadanie 5a: Testowa funkcja celu --- 
	cout << "--- Zadanie 5a: Funkcja testowa ---" << endl;

	// Parametry
	double epsilon = 1e-5;
	int Nmax = 1000;
	int N_runs = 100; 
	double alphas[] = { 1.2, 1.5, 2.0 }; 
	double d_exp = 1.0; // Krok dla ekspansji

	ofstream f_test("lab1_test_results.csv");
	f_test << "alpha;run;x0;a_exp;b_exp;f_calls_exp;x_fib;y_fib;f_calls_fib;x_lag;y_lag;f_calls_lag\n";
	f_test << fixed << setprecision(6);

	// Pętla po współczynnikach ekspansji 
	for (double alpha : alphas)
	{
		cout << "Testuje dla alpha = " << alpha << "..." << endl;
		for (int i = 0; i < N_runs; ++i)
		{
			solution::clear_calls();

			double x0 = dis(gen); // Losowy punkt startowy
			

			// 1. Ekspansja
			double* p = expansion(ff1T, x0, d_exp, alpha, Nmax);
			double a = p[0], b = p[1];
			int f_exp = solution::f_calls;
			delete[] p;

			// 2. Fibonacci
			int f_before_fib = solution::f_calls;
			solution sol_fib = fib(ff1T, a, b, epsilon);
			int f_fib = solution::f_calls - f_before_fib;

			// 3. Lagrange
			int f_before_lag = solution::f_calls;
			solution sol_lag = lag(ff1T, a, b, epsilon, epsilon / 10.0, Nmax);
			int f_lag = solution::f_calls - f_before_lag;

			// Zapis do pliku
			f_test << alpha << ";" << i << ";" << x0 << ";"
				<< a << ";" << b << ";" << f_exp << ";"
				<< sol_fib.x(0) << ";" << sol_fib.y(0) << ";" << f_fib << ";"
				<< sol_lag.x(0) << ";" << sol_lag.y(0) << ";" << f_lag << "\n";
		}
	}
	f_test.close();
	cout << "Wyniki testow (z ekspansja) zapisano do lab1_test_results.csv" << endl;

	// Zadanie 5a: Bez wstępnego zawężania 
	cout << "\nTestuje bez wstepnego zawezania (przedzial [-100, 100])..." << endl;
	double a_full = -100.0, b_full = 100.0;

	solution::clear_calls();
	solution sol_fib_full = fib(ff1T, a_full, b_full, epsilon);
	cout << "Fibonacci (pelny przedzial):\n" << sol_fib_full << endl;

	solution::clear_calls();
	solution sol_lag_full = lag(ff1T, a_full, b_full, epsilon, epsilon / 10.0, Nmax);
	cout << "Lagrange (pelny przedzial):\n" << sol_lag_full << endl;
	cout << "Wykresy zbieznosci zapisano do fib_intervals.csv i lag_intervals.csv" << endl; 

	// --- Zadanie 5b: Problem rzeczywisty --- 
	cout << "\n--- Zadanie 5b: Problem rzeczywisty (zbiorniki) ---" << endl;

	// Sprawdzenie poprawności implementacji modelu
	solution::clear_calls();
	solution check_sol(50.0); 
	check_sol.fit_fun(ff1R);
	
	cout << "Sprawdzenie (DA=50): y = " << check_sol.y(0)
		<< " (oczekiwano ok. 12.5, jesli T_max=62.5)" << endl;

	// Parametry optymalizacj
	double a_real = 1.0;  
	double b_real = 100.0; 
	double eps_real = 1e-2;
	int Nmax_real = 100;

	// Optymalizacja 
	solution::clear_calls();
	solution sol_fib_real = fib(ff1R, a_real, b_real, eps_real);
	cout << "\nFibonacci (problem rzeczywisty):\n" << sol_fib_real << endl;

	solution::clear_calls();
	solution sol_lag_real = lag(ff1R, a_real, b_real, eps_real, eps_real / 10.0, Nmax_real);
	cout << "Lagrange (problem rzeczywisty):\n" << sol_lag_real << endl;

	// --- Symulacja dla optymalnego DA ---
	cout << "\nGenerowanie symulacji porownawczej dla obu metod..." << endl;

	// Warunki początkowe (takie same dla obu)
	matrix Y0(3, 1);
	Y0(0) = 5.0; // VA0
	Y0(1) = 1.0; // VB0
	Y0(2) = 20.0; // TB0

	// 1. Uruchom symulację dla wyniku Fibonacciego
	matrix DA_fib_m2 = matrix(sol_fib_real.x(0) / 10000.0);
	matrix* Y_fib = solve_ode(df1, 0, 1.0, 2000.0, Y0, NAN, DA_fib_m2);

	// 2. Uruchom symulację dla wyniku Lagrange'a
	matrix DA_lag_m2 = matrix(sol_lag_real.x(0) / 10000.0);
	matrix* Y_lag = solve_ode(df1, 0, 1.0, 2000.0, Y0, NAN, DA_lag_m2);

	// 3. Połącz wyniki w jedną macierz

	// Pobierz kolumnę czasu (jest taka sama dla obu)
	matrix T = Y_fib[0];

	// Pobierz poszczególne kolumny stanu dla obu symulacji
	matrix VA_fib = get_col(Y_fib[1], 0);
	matrix VB_fib = get_col(Y_fib[1], 1);
	matrix TB_fib = get_col(Y_fib[1], 2);

	matrix VA_lag = get_col(Y_lag[1], 0);
	matrix VB_lag = get_col(Y_lag[1], 1);
	matrix TB_lag = get_col(Y_lag[1], 2);

	// Zbuduj macierz wyjściową zgodnie z szablonem Excela
	matrix C = hcat(T, VA_fib);
	C = hcat(C, VA_lag);
	C = hcat(C, VB_fib);
	C = hcat(C, VB_lag);
	C = hcat(C, TB_fib);
	C = hcat(C, TB_lag);

	// 4. Zapisz do pliku CSV
	ofstream f_sim("lab1_symulacja_opt.csv");
	// Nagłówek pasujący do szablonu
	f_sim << "t;VA_Fibonacci;VA_Lagrange;VB_Fibonacci;VB_Lagrange;TB_Fibonacci;TB_Lagrange\n";

	// Użyj operatora << do zapisu całej macierzy (dodałem formatowanie dla czytelności)
	f_sim << fixed << setprecision(6) << C;
	f_sim.close();

	// 5. Posprzątaj pamięć
	delete[] Y_fib;
	delete[] Y_lag;

	cout << "Symulacje porownawcze zapisano do lab1_symulacja_opt.csv" << endl;
}

void lab2() {
	ofstream Sout("wyniki_lab2.csv");
	if (!Sout.is_open())
	{
		cerr << "Nie mozna otworzyc pliku 'wyniki_lab2.csv' do zapisu!" << endl;
		return;
	}

	cout << "Rozpoczynam optymalizacje... Wyniki beda zapisane do 'wyniki_lab2.csv'" << endl;

	Sout << fixed << setprecision(6); // Ustawienie formatowania liczb w pliku
	cout << fixed << setprecision(6); // Ustawienie formatowania liczb na konsoli

	int N_runs = 100; //
	int Nmax = 20000;
	double epsilon = 1e-3;
	matrix lb(2, 1, -1.0); //
	matrix ub(2, 1, 1.0); //

	// --- Generowanie 100 tych samych punktów startowych ---
	std::vector<matrix> start_points;
	for (int i = 0; i < N_runs; ++i)
	{
		matrix x0 = rand_mat(2, 1); // [0, 1]
		x0(0) = (ub(0) - lb(0)) * x0(0) + lb(0); // Skalowanie do [-1, 1]
		x0(1) = (ub(1) - lb(1)) * x0(1) + lb(1); // Skalowanie do [-1, 1]
		start_points.push_back(x0);
	}
	cout << "Wygenerowano " << start_points.size() << " losowych punktow startowych." << endl;
	// ---------------------------------------------------------


	// Parametry algorytmów
	double s_hj_all[] = { 1.0, 0.1, 0.01 }; //
	double alpha_hj = 0.5; //

	matrix s0_rosen_all[] = { matrix(2, 1, 0.5), matrix(2, 1, 0.1), matrix(2, 1, 0.05) }; //
	double alpha_rosen = 2.0; //
	double beta_rosen = 0.5; //

	// --- Metoda Hooke'a-Jeevesa ---
	Sout << "==========================================\n";
	Sout << "Metoda Hooke'a-Jeevesa\n";
	Sout << "==========================================\n";

	for (double s_hj : s_hj_all)
	{
		Sout << "\n--- Krok startowy s = " << s_hj << " ---\n";
		Sout << "--- Tabela 1 (dane surowe) ---\n";
		Sout << "x1_start;x2_start;x1_opt;x2_opt;y_opt;f_calls;flag\n";

		int global_found = 0;
		double f_calls_sum = 0;

		cout << "Uruchamiam HJ z krokiem s = " << s_hj << "..." << endl;

		// Iteracja po wczesniej zdefiniowanych punktach startowych
		for (const auto& x0 : start_points)
		{
			solution::clear_calls();

			solution opt = HJ(ff_lab2_T, x0, s_hj, alpha_hj, epsilon, Nmax);

			// Zapis wyników do pliku (Tabela 1)
			Sout << x0(0) << ";" << x0(1) << ";"
				<< opt.x(0) << ";" << opt.x(1) << ";"
				<< opt.y << ";" << solution::f_calls << ";" << opt.flag << "\n";

			// Sprawdzenie, czy znaleziono minimum globalne (w x=[0,0], y=0)
			if (opt.flag == 1 && norm(opt.x) < 0.1 && opt.y < 0.01)
			{
				global_found++;
				f_calls_sum += solution::f_calls;
			}
		}

		// Zapis podsumowania do pliku (Tabela 2)
		Sout << "\n--- Tabela 2 (podsumowanie) dla s = " << s_hj << " ---\n";
		Sout << "Liczba znalezionych minimow globalnych:;" << global_found << "/" << N_runs << "\n";
		if (global_found > 0)
		{
			Sout << "Srednia liczba wywolan f. celu (dla globalnych):;" << f_calls_sum / global_found << "\n";
		}
		Sout << "------------------------------------------\n";

		cout << "Zakonczono HJ dla s = " << s_hj << ". Znaleziono " << global_found << " globalnych minimow." << endl;
	}


	// --- Metoda Rosenbrocka ---
	Sout << "\n\n==========================================\n";
	Sout << "Metoda Rosenbrocka\n";
	Sout << "==========================================\n";

	for (const auto& s0_rosen : s0_rosen_all)
	{
		double s_val = s0_rosen(0); // Zakładamy, że wszystkie kroki są równe
		Sout << "\n--- Krok startowy s = " << s_val << " ---\n";
		Sout << "--- Tabela 1 (dane surowe) ---\n";
		Sout << "x1_start;x2_start;x1_opt;x2_opt;y_opt;f_calls;flag\n";

		int global_found = 0;
		double f_calls_sum = 0;

		cout << "Uruchamiam Rosenbrocka z krokiem s = " << s_val << "..." << endl;

		// Iteracja po wczesniej zdefiniowanych punktach startowych
		for (const auto& x0 : start_points)
		{
			solution::clear_calls();

			solution opt = Rosen(ff_lab2_T, x0, s0_rosen, alpha_rosen, beta_rosen, epsilon, Nmax);

			// Zapis wyników do pliku (Tabela 1)
			Sout << x0(0) << ";" << x0(1) << ";"
				<< opt.x(0) << ";" << opt.x(1) << ";"
				<< opt.y << ";" << solution::f_calls << ";" << opt.flag << "\n";

			// Sprawdzenie, czy znaleziono minimum globalne (w x=[0,0], y=0)
			if (opt.flag == 1 && norm(opt.x) < 0.1 && opt.y < 0.01)
			{
				global_found++;
				f_calls_sum += solution::f_calls;
			}
		}

		// Zapis podsumowania do pliku (Tabela 2)
		Sout << "\n--- Tabela 2 (podsumowanie) dla s = " << s_val << " ---\n";
		Sout << "Liczba znalezionych minimow globalnych:;" << global_found << "/" << N_runs << "\n";
		if (global_found > 0)
		{
			Sout << "Srednia liczba wywolan f. celu (dla globalnych):;" << f_calls_sum / global_found << "\n";
		}
		Sout << "------------------------------------------\n";

		cout << "Zakonczono Rosenbrocka dla s = " << s_val << ". Znaleziono " << global_found << " globalnych minimow." << endl;
	}

	// Zamknięcie pliku
	Sout.close();
	cout << "\nZakonczono. Zapisano wszystkie wyniki do pliku 'wyniki_lab2.csv'." << endl;
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
