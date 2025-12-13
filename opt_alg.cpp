#include"opt_alg.h"
#include<cmath>
#include<fstream> // Niezbędne do zapisu plików

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}



const int FIB_MAX = 2000; // Przyjmujemy maksymalną liczbę iteracji N=20 (można zmienić)
// Tablica liczb Fibonacciego F0 do F20: 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765
double F[FIB_MAX + 1];

void calculate_fib() {
	F[0] = 0;
	F[1] = 1;
	for (int i = 2; i <= FIB_MAX; ++i) {
		F[i] = F[i - 1] + F[i - 2];
	}
}

// ---------------------- METODA EKSPANSJI ----------------------
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2) // throw (string)
{
	// Zmieniona na "zmodyfikowaną metodę ekspansji" (zgodnie z instrukcją)
	try {
		double* p = new double[2] { 0, 0 };
		solution::clear_calls();

		solution X0(x0);
		X0.fit_fun(ff, ud1, ud2);

		solution X1(x0 + d);
		X1.fit_fun(ff, ud1, ud2);

		if (X1.y(0) > X0.y(0)) {
			d = -d;
			X1.x(0) = x0 + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y(0) > X0.y(0)) {
				p[0] = x0 + d;
				p[1] = x0 - d;
				return p;
			}
		}

		double a, b;
		solution X2(x0);

		a = x0;
		b = x0 + d;
		d = alpha * d;

		X2.x(0) = x0 + d;
		X2.fit_fun(ff, ud1, ud2);

		while (X2.y(0) < X1.y(0) && solution::f_calls < Nmax) {
			a = X1.x(0);
			b = X2.x(0);

			X1 = X2;
			d = alpha * d;
			X2.x(0) = X1.x(0) + d;
			X2.fit_fun(ff, ud1, ud2);
		}

		// Zwracamy przedział [a, b] posortowany rosnąco
		p[0] = min(a, b);
		p[1] = max(a, b);
		return p;

	}
	catch (string ex_info) {
		throw("double* expansion(...):\n" + ex_info);
	}
}


// ---------------------- METODA FIBONACCIEGO ----------------------
solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) // throw (string)
{
	try {
		calculate_fib(); // Upewnienie się, że tablica F jest obliczona

		solution Xopt;
		solution::clear_calls();

		double L = b - a;
		int N = 0;

		// Znajdź N takie, że F_N >= L/epsilon
		while (F[N] < L / epsilon && N <= FIB_MAX) {
			N++;
		}
		if (N > FIB_MAX) {
			// Zbyt mała dokładność lub za duży przedział
			throw string("fib(...):\nNmax osiagniete lub blad w obliczeniach F_N.");
		}

		double a_curr = a;
		double b_curr = b;
		double x1, x2, f1, f2;

		x1 = a_curr + (F[N - 2] / F[N]) * (b_curr - a_curr);
		x2 = a_curr + (F[N - 1] / F[N]) * (b_curr - a_curr);

		matrix mX1(1, 1, x1);
		matrix mX2(1, 1, x2);
		solution X1(mX1), X2(mX2); // <-- Użycie konstruktora skalara
		f1 = X1.fit_fun(ff, ud1, ud2)(0); // <-- Błąd występuje tutaj
		f2 = X2.fit_fun(ff, ud1, ud2)(0);
		for (int k = 1; k < N - 1; ++k) {
			L = b_curr - a_curr;

			if (f1 > f2) {
				a_curr = x1;
				x1 = x2;
				f1 = f2;
				x2 = a_curr + (F[N - k - 1] / F[N - k]) * (b_curr - a_curr);

				// Ostatnia iteracja - zbieżność do epsilon
				if (k == N - 2) {
					x2 = x1 + epsilon;
				}

				X2.x(0) = x2;
				f2 = X2.fit_fun(ff, ud1, ud2)(0);
			}
			else {
				b_curr = x2;
				x2 = x1;
				f2 = f1;
				x1 = a_curr + (F[N - k - 2] / F[N - k]) * (b_curr - a_curr);

				// Ostatnia iteracja
				if (k == N - 2) {
					x1 = x2 - epsilon;
				}

				X1.x(0) = x1;
				f1 = X1.fit_fun(ff, ud1, ud2)(0);
			}
		}

		// Wynik to środek ostatniego przedziału
		Xopt.x(0) = (a_curr + b_curr) / 2.0;
		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
		Xopt.ud = matrix(2, 1);
		Xopt.ud(0) = a_curr; // dolne ograniczenie
		Xopt.ud(1) = b_curr; // górne ograniczenie
		return Xopt;

	}
	catch (string ex_info) {
		throw("solution fib(...):\n" + ex_info);
	}
}


// ---------------------- METODA LAGRANGE'A (KWADRATOWA) ----------------------
solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) // throw (string)
{
	try {
		solution Xopt;
		solution::clear_calls();

		// c jest losowym punktem wewnątrz [a, b]
		double c = a + gamma * (b - a);
		if (gamma < 0 || gamma > 1) { // Lepszym podejściem jest użycie gamma jako stałej np. 0.5
			gamma = 0.5;
			c = a + gamma * (b - a);
		}

		double a_curr = a;
		double b_curr = b;
		double c_curr = c;
		double d_curr; // d to punkt d^{(i)} z pseudokodu

		solution Xa(a_curr), Xb(b_curr), Xc(c_curr);
		double fa, fb, fc;

		fa = Xa.fit_fun(ff, ud1, ud2)(0);
		fb = Xb.fit_fun(ff, ud1, ud2)(0);
		fc = Xc.fit_fun(ff, ud1, ud2)(0);

		int i = 0;
		while ((b_curr - a_curr) > epsilon && i < Nmax) {
			i++;

			// Wyznaczanie d_curr (d^{(i)} z pseudokodu)
			double licznik = fa * (pow(b_curr, 2) - pow(c_curr, 2)) + fb * (pow(c_curr, 2) - pow(a_curr, 2)) + fc * (pow(a_curr, 2) - pow(b_curr, 2));
			double mianownik = fa * (b_curr - c_curr) + fb * (c_curr - a_curr) + fc * (a_curr - b_curr);

			if (abs(mianownik) < 1e-12) {
				// Mianownik bliski zeru, funkcja prawie liniowa lub błąd
				// Zmniejszamy przedział o stały krok
				if (fa < fb) b_curr = c_curr;
				else a_curr = c_curr;
				c_curr = a_curr + (b_curr - a_curr) / 2.0;
			}
			else {
				d_curr = 0.5 * licznik / mianownik;

				// Sprawdzenie warunków zbieżności
				if (d_curr < a_curr || d_curr > b_curr) {
					// d_curr poza przedziałem [a, b], przechodzimy na prostszą metodę
					if (fa < fb) b_curr = c_curr;
					else a_curr = c_curr;
					c_curr = a_curr + (b_curr - a_curr) / 2.0;
				}
				else {
					solution Xd(d_curr);
					double fd = Xd.fit_fun(ff, ud1, ud2)(0);

					if (fd < fc) {
						if (d_curr < c_curr) {
							b_curr = c_curr;
							c_curr = d_curr;
							fc = fd;
						}
						else {
							a_curr = c_curr;
							c_curr = d_curr;
							fc = fd;
						}
					}
					else { // fd >= fc
						if (d_curr < c_curr) {
							a_curr = d_curr;
							fa = fd;
						}
						else {
							b_curr = d_curr;
							fb = fd;
						}
					}
				}
			}

			// Ponowne obliczenie fa i fb dla nowych a i b
			Xa.x(0) = a_curr; fa = Xa.fit_fun(ff, ud1, ud2)(0);
			Xb.x(0) = b_curr; fb = Xb.fit_fun(ff, ud1, ud2)(0);
		}

		Xopt.x(0) = c_curr;
		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
		Xopt.ud = matrix(2, 1);
		Xopt.ud(0) = a_curr;
		Xopt.ud(1) = b_curr;
		return Xopt;

	}
	catch (string ex_info) {
		throw("solution lag(...):\n" + ex_info);
	}
}




solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_dim(XB);
		for (int j = 0; j < n; ++j)
		{
			solution X_test = XB;
			X_test.x(j) = XB.x(j) + s;
			X_test.fit_fun(ff, ud1, ud2);
			if (X_test.y < XB.y)
			{
				XB = X_test;
			}
			else
			{
				X_test.x(j) = XB.x(j) - s;
				X_test.fit_fun(ff, ud1, ud2);
				if (X_test.y < XB.y)
				{
					XB = X_test;
				}
			}
		}
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// --- LOGIKA ZAPISU DO PLIKU ---
		bool record = false;
		// Jeśli ud1 nie jest puste i ma wartość 1, włączamy nagrywanie
		if (get_len(ud1) > 0 && !isnan(ud1(0, 0)) && ud1(0, 0) == 1.0) record = true;

		ofstream Traj;
		if (record) Traj.open("trajektoria_HJ.csv");
		// ------------------------------

		solution X(x0), XB(x0);
		X.fit_fun(ff, ud1, ud2);
		XB.fit_fun(ff, ud1, ud2);

		while (true)
		{
			// Zapis punktu bieżącego
			if (record) Traj << XB.x(0) << ";" << XB.x(1) << endl;

			XB = X;
			X = HJ_trial(ff, XB, s, ud1, ud2);

			if (X.y < XB.y)
			{
				while (true)
				{
					// Zapis punktu w pętli pattern move
					if (record) Traj << XB.x(0) << ";" << XB.x(1) << endl;

					solution X_prev = XB;
					XB = X;
					X.x = 2.0 * XB.x - X_prev.x;
					X.fit_fun(ff, ud1, ud2);
					X = HJ_trial(ff, X, s, ud1, ud2);

					if (solution::f_calls > Nmax)
					{
						X.flag = 0;
						if (record) Traj.close();
						return X;
					}
					if (X.y >= XB.y) break;
				}
				X = XB;
			}
			else
			{
				s = alpha * s;
			}

			if (solution::f_calls > Nmax)
			{
				X.flag = 0;
				break;
			}
			if (s < epsilon)
			{
				X.flag = 1;
				break;
			}
		}

		if (record) Traj.close();
		return X;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// --- LOGIKA ZAPISU DO PLIKU ---
		bool record = false;
		if (get_len(ud1) > 0 && !isnan(ud1(0, 0)) && ud1(0, 0) == 1.0) record = true;

		ofstream Traj;
		if (record) Traj.open("trajektoria_Rosen.csv");
		// ------------------------------

		int n = get_len(x0);
		matrix D = ident_mat(n);
		matrix s = s0;
		matrix lambda(n, 1);
		matrix p(n, 1);

		solution X(x0);
		X.fit_fun(ff, ud1, ud2);
		solution XB = X;

		while (true)
		{
			// Zapis punktu
			if (record) Traj << X.x(0) << ";" << X.x(1) << endl;

			matrix s_prev = s;
			XB = X;

			for (int j = 0; j < n; ++j)
			{
				matrix dj = get_col(D, j);
				solution X_test = X;
				X_test.x = X.x + s(j) * dj;
				X_test.fit_fun(ff, ud1, ud2);

				if (X_test.y < X.y)
				{
					X = X_test;
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else
				{
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}
			}

			if (solution::f_calls > Nmax)
			{
				X.flag = 0;
				break;
			}

			bool change_basis = true;
			for (int j = 0; j < n; ++j)
				if (lambda(j) == 0.0 || p(j) == 0.0)
				{
					change_basis = false;
					break;
				}

			if (change_basis)
			{
				matrix Q(n, n);
				for (int j = 0; j < n; ++j)
				{
					matrix Qj(n, 1);
					for (int k = j; k < n; ++k)
						Qj = Qj + get_col(D, k) * lambda(k);
					Q.set_col(Qj, j);
				}

				matrix D_new(n, n);
				matrix v = get_col(Q, 0);
				D_new.set_col(v / norm(v), 0);

				for (int j = 1; j < n; ++j)
				{
					v = get_col(Q, j);
					for (int k = 0; k < j; ++k)
					{
						matrix dk = get_col(D_new, k);
						v = v - (trans(get_col(Q, j)) * dk) * dk;
					}
					D_new.set_col(v / norm(v), j);
				}
				D = D_new;
				lambda = matrix(n, 1);
				p = matrix(n, 1);
				s = s0;
			}

			double max_s = 0;
			for (int j = 0; j < n; ++j)
				if (abs(s_prev(j)) > max_s) max_s = abs(s_prev(j));

			if (max_s < epsilon)
			{
				X.flag = 1;
				break;
			}
		}

		if (record) Traj.close();
		return X;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}


// ---------------------- FUNKCJA KARY ----------------------
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		double alpha = dc;
		solution X(x0);
		solution X_prev(x0);
		matrix current_ud2 = ud2;
		current_ud2(0) = c;

		while (true)
		{
			// Algorytm Neldera-Meada
			X = sym_NM(ff, X_prev.x, 0.5, 1.0, 0.5, 2.0, 0.5, epsilon, Nmax, ud1, current_ud2);

			if (solution::f_calls > Nmax)
			{
				X.flag = 0;
				return X;
			}
			if (norm(X.x - X_prev.x) < epsilon)
			{
				X.flag = 1;
				return X;
			}

			X_prev = X;
			c = alpha * c;
			current_ud2(0) = c;
		}
	}
	catch (string ex_info) { throw ("solution pen(...):\n" + ex_info); }
}

// ---------------------- SYMPLEKS NELDERA-MEADA ----------------------
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(x0);
		solution* p = new solution[n + 1];
		p[0].x = x0;
		p[0].fit_fun(ff, ud1, ud2);

		for (int i = 1; i <= n; ++i)
		{
			matrix e = matrix(n, 1);
			e(i - 1) = 1.0;
			p[i].x = p[0].x + s * e;
			p[i].fit_fun(ff, ud1, ud2);
		}

		while (true)
		{
			// Indeksy min i max
			int i_min = 0, i_max = 0;
			for (int i = 1; i <= n; ++i)
			{
				if (p[i].y < p[i_min].y) i_min = i;
				if (p[i].y > p[i_max].y) i_max = i;
			}

			// Środek ciężkości bez max
			matrix p_bar(n, 1);
			for (int i = 0; i <= n; ++i)
				if (i != i_max) p_bar = p_bar + p[i].x;
			p_bar = p_bar / (double)n;
			solution P_bar(p_bar);

			// Odbicie
			solution P_odb(P_bar.x + alpha * (P_bar.x - p[i_max].x));
			P_odb.fit_fun(ff, ud1, ud2);

			if (P_odb.y < p[i_min].y)
			{
				// Ekspansja
				solution P_e(P_bar.x + gamma * (P_odb.x - P_bar.x));
				P_e.fit_fun(ff, ud1, ud2);
				if (P_e.y < P_odb.y)
					p[i_max] = P_e;
				else
					p[i_max] = P_odb;
			}
			else
			{
				if (P_odb.y < p[i_max].y)
				{
					// Punkt odbity jest lepszy od najgorszego (akceptacja)
					p[i_max] = P_odb;
				}
				else
				{
					// Zawężenie (Contraction) - punkt odbity nie poprawił najgorszego
					solution P_z(P_bar.x + beta * (p[i_max].x - P_bar.x));
					P_z.fit_fun(ff, ud1, ud2);

					if (P_z.y < p[i_max].y)
					{
						p[i_max] = P_z;
					}
					else
					{
						// Redukcja (Shrink)
						for (int i = 0; i <= n; ++i)
						{
							if (i != i_min)
							{
								p[i].x = delta * (p[i].x + p[i_min].x);
								p[i].fit_fun(ff, ud1, ud2);
							}
						}
					}
				}
			}

			if (solution::f_calls > Nmax)
			{
				solution ret = p[i_min];
				delete[] p;
				ret.flag = 0;
				return ret;
			}

			double max_dist = 0;
			for (int i = 0; i <= n; ++i)
			{
				double d = norm(p[i_min].x - p[i].x);
				if (d > max_dist) max_dist = d;
			}
			if (max_dist < epsilon)
			{
				solution ret = p[i_min];
				delete[] p;
				ret.flag = 1;
				return ret;
			}
		}
	}
	catch (string ex_info) { throw ("solution sym_NM(...):\n" + ex_info); }
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }