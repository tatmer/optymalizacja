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

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
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

// Pozostałe puste funkcje
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }
solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2) { return solution(); }