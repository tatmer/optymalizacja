#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

// --- Laboratorium 1 ---

// Testowa funkcja celu 
matrix ff1T(matrix x_m, matrix ud1, matrix ud2)
{
	double x = x_m(0);
	double t = 0.1 * x;
	double y = -cos(t) * exp(-pow(t - 2.0 * M_PI, 2.0)) + 0.002 * pow(t, 2.0);
	return matrix(y);
}

// Równania różniczkowe dla problemu zbiorników 
matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	// Stałe
	double PA = 2.0;        
	double PB = 1.0;         
	double TA0 = 95.0;      
	double TBin = 20.0;     
	double FBin = 10.0 / 1000.0; 
	double DB = 36.5665 / 10000.0; 
	double a = 0.98;         
	double b = 0.63;        
	double g = 9.81;         
	double abg_sqrt = a * b * sqrt(2.0 * g);

	double DA = ud2(0); 

	
	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);

	matrix dY(3, 1);

	
	double F_A_out = 0;
	if (VA > 1e-6) // Zabezpieczenie przed ujemną objętością
		F_A_out = abg_sqrt * DA * sqrt(VA / PA); 
	
	double F_B_out = 0;
	if (VB > 1e-6)
		F_B_out = abg_sqrt * DB * sqrt(VB / PB); 

	// Równania różniczkowe
	// dV_A / dt
	dY(0) = -F_A_out;

	// dV_B / dt
	dY(1) = F_A_out + FBin - F_B_out;

	// dT_B / dt (bilans ciepła)
	if (VB > 1e-6)
	{
		dY(2) = (F_A_out * (TA0 - TB) + FBin * (TBin - TB)) / VB; 
	}
	else
	{
		dY(2) = 0;
	}
	
	return dY;
}
// Rzeczywista funkcja celu (zbiorniki) 
matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	double DA_cm2 = x(0);
	if (DA_cm2 < 0) return 1e10; 
	matrix DA_m2 = matrix(DA_cm2 / 10000.0); 
	matrix Y0(3, 1);
	Y0(0) = 5.0; 
	Y0(1) = 1.0; 
	Y0(2) = 20.0; 
	matrix* Y = solve_ode(df1, 0, 1.0, 2000.0, Y0, NAN, DA_m2);
	int n = get_len(Y[0]);
	double TB_max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
	{
		if (Y[1](i, 2) > TB_max)
			TB_max = Y[1](i, 2);
	}

	delete[] Y;

	double y = abs(TB_max - 50.0);
	return matrix(y);
}

matrix ff_lab2_T(matrix x, matrix ud1, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	double y_val = x1 * x1 + x2 * x2 - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2.0;
	return matrix(y_val);
}

matrix df_lab2_R(double t, matrix Y, matrix ud1, matrix ud2)
{
	double mr = 1.0;   // masa ramienia [kg]
	double mc = 5.0;   // masa ciężarka [kg]
	double l = 2.0;    // długość ramienia [m]
	double b = 0.25;   // współczynnik tarcia [Nms]

	// Moment bezwładności [cite: 45]
	double I = (1.0 / 3.0) * mr * pow(l, 2) + mc * pow(l, 2);

	double k1 = ud1(0);
	double k2 = ud1(1);

	double alpha_ref = M_PI; // pi rad
	double omega_ref = 0.0;

	double alpha = Y(0);
	double omega = Y(1);

	// Moment sterujący M(t) [cite: 47]
	double M_t = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

	matrix dY(2, 1);
	dY(0) = omega;
	// Równanie dynamiki: I * alpha'' + b * alpha' = M(t) -> alpha'' = (M(t) - b*omega) / I [cite: 40]
	dY(1) = (M_t - b * omega) / I;

	return dY;
}

// Funkcja celu Q
// x(0) = k1, x(1) = k2
matrix ff_lab2_R(matrix x, matrix ud1, matrix ud2)
{
	matrix Y0(2, 1); // Warunki początkowe: alpha=0, omega=0

	// Parametry symulacji
	double t0 = 0;
	double dt = 0.1;
	double tend = 100;

	// Rozwiązanie ODE. Przekazujemy x (czyli k1, k2) jako ud1 do funkcji df_lab2_R
	matrix* Y = solve_ode(df_lab2_R, t0, dt, tend, Y0, x);

	int n = get_len(Y[0]); // liczba kroków czasowych
	double Q = 0;

	double alpha_ref = M_PI;
	double omega_ref = 0.0;

	// Całkowanie metodą prostokątów funkcji podcałkowej 
	// Q = Integral( 10*(ref-a)^2 + (ref-w)^2 + M(t)^2 ) dt
	for (int i = 0; i < n; ++i)
	{
		double alpha = Y[1](i, 0);
		double omega = Y[1](i, 1);

		// Odtworzenie momentu M(t) dla danego kroku (ponieważ nie jest zwracany przez solve_ode)
		double k1 = x(0);
		double k2 = x(1);
		double M_t = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

		double component = 10 * pow(alpha_ref - alpha, 2) + pow(omega_ref - omega, 2) + pow(M_t, 2);

		Q = Q + component * dt;
	}

	// Czyszczenie pamięci po solverze
	Y[0].~matrix();
	Y[1].~matrix();

	return matrix(Q);
}

// --- Laboratorium 3 ---
matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	double a = (!isnan(ud1(0))) ? ud1(0) : 4.0;

	// Funkcja celu
	double arg = sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2));
	double y = (abs(arg) < 1e-9) ? 1.0 : sin(M_PI * arg) / (M_PI * arg);

	if (isnan(ud2(0)) || isnan(ud2(1))) return matrix(y);

	double c = ud2(0);
	int type = (int)ud2(1);
	double penalty = 0;

	// Ograniczenia
	double g1 = -x1 + 1.0;
	double g2 = -x2 + 1.0;
	double g3 = sqrt(x1 * x1 + x2 * x2) - a;

	if (type == 1) // Zewnetrzna
	{
		penalty += pow(max(0.0, g1), 2);
		penalty += pow(max(0.0, g2), 2);
		penalty += pow(max(0.0, g3), 2);
		return matrix(y + c * penalty);
	}
	else if (type == 2) // Wewnetrzna
	{
		if (g1 >= -1e-7 || g2 >= -1e-7 || g3 >= -1e-7) return matrix(1e10);
		penalty -= 1.0 / g1;
		penalty -= 1.0 / g2;
		penalty -= 1.0 / g3;
		return matrix(y + c * penalty);
	}
	return matrix(y);
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	double m = 0.6, r = 0.12, C = 0.47, rho = 1.2, g = 9.81;
	double S = M_PI * r * r;
	double omega = (!isnan(ud2(0, 0))) ? ud2(0) : ((!isnan(ud1(0, 0))) ? ud1(0) : 0.0);
	double vx = Y(1), vy = Y(3);

	// Siły oporu (przeciwne do zwrotu prędkości, stąd vx * abs(vx) w mianowniku daje znak, a minus we wzorze na dY)
	double Dx = 0.5 * C * rho * S * vx * abs(vx);
	double Dy = 0.5 * C * rho * S * vy * abs(vy);

	// Siły Magnusa
	double FMx = rho * vy * omega * M_PI * pow(r, 3);
	double FMy = rho * vx * omega * M_PI * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = vx;
	dY(1) = -(Dx + FMx) / m;
	dY(2) = vy;
	dY(3) = -(Dy + FMy) / m - g;
	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	double v0x = x(0);
	double omega = x(1);
	double c = (!isnan(ud2(0))) ? ud2(0) : 0.0;

	matrix Y0(4, 1);
	Y0(0) = 0.0; Y0(1) = v0x; Y0(2) = 100.0; Y0(3) = 0.0;
	matrix omega_mat(1, 1); omega_mat(0) = omega;

	matrix* S = solve_ode(df3, 0, 0.01, 7.0, Y0, NAN, omega_mat);
	matrix Y = S[1];
	int n = get_len(S[0]);

	double x_at_50 = 0.0;
	double x_end = 0.0;
	bool found_50 = false, found_0 = false;

	for (int i = 0; i < n - 1; ++i)
	{
		double y_curr = Y(i, 2), y_next = Y(i + 1, 2);
		double x_curr = Y(i, 0), x_next = Y(i + 1, 0);

		if (!found_50 && y_curr >= 50.0 && y_next <= 50.0) {
			x_at_50 = x_curr + (x_next - x_curr) * (50.0 - y_curr) / (y_next - y_curr);
			found_50 = true;
		}
		if (!found_0 && y_curr >= 0.0 && y_next <= 0.0) {
			x_end = x_curr + (x_next - x_curr) * (0.0 - y_curr) / (y_next - y_curr);
			found_0 = true;
		}
	}
	if (!found_0) x_end = Y(n - 1, 0);
	delete[] S;

	double penalty = 0;
	// v0x [-10, 10]
	if (v0x < -10) penalty += pow(-10 - v0x, 2);
	else if (v0x > 10) penalty += pow(v0x - 10, 2);
	// omega [-10, 10]
	if (omega < -10) penalty += pow(-10 - omega, 2);
	else if (omega > 10) penalty += pow(omega - 10, 2);
	// x na 50m w odleglosci 2m od x=5
	if (abs(x_at_50 - 5.0) > 2.0) penalty += pow(abs(x_at_50 - 5.0) - 2.0, 2);

	return matrix(-x_end + c * penalty);
}