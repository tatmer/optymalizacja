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

// Testowa funkcja celu [cite: 10]
matrix ff1T(matrix x_m, matrix ud1, matrix ud2)
{
	double x = x_m(0);
	double t = 0.1 * x;
	double y = -cos(t) * exp(-pow(t - 2.0 * M_PI, 2.0)) + 0.002 * pow(t, 2.0);
	return matrix(y);
}

// Równania różniczkowe dla problemu zbiorników [cite: 31-41]
matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	// Stałe
	double PA = 2.0;         // m^2 [cite: 31]
	double PB = 1.0;         // m^2 [cite: 32]
	double TA0 = 95.0;       // st. C [cite: 31]
	double TBin = 20.0;      // st. C [cite: 34]
	double FBin = 10.0 / 1000.0; // L/s -> m^3/s [cite: 34]
	double DB = 36.5665 / 10000.0; // cm^2 -> m^2 [cite: 35]
	double a = 0.98;         // [cite: 38]
	double b = 0.63;         // [cite: 38]
	double g = 9.81;         // [cite: 38]
	double abg_sqrt = a * b * sqrt(2.0 * g);

	// Zmienna decyzyjna (pole przekroju D_A)
	// Musi być podana w m^2
	double DA = ud2(0); 

	// Wektor stanu Y
	// Y(0) = VA (objętość w A)
	// Y(1) = VB (objętość w B)
	// Y(2) = TB (temperatura w B)
	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);

	matrix dY(3, 1);

	// Przepływy
	double F_A_out = 0;
	if (VA > 1e-6) // Zabezpieczenie przed ujemną objętością
		F_A_out = abg_sqrt * DA * sqrt(VA / PA); // [cite: 37]
	
	double F_B_out = 0;
	if (VB > 1e-6)
		F_B_out = abg_sqrt * DB * sqrt(VB / PB); // [cite: 37]

	// Równania różniczkowe
	// dV_A / dt
	dY(0) = -F_A_out;

	// dV_B / dt
	dY(1) = F_A_out + FBin - F_B_out;

	// dT_B / dt (bilans ciepła)
	if (VB > 1e-6)
	{
		// d(VB*TB)/dt = F_A_out*TA0 + FBin*TBin - F_B_out*TB
		// (dVB/dt)*TB + VB*(dTB/dt) = ...
		// VB*(dTB/dt) = F_A_out*TA0 + FBin*TBin - F_B_out*TB - (dVB/dt)*TB
		// VB*(dTB/dt) = F_A_out*TA0 + FBin*TBin - F_B_out*TB - (F_A_out + FBin - F_B_out)*TB
		// VB*(dTB/dt) = F_A_out*(TA0 - TB) + FBin*(TBin - TB)
		dY(2) = (F_A_out * (TA0 - TB) + FBin * (TBin - TB)) / VB; // [cite: 40]
	}
	else
	{
		dY(2) = 0;
	}
	
	return dY;
}

// Rzeczywista funkcja celu (zbiorniki) [cite: 42]
matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	// x(0) to D_A w cm^2 [cite: 43]
	double DA_cm2 = x(0);
	if (DA_cm2 < 0) return 1e10; // Karcimy za wartości ujemne
	matrix DA_m2 = matrix(DA_cm2 / 10000.0); // Przeliczenie na m^2

	// Warunki początkowe [cite: 31, 32]
	matrix Y0(3, 1);
	Y0(0) = 5.0; // VA0 = 5 m^3
	Y0(1) = 1.0; // VB0 = 1 m^3
	Y0(2) = 20.0; // TB0 = 20 st. C

	// Rozwiązanie ODE [cite: 43]
	matrix* Y = solve_ode(df1, 0, 1.0, 2000.0, Y0, NAN, DA_m2);

	// Y[0] - wektor czasu
	// Y[1] - macierz wyników (N x 3), gdzie kolumna 2 (indeks 2) to TB
	int n = get_len(Y[0]);
	double TB_max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
	{
		if (Y[1](i, 2) > TB_max)
			TB_max = Y[1](i, 2);
	}

	// Zwolnienie pamięci
	delete[] Y;

	// Funkcja celu: |TB_max - 50| [cite: 42]
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