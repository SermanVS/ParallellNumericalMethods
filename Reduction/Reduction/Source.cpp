#include <iostream>
#include <vector>

using namespace std;

class heat_task {
public:
	double T; // момент времени, в который необходимо аппроксимировать u(x, t)
	double L; // длина стержня
	int n; // размер сетки по x
	int m; // размер сетки по t
	double initial_condition(double x); // функция, задающая начальное условие
	double left_condition(double t); // функция, задающая граничное условие при x = 0
	double right_condition(double t); // функция, задающая граничное условие при x = L
	double f(double x, double t); // функция, задающая внешнее воздействие
};






void heat_equation_crank_nicolson(heat_task task, double* v)
{

}
int main()
{

}