#include <iostream>
#include <vector>

using namespace std;

class heat_task {
public:
	double T; // ������ �������, � ������� ���������� ���������������� u(x, t)
	double L; // ����� �������
	int n; // ������ ����� �� x
	int m; // ������ ����� �� t
	double initial_condition(double x); // �������, �������� ��������� �������
	double left_condition(double t); // �������, �������� ��������� ������� ��� x = 0
	double right_condition(double t); // �������, �������� ��������� ������� ��� x = L
	double f(double x, double t); // �������, �������� ������� �����������
};






void heat_equation_crank_nicolson(heat_task task, double* v)
{

}
int main()
{

}