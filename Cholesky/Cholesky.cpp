#include <iostream>
#include <vector>
#include <random>
#include <ctime>
using namespace std;


struct matrix
{
  int i, j, n, m;
  matrix(int _i, int _j, int _n, int _m) : i(_i), j(_j), n(_n), m(_m) {};
};

int ind(int i, int j, int N)
{
  return i * N + j;
}
void Cholesky_Decomposition(double* A, double* L, int n);
/*
  Generates a symmetrical and positive-definite matrix.
  N [in] - Size of matrix (N x N).
  Returns: Matrix of size N x N that is symmetrical and positive-definite.
*/
double* generate_matrix(int N)
{
  // Generate a random matrix M
  double* matrix = new double [N * N];
  memset(matrix, 0, (N * N) * sizeof(*matrix));
  for (int i = 0; i < N; ++i)
  {
    for (int j = i + 1; j < N; ++j)
    {
      matrix[i * N + j] = 1 + rand() % 100;
      matrix[j * N + i] = matrix[i * N + j];
    }
  }

  for (int i = 0; i < N; ++i)
  {
    double sum = 0;
    matrix[i * N + i] = 0;
    for (int j = 0; j < N; ++j)
    {
      sum += matrix[i * N + j];
    }
    matrix[i * N + i] = sum;
  }

  return matrix;
}

void pr(double* A, int N, int M)
{
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < M; ++j)
    {
      cout << A[i * M + j] << " ";
    }
    cout << endl;
  }
}

/*
  Standard Cholesky decomposition.
  matrix [in]- A matrix of the system.
  L [out] - Output matrix L.
*/
void Cholesky_Decomposition_standard(double* A, double* L, int N)
{
  for (int i = 0; i < N; ++i)
  {
    L[i * N + i] = A[i * N + i];
    for (int k = 0; k < i; ++k)
    {
      L[i * N + i] -= L[i * N + k] * L[i * N + k];
    }
    L[i * N + i] = sqrt(L[i * N + i]);

    for (int j = i + 1; j < N; ++j)
    {
      L[j * N + i] = A[j * N + i];
      for (int k = 0; k < i; ++k)
      {
        L[j * N + i] -= L[i * N + k] * L[j * N + k];
      }
      L[j * N + i] /= L[i * N + i];
    }
  }
}

/*
  Cuts matrix A into blocks according to size r.
  A [in] - A matrix of the system.
  A11 [out] - Top left block.
  A21 [out] - Bottom left block.
  A22 [out] - Bottom right block.
  r [in] - B.lock size.
*/
void get_slices(const double* A, double* A11, double* A21, double* A22, int N, int r)
{
  int len1 = 0, len2 = 0, len3 = 0;
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < r; j++)
    {
      if (i < r)
      {
        A11[len1++] = A[i * N + j];
      }
      else
      {
        A21[len2++] = A[i * N + j];
      }   
    }
    if (i >= r)
    {
      for (int j = r; j < N; j++)
      {
        A22[len3++] = A[i * N + j];
      }
    }
  }        
}

/*
  Calculates L11 block using standard Cholesky decomposition.
  A11 [in] - A11 block of the A matrix, size is r*r.
  r [in] - Block size.
  Returns: L11 block of the L matrix.
*/
double* get_L11(double* A11, int r)
{
  double* L11 = new double [r * r];
  memset(L11, 0, (r * r) * sizeof(*L11));
  Cholesky_Decomposition_standard(A11, L11, r);
  return L11;
}

/*
  Inverts a lower triangulate matrix.
  L [in] - Lower triangular matrix to invert.
  r [in] - Block size, here it means number of elements in a row.
  Returns: M - inverted matrix L.
*/
double* invertLT(const double* L, int r)
{
  double* M = new double [r * r];
  memset(M, 0, (r * r) * sizeof(*M));
  for (int i = 0; i < r; i++)
  {
    M[i * r + i] = 1.0 / (L[i * r + i]);
    for (int j = 0; j < i; j++)
    {
      for (int k = j; k < i; k++)
      {
        M[i* r + j] += L[i * r + k] * M[k * r + j];
      }
      M[i * r + j] = -M[i * r + j] / L[i * r + i];
    }
  }
  return M;
}

/*
  Transposes a matrix.
  L [in] - Matrix to transpose.
  n [in] - Number of rows in matrix.
  m [in] - Number of cols in matrix.
  Returns: M - transposed matrix L.
*/
double* transpose(const double* L, int n, int m)
{
  double* M = new double [n * m];
  memset(M, 0, n * m * sizeof(*M));

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      M[j * m + i] = L[i * m + j];
    }
  }
  return M;
}


/*
  Calculates the L21 block.
  L11 [in] - L11 block of the L matrix.
  A21 [in] - A21 block of the A matrix.
  r [in] - Block size.
  Returns: L21 block of the L matrix.
  Remarks: It's written as xA=b. So, x = b * (A^-1). We need to calculate
  an inverse matrix for our lower triangular matrix.
  Here we need to decide what's faster - inversion of one matrix, or
  transposition of two matrices. We'll check that on some tests.
*/
double* get_L21(double* L11, double* A21, int N, int r)
{
  int Ns = N - r;
  double* L21 = new double[Ns * r];
  memset(L21, 0, (Ns * r) * sizeof(*L21));
  double* transp_L11 = transpose(L11, r, r);
  double* inv_transp_L11 = invertLT(transp_L11, r);
  
  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < r; ++j)
    {
      for (int k = 0; k < r; ++k)
      {
        L21[i * r + j] += A21[i * r + k] * inv_transp_L11[k * r + j];
      }
    }
  }
  
  return L21;
}

/*
  Calculates the L21 block.
  L11 [in] - L11 block of the L matrix.
  A21 [in] - A21 block of the A matrix.
  r [in] - Block size.
  Returns: L21 block of the L matrix.

double* get_L22(double* L21, double* A22, int N, int r)
{
  double* L22;
  double* product;
  int Ns = N - r;

  double* transp_L21 = transpose(L21, Ns, r);
  memset(L22, 0, (Ns * Ns));
  memset(product, 0, (Ns * Ns));

  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < Ns; ++j)
    {
      for (int k = 0; k < r; ++k)
      {
        product[i * Ns + j] += L21[i * r + k] * transp_L21[k * Ns + j];
      }
    }
  }

  double* tilde_A22;
  memset(tilde_A22, 0, (Ns * Ns));
  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < Ns; ++j)
    {
      tilde_A22[i * Ns + j] = A22[i * Ns + j] - product[i * Ns + j];
    }
  }
}
*/

/*
  Composes L out of blocks.
*/
double* composeL(double* L, double* L11, double* L21, int r, int N, int iter)
{
  int len1 = 0, len2 = 0;
  pr(L11, r, r);
  cout << endl;
  pr(L21, N - r, r);
  cout << endl;

  cout << "r = " << r << ", N = " << N << endl;
  for (int i = iter * r; i < iter * r + r; ++i)
  {
    for (int j = iter * r; j < iter * r + r; ++j)
    {
      L[i * N + j] = L11[len1++];
    }
  }
  for (int i = r; i < N; ++i)
  {
    for (int j = 0; j < r; ++j)
    {
      L[i * N + j] = L21[len2++];
    }
  }
  return L;
}

/*
  Block Cholesky decomposition
  A [in] - A matrix of the system.
  Returns: Matrix L.
*/
void Cholesky_Decomposition(const double* A, double* L, int N)
{
  int r = 2;
  int Ns = 0;
  int numOfIterations = N / r;

  for (int i = 0; i < numOfIterations; ++i)
  {   
    matrix A11(i * r, i * r, r, r);
    matrix A21(i * )
    Ns = N - r * i;
    if (Ns == r) 
    {
      break;
    }

    get_slices(C, A11, A21, A22, Ns, r);

    cout << "Look here" << endl;
    pr(C, Ns, Ns);
    cout << endl;
    pr(A11, r, r);
    cout << endl;
    pr(A21, Ns - r, r);
    cout << endl;
    pr(A22, Ns - r, Ns - r);
    cout << "Stop looking" << endl;

    L11 = get_L11(A11, r);
    L21 = get_L21(L11, A21, Ns, r);

    composeL(L, L11, L21, r, Ns, i);
    C = A22;
    pr(L, N, N);
    cout << endl;
  }
  //L22 = new double[(Ns - r) * (Ns - r)];
  //memset(L22, 0, (Ns - r) * (Ns - r) * sizeof(*L22));
  // Надо записывать в L по ходу итераций, мб даже менять А для сохранения памяти,
  // а потом применить на последнем шаге стандартный метод Холецкого, и записать его аутпут в матрицу
  //Cholesky_Decomposition_standard(A22, L22, Ns);
  //for (int i = Ns; i < N; ++i)
  //{
  //  for (int j = Ns; j < N; ++i)
  //  {
      //L[i * N + j] = L22[i * Ns + j];
  //  }
  //}
}

/*
  Prints the original matrix and matrix L.
  matrix [in] - Original A matrix.
  L [in] - Resulting L matrix.
*/
void print_result(const double* matrix, int N, double* L)
{

  // Print original matrix A
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      cout << matrix[i * N + j] << " ";
    }
    cout << endl;
  }

  cout << endl;

  // Print matrix L.
  
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      cout << L[i * N + j] << " ";
    }
    cout << endl;
  }

  cout << endl;
}


// Driver Code
int main()
{
  srand(time(NULL));
  const int N = 6;
  const double* matrix = generate_matrix(N);
  double* L = new double [(N * N)];
  memset(L, 0, (N*N) * sizeof * L);
  time_t begin, end;
  print_result(matrix, N, L);
  
  begin = clock();
  Cholesky_Decomposition(matrix, L, N);
  end = clock();

  print_result(matrix, N, L);
  cout << "Total time: " << (end - begin) / 1000.0 << endl;

  delete[] L;
  delete[] matrix;
  return 0;
}