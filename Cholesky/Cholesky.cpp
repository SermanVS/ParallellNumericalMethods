#include <iostream>
#include <vector>
#include <random>
#include <ctime>
using namespace std;

void Cholesky_Decomposition(double* A, double* L, int n);

/*
  Transposes a matrix.
  L [in] - Matrix to transpose.
  n [in] - Number of rows in matrix.
  m [in] - Number of cols in matrix.
  Returns: M - transposed matrix L.
*/
double* transpose(double* L, int n, int m)
{
  double* M = new double[m * n];
  memset(M, 0, n * m * sizeof(*M));

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      M[j * n + i] = L[i * m + j];
    }
  }
  return M;
}

/*
  Calculates product of two matrices.
  A [in] - First matrix.
  B [in] - Second matrix.
  product [out] - Product of A and B.
  Ns [in] - Number of rows in resulting matrix.
  r [in] - Number of cols in resulting matrix.
*/
void prod(double* A, double* B, double* product, int Ns, int r)
{
  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < Ns; ++j)
    {
      for (int k = 0; k < r; ++k)
      {
        product[i * Ns + j] += A[i * r + k] * B[k * Ns + j];
      }
    }
  }
}

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

/*
  Prints matrix A.
  A [in] - Matrix to print.
  N [in] - Rows in A.
  M [in] - Cols in A.
*/
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
  N [in] - Dimension of square matrix.
*/
void Cholesky_Decomposition_standard(double* A, double* L, int N)
{
  memset(L, 0, N * N * sizeof(*L));
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
  Calculates L11 block using standard Cholesky decomposition.
  A11 [in] - A11 block of the A matrix, size is r*r.
  r [in] - Block size.
  Returns: L11 block of the L matrix.
*/
double* get_L11(double* A11, int r)
{
  double* L11 = new double [r * r];
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
  Calculates the L21 block from formula L21 * L11^T = A21.
  L11 [in] - L11 block of the L matrix.
  A21 [in] - A21 block of the A matrix.
  Ns [in] - Number of rows in L21.
  r [in] - Block size.
  Returns: L21 block of the L matrix.
  Remarks: It's written as xA=b. So, x = b * (A^-1). We need to calculate
  an inverse matrix for our lower triangular matrix.
  Here we need to decide what's faster - inversion of one matrix, or
  transposition of two matrices. We'll check that on some tests.
*/
double* get_L21(double* L11, double* A21, int Ns, int r)
{
  double* L21 = new double[Ns * r];
  memset(L21, 0, (Ns * r) * sizeof(*L21));
  double* inv_L11 = invertLT(L11, r);
  double* inv_transp_L11 = transpose(inv_L11, r, r);
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
  delete[] inv_L11;
  delete[] inv_transp_L11;
  return L21;
}

/*
  Cuts matrix A into blocks according to size r.
  A [in] - A matrix of the system.
  A11 [out] - Top left block.
  A21 [out] - Bottom left block.
  A22 [out] - Bottom right block.
  N [in] - Dimension of square matrix A.
  r [in] - B.lock size.
  i0 [in] - Starting row to count from.
  j0 [in] - Starting col to count from.
  Remarks:
    We need starting row and col because we update
    the matrix A on the go and we need to work only with the
    updated part.
    
*/
void get_slices(const double* A, double* A11, double* A21, double* A22, int N, int r, int i0, int j0)
{
  int len1 = 0, len2 = 0, len3 = 0;
  for (int i = i0; i < N; i++)
  {
    for (int j = j0; j < j0 + r; j++)
    {
      if (i < i0 + r)
      {
        A11[len1++] = A[i * N + j];
      }
      else
      {
        A21[len2++] = A[i * N + j];
      }
    }
    if (i >= i0 + r)
    {
      for (int j = j0 + r; j < N; j++)
      {
        A22[len3++] = A[i * N + j];
      }
    }
  }
}

/*
  Composes L out of blocks.
  L [out] - L matrix to compose.
  L11 [in] - L11 block (top left).
  L21 [in] - L21 block (bottom left).
  r [in] - Block size.
  N [in] - Dimension of square matrix L.
  i0 [in] - Starting row to count from.
  j0 [in] - Starting col to count from.
  Remarks:
    We need starting row and col because we update
    the matrix L on the go and we need to work only with the
    part that hasn't been updated.
*/
void composeL(double* L, double* L11, double* L21, int r, int N, int i0, int j0)
{
  int len1 = 0, len2 = 0;
  for (int i = i0; i < i0 + r; ++i)
  {
    for (int j = j0; j < j0 + r; ++j)
    {
      L[i * N + j] = L11[len1++];
    }
  }
  for (int i = i0 + r; i < N; ++i)
  {
    for (int j = j0; j < j0 + r; ++j)
    {
      L[i * N + j] = L21[len2++];
    }
  }
}



/*
  Calculates A22 tilde matrix. (A22 - L21*L21^T).
  A [in] - Matrix A to take elements from.
  A22 [in] - A22 matrix.
  L21 [in] - L21 matrix.
  N [in] - Dimension of the square matrix A.
  Ns [in] - Number of rows in L21.
  r [in] - Number of cols in L21.
  i0 [in] - Starting row to count from.
  j0 [in] - Starting col to count from.
  Remarks:
    We need starting row and col because we update
    the matrix A on the go and we need to work only with the
    updated part.
*/
void getA22tilde(double* A, double* A22, double* L21, int N, int Ns, int r, int i0, int j0)
{
  double* L21_transp = transpose(L21, Ns, r);
  double* product = new double[Ns * Ns];
  memset(product, 0, Ns * Ns * sizeof(*product));
  double* A22_tilde = new double[Ns * Ns];
  prod(L21, L21_transp, product, Ns, r);

  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < Ns; ++j)
    {
      A22_tilde[i * Ns + j] = A22[i * Ns + j] - product[i * Ns + j];
    }
  }
  int len1 = 0;
  for (int i = i0; i < N; ++i)
  {
    for (int j = j0; j < N; ++j)
    {
      A[i * N + j] = A22_tilde[len1++];
    }
  }
  
  delete[] A22_tilde;
  delete[] product;
  delete[] L21_transp;
}

bool is_lower(double* L, int N)
{
  bool isLower = true;
  for (int row = 0; row < N; row++)
  {
    for (int col = 0; col < N; col++)
    {
      if (col > row && L[row * N + col] != 0)
      {
        isLower = false;
        break;
      }
    }
  }
  return isLower;
}
/*
  Block Cholesky decomposition. (A = L * L^T)
  A [in] - A matrix of the system.
  L [out] - Resulting L matrix.
  N [in] - Dimension of square matrix A.
*/
void Cholesky_Decomposition(double* A, double* L, int N)
{
  double* A11 = nullptr, * A21 = nullptr, * A22 = nullptr, * L11 = nullptr, * L21 = nullptr, * L22 = nullptr;
  int r = 500;
  int Ns = N;
  int numOfIterations = N / r;
  int i0 = 0;
  int j0 = 0;

  if ((r >= N) || (N - r < 2))
  {
    Cholesky_Decomposition_standard(A, L, N);
    return;
  }

  for (int i = 0; i < numOfIterations; ++i)
  { 
    A11 = new double[r * r];
    memset(A11, 0, (r * r) * sizeof(*A11));
    A21 = new double[(Ns - r) * r];
    memset(A21, 0, ((Ns - r) * r) * sizeof(*A21));
    A22 = new double[(Ns - r) * (Ns - r)];
    memset(A22, 0, ((Ns - r) * (Ns - r)) * sizeof(*A22));
    L11 = new double[r * r];
    memset(L11, 0, (r * r) * sizeof(*L11));
    L21 = new double[(Ns - r) * r];
    memset(L21, 0, ((Ns - r) * r) * sizeof(*L21));
    
    i0 = i * r;
    j0 = i * r;
    get_slices(A, A11, A21, A22, N, r, i0, j0);

    L11 = get_L11(A11, r);
    L21 = get_L21(L11, A21, Ns - r, r);

    composeL(L, L11, L21, r, N, i0, j0); 

    getA22tilde(A, A22, L21, N, Ns-r, r, i0+r, j0+r);
    Ns = N - r * (i + 1);
  }
  
  delete[] A11;
  delete[] A21;
  delete[] L11;

  L22 = new double[Ns * Ns];
  double* A22_tilde = new double[Ns * Ns];
  getA22tilde(A22_tilde, A22, L21, Ns, Ns, r, 0, 0);
  Cholesky_Decomposition_standard(A22_tilde, L22, Ns);

  i0 += r;
  j0 += r;
  int len1 = 0;
  for (int i = i0; i < N; ++i)
  {
    for (int j = j0; j < N; ++j)
    {
      L[i * N + j] = L22[len1++];
    }
  }
  delete[] A22_tilde;
  delete[] A22;
  delete[] L21;
  delete[] L22;
}

/*
  Prints the original matrix and matrix L.
  matrix [in] - Original A matrix.
  N [in] - Dimension of square matrix.
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

bool compare(double* A1, double* A2, int N)
{
  bool equal = true;
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      if (fabs(A1[i * N + j] - A2[i * N + j]) > 0.001)
      {
        equal = false;
      }
    }
  }
  return equal;
}
// Driver Code
int main()
{
  srand(time(NULL));
  const int N = 1000;
  time_t begin, end;

  double* matrix = generate_matrix(N);
  double* matrix2 = new double[N * N];
  copy(matrix, matrix + N * N, matrix2);
  double* L = new double [N * N];
  double* L2 = new double[N * N];
  memset(L, 0, (N*N) * sizeof(* L));
  memset(L2, 0, (N * N) * sizeof(*L2));
  begin = clock();
  Cholesky_Decomposition(matrix, L, N);
  end = clock();

  //print_result(matrix2, N, L);
  cout << "Total time: " << (end - begin) / 1000.0 << endl;
  cout << "Compare with standard" << endl;

  Cholesky_Decomposition_standard(matrix2, L2, N);

  cout << "Answers are equal: " << compare(L, L2, N) << endl;
  delete[] L;
  delete[] matrix;
  delete[] matrix2;
  return 0;
}