#include <cmath>
#include <iostream>
#include <iomanip>

#define DIMENSION 4

void decomposition();
void gaussian_elimination();
bool zero(double val);
void swap_with_max(int index);
bool zero(double val);
int max_in_column(int col);
void print_matrix(double matrix[DIMENSION][DIMENSION]);
void show_results();
void solveAxB();

double matrix_A[DIMENSION][DIMENSION] = {
    1.00, -20.00, 30.00, -4.00,
    2.00, -40.00, -6.00, 50.0,
    9.00, -180.0, 11.00, -12.0,
    -16.0, 15.00, -140.0, 13.00};

double vector_B[DIMENSION] = {
    35.0,
    104.0,
    -366.0,
    -354.0};

int rows_order[] = {0, 1, 2, 3};

double matrix_L[DIMENSION][DIMENSION] = {{0.0}};
double matrix_U[DIMENSION][DIMENSION] = {{0.0}};

double vector_X[DIMENSION];
double vector_Y[DIMENSION];

int main()
{
    decomposition();
    solveAxB();
}

void solveAxB()
{

    for (int i = 0; i < DIMENSION; i++)
    {
        double sum = 0.0;
        for (int j = 0; j <= i - 1; j++)
        {
            sum += matrix_A[rows_order[i]][j] * vector_Y[rows_order[j]];
        }
        vector_Y[rows_order[i]] = (vector_B[rows_order[i]] - sum);
    }

    for (int i = DIMENSION - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < DIMENSION; j++)
        {
            sum += matrix_A[rows_order[i]][j] * vector_X[rows_order[j]];
        }
        vector_X[rows_order[i]] = (vector_Y[rows_order[i]] - sum) / matrix_A[rows_order[i]][i];
    }


    std::cout<<std::endl<<"Wektor Y: "<<std::endl<<std::endl;

    for (int i = 0; i < DIMENSION; i++)
            std::cout << vector_Y[rows_order[i]] << std::endl;


    std::cout<<std::endl<<"Wektor X: "<<std::endl<<std::endl;

    for (int i = 0; i < DIMENSION; i++)
        std::cout << vector_X[rows_order[i]] << std::endl;
}

void decomposition()
{

    std::cout
        << std::endl
        << std::endl
        << "Macierz A na samym poczatku: ";

    print_matrix(matrix_A);

    std::cout
        << "Zaczynam dekompozycje macierzy LU metoda Gaussa...";

    gaussian_elimination();

    std::cout
        << "Zakonczylem dekompozycje macierzy LU metoda Gaussa.";

    // Wypelniam macierz L jedynkami na przekatnej:

    for (int i = 0; i < DIMENSION; i++)
        matrix_L[rows_order[i]][i] = 1;

    // Wyprowadzam wyniki:
    show_results();
}

void gaussian_elimination()
{
    for (int i = 0; i < DIMENSION; i++)
    {
        print_matrix(matrix_A);

        int row = rows_order[i];
        int column = i;

        if (zero(matrix_A[row][column]))
        {
            swap_with_max(i);
            row = rows_order[i];
        }

        for (int j = i + 1; j < DIMENSION; j++)
        {

            int next_row = rows_order[j];

            if (zero(matrix_A[next_row][column]))
            {
                matrix_L[next_row][column] = 0;
            }
            else
            {
                double quotient = matrix_A[next_row][column] / matrix_A[row][column];

                matrix_A[next_row][column] = 0.0;
                matrix_L[next_row][column] = quotient;
                vector_B[next_row] = vector_B[next_row] - quotient * vector_B[row];

                for (int k = i + 1; k < DIMENSION; k++)
                {
                    int next_col = k;
                    matrix_A[next_row][next_col] = matrix_A[next_row][next_col] - (quotient * matrix_A[row][next_col]);
                }
            }
        }
    }
}

bool zero(double val)
{
    return (val == 0);
}

void swap_with_max(int index)
{
    int max = max_in_column(index);
    std::swap(rows_order[index], rows_order[max]);
}

int max_in_column(int col)
{
    int max = col + 1;
    for (int i = col + 2; i < DIMENSION; i++)
    {
        if (fabs(matrix_A[rows_order[max]][col]) < fabs(matrix_A[i][col]))
            max = i;
    }

    return max;
}

void print_matrix(double matrix[DIMENSION][DIMENSION])
{

    std::cout << std::endl
              << std::endl;

    for (int i = 0; i < DIMENSION; i++)
    {
        int row = rows_order[i];

        for (int j = 0; j < DIMENSION; j++)
        {
            int column = j;
            std::cout << std::setw(10) << matrix[row][column] << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl
              << std::endl;
}

void show_results()
{

    std::cout
        << std::endl
        << std::endl
        << "================   Wyniki   ================= "
        << std::endl
        << std::endl;

    std::cout
        << std::endl
        << "Macierz L: "
        << std::endl
        << std::endl;
    print_matrix(matrix_L);

    std::cout
        << std::endl
        << "Macierz U: "
        << std::endl
        << std::endl;
    print_matrix(matrix_A);

    std::cout<< "Wektor B: "<<std::endl<<std::endl;
    for (int i = 0; i < DIMENSION; i++)
        std::cout << vector_B[rows_order[i]] << std::endl;
}