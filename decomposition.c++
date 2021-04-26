#include <iostream>
#include <cmath>
#include <iomanip>

#define DIMENSION 4

double *solveAxB(double matrix_A[DIMENSION][DIMENSION], double matrix_L[DIMENSION][DIMENSION], int *rows_order, double vector_B[DIMENSION]);
void decomposition_LU(double matrix_A[DIMENSION][DIMENSION], double matrix_L[DIMENSION][DIMENSION], int *rows_order, double vector_B[DIMENSION]);
void gaussian_elimination(double matrix_A[DIMENSION][DIMENSION], double matrix_L[DIMENSION][DIMENSION], int *rows_order, double vector_B[DIMENSION]);
void swap_with_max(int index, double matrix_A[DIMENSION][DIMENSION], int *rows_order);
bool zero(double val);
int max_in_column(int col, double matrix[DIMENSION][DIMENSION], int *rows_order);
void print_vector(double vector_B[DIMENSION], int *rows_order);
void print_matrix(double matrix[DIMENSION][DIMENSION], int *rows_order);

int main()
{
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

    double matrix_L[DIMENSION][DIMENSION] = {{0.0}};
    double matrix_U[DIMENSION][DIMENSION] = {{0.0}};

    int rows_order[] = {0, 1, 2, 3};

    solveAxB(matrix_A, matrix_L, rows_order, vector_B);
}

double *solveAxB(double matrix_A[DIMENSION][DIMENSION], double matrix_L[DIMENSION][DIMENSION], int *rows_order, double vector_B[DIMENSION])
{
    decomposition_LU(matrix_A, matrix_L, rows_order, vector_B);

    double vector_X[DIMENSION];
    double vector_Y[DIMENSION];

    for (int i = 0; i < DIMENSION; i++)
    {
        double sum = 0.0;
        for (int j = 0; j <= i - 1; j++)
        {
            sum += matrix_A[i][j];
        }
        vector_Y[i] = (vector_B[i] - sum);
    }

    


    std::cout<<" "<<vector_B[0]<<" "<<vector_B[1]<<" "<<vector_B[2]<<" "<<vector_B[3];

    return NULL;
}

void decomposition_LU(double matrix_A[DIMENSION][DIMENSION], double matrix_L[DIMENSION][DIMENSION], int *rows_order, double vector_B[DIMENSION])
{

    gaussian_elimination(matrix_A, matrix_L, rows_order, vector_B);
    print_matrix(matrix_A, rows_order);
}

void gaussian_elimination(double matrix_A[DIMENSION][DIMENSION], double matrix_L[DIMENSION][DIMENSION], int *rows_order, double vector_B[DIMENSION])
{
    for (int i = 0; i < DIMENSION; i++)
    {

        int row = rows_order[i];
        int column = i;

        if (zero(matrix_A[row][column]))
        {
            swap_with_max(i, matrix_A, rows_order);
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

                for (int next_col = column + 1; next_col < DIMENSION; next_col++)
                {
                    matrix_A[next_row][next_col] = matrix_A[next_row][next_col] - (quotient * matrix_A[row][next_col]);
                }
            }
        }
    }
}

void swap_with_max(int index, double matrix_A[DIMENSION][DIMENSION], int *rows_order)
{
    int max = max_in_column(index, matrix_A, rows_order);
    std::swap(rows_order[index], rows_order[max]);
}

bool zero(double val)
{
    return (val == 0);
}

int max_in_column(int col, double matrix[DIMENSION][DIMENSION], int *rows_order)
{
    int max = col + 1;
    for (int i = col + 2; i < DIMENSION; i++)
    {
        if (fabs(matrix[rows_order[max]][col]) < fabs(matrix[i][col]))
            max = i;
    }

    return max;
}

void print_vector(double vector_B[DIMENSION], int *rows_order)
{

    for (int i = 0; i < DIMENSION; i++)
        std::cout << vector_B[rows_order[i]] << std::endl;
}

void print_matrix(double matrix[DIMENSION][DIMENSION], int *rows_order)
{
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
}