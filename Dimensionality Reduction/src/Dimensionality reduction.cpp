#include <iostream>
#include <fstream>
#include <vector>

#include "poissonsolver.h"

using namespace std;

void TestCase1()
{
    // Initializing the grid
    UniformGrid grid;
    grid.LX = 1.0;
    grid.LY = 1.0;
    grid.nCellsX = 200;
    grid.nCellsY = 200;
    grid.nRegsX = 50;
    grid.nRegsY = 50;

    double hx = grid.LX / grid.nCellsX;
    double hy = grid.LY / grid.nCellsY;

    // Setting the boundary conditions in the cell centers
    BCTypes bcTypes;
    bcTypes.alongX = BCType::Dirichlet;
    bcTypes.alongY = BCType::Dirichlet;

    BoundaryConditions bc;

    bc.left = vector<double>(grid.nCellsY);
    bc.right = vector<double>(grid.nCellsY);
    for (int i = 0; i < grid.nCellsY; i++)
    {
        double y = hy * (i + 0.5);

        bc.left[i] = sin(pi * y);
        bc.right[i] = exp(pi * grid.LX) * sin(pi * y) + 0.5 * y * y;
    }

    bc.bottom = vector<double>(grid.nCellsX);
    bc.top = vector<double>(grid.nCellsX);
    for (int i = 0; i < grid.nCellsX; i++)
    {
        double x = hx * (i + 0.5);

        bc.bottom[i] = 0;
        bc.top[i] = 0.5 * x * x;
    }

    // Setting the RHS in the cell centers
    vector<double> rhsUnfolded(grid.nCellsX * grid.nCellsY);
    for (int j = 0; j < grid.nCellsY; j++)
    {
        for (int i = 0; i < grid.nCellsX; i++)
        {
            double x = hx * (i + 0.5);
            double y = hy * (j + 0.5);
            int ind = i + j * grid.nCellsX;

            rhsUnfolded[ind] = x * x + y * y;
        }
    }

    // Solving the equation
    PoissonSolver solver(grid, bcTypes, bc);
    vector<double> solutionUnfolded = solver.Solve(rhsUnfolded);

    // Exporting the solution
    std::ofstream output;
    output.open("output (case 1).txt");

    for (int j = 0; j < grid.nCellsY; j++)
    {

        for (int i = 0; i < grid.nCellsX; i++)
        {
            double x = hx * (i + 0.5);
            double y = hy * (j + 0.5);
            int ind = i + j * grid.nCellsX;

            output << x << " " << y << " " << solutionUnfolded[ind] << endl;
        }
    }
    output.close();
}

void TestCase2()
{
    // Initializing the grid
    UniformGrid grid;
    grid.LX = 1.0;
    grid.LY = 1.0;
    grid.nCellsX = 200;
    grid.nCellsY = 200;
    grid.nRegsX = 20;
    grid.nRegsY = 20;

    double hx = grid.LX / grid.nCellsX;
    double hy = grid.LY / grid.nCellsY;

    // Setting the boundary conditions in the cell centers
    BCTypes bcTypes;
    bcTypes.alongX = BCType::Dirichlet;
    bcTypes.alongY = BCType::Periodic;

    BoundaryConditions bc;

    bc.left = vector<double>(grid.nCellsY, 0);
    bc.right = vector<double>(grid.nCellsY, 0);

    // Setting the charge density in the cell centers
    vector<double> rhsUnfolded(grid.nCellsX * grid.nCellsY);
    for (int j = 0; j < grid.nCellsY; j++)
    {

        for (int i = 0; i < grid.nCellsX; i++)
        {
            double x = hx * (i + 0.5);
            double y = hy * (j + 0.5);
            int ind = i + j * grid.nCellsX;

            rhsUnfolded[ind] = y * sin(5 * pi * x) + exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 0.02);
        }
    }

    // Solving the equation
    PoissonSolver solver(grid, bcTypes, bc);
    vector<double> solutionUnfolded = solver.Solve(rhsUnfolded);

    // Exporting the solution
    std::ofstream output;
    output.open("output (case 2).txt");

    for (int j = 0; j < grid.nCellsY; j++)
    {

        for (int i = 0; i < grid.nCellsX; i++)
        {
            double x = hx * (i + 0.5);
            double y = hy * (j + 0.5);
            int ind = i + j * grid.nCellsX;

            output << x << " " << y << " " << solutionUnfolded[ind] << endl;
        }
    }
    output.close();
}

void TestCase3()
{
    // Initializing the grid
    UniformGrid grid;
    grid.LX = 1.0;
    grid.LY = 1.0;
    grid.nCellsX = 200;
    grid.nCellsY = 200;
    grid.nRegsX = 50;
    grid.nRegsY = 50;

    double hx = grid.LX / grid.nCellsX;
    double hy = grid.LY / grid.nCellsY;

    // Setting the boundary conditions in the cell centers
    BCTypes bcTypes;
    bcTypes.alongX = BCType::Dirichlet;
    bcTypes.alongY = BCType::Dirichlet;

    BoundaryConditions bc;

    bc.left = vector<double>(grid.nCellsY);
    bc.right = vector<double>(grid.nCellsY);
    for (int i = 0; i < grid.nCellsY; i++)
    {
        double y = hy * (i + 0.5);

        bc.left[i] = sin(3 * pi * y);
        bc.right[i] = sin(5 * pi + 3 * pi * y);
    }

    bc.bottom = vector<double>(grid.nCellsX);
    bc.top = vector<double>(grid.nCellsX);
    for (int i = 0; i < grid.nCellsX; i++)
    {
        double x = hx * (i + 0.5);

        bc.bottom[i] = sin(5 * pi * x);
        bc.top[i] = sin(5 * pi * x + 3 * pi);
    }

    // Setting the charge density in the cell centers
    vector<double> rhsUnfolded(grid.nCellsX * grid.nCellsY);
    for (int j = 0; j < grid.nCellsY; j++)
    {
        for (int i = 0; i < grid.nCellsX; i++)
        {
            double x = hx * (i + 0.5);
            double y = hy * (j + 0.5);
            int ind = i + j * grid.nCellsX;

            rhsUnfolded[ind] = -34 * pi * pi * sin(5 * pi * x + 3 * pi * y);
        }
    }

    // Solving the equation
    PoissonSolver solver(grid, bcTypes, bc);
    vector<double> solutionUnfolded = solver.Solve(rhsUnfolded);

    // Exporting the solution
    std::ofstream output;
    output.open("output (case 3).txt");

    for (int j = 0; j < grid.nCellsY; j++)
    {

        for (int i = 0; i < grid.nCellsX; i++)
        {
            double x = hx * (i + 0.5);
            double y = hy * (j + 0.5);
            int ind = i + j * grid.nCellsX;

            output << x << " " << y << " " << solutionUnfolded[ind] << endl;
        }
    }
    output.close();
}

int main()
{
    try
    {
        TestCase1();
        TestCase2();
        TestCase3();
    }
    catch (const std::exception& e)
    {
        std::cout << e.what();
        exit(-1);
    }
}