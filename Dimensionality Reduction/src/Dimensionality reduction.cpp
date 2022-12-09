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

    // Setting the boundary conditions in the nodes
    BCTypes bcTypes;
    bcTypes.alongX = BCType::Dirichlet;
    bcTypes.alongY = BCType::Dirichlet;

    BoundaryConditions bc;

    bc.left = vector<double>(grid.nCellsY + 1);
    bc.right = vector<double>(grid.nCellsY + 1);
    for (int i = 0; i < grid.nCellsY + 1; i++)
    {
        double y = hy * i;

        bc.left[i] = sin(pi * y);
        bc.right[i] = exp(pi * grid.LX) * sin(pi * y) + 0.5 * y * y;
    }

    bc.bottom = vector<double>(grid.nCellsX + 1);
    bc.top = vector<double>(grid.nCellsX + 1);
    for (int i = 0; i < grid.nCellsX + 1; i++)
    {
        double x = hx * i;

        bc.bottom[i] = 0;
        bc.top[i] = 0.5 * x * x;
    }

    // Setting the RHS in the nodes
    vector<double> rhsUnfolded((grid.nCellsX + 1) * (grid.nCellsY + 1));
    for (int j = 0; j < grid.nCellsY + 1; j++)
    {
        for (int i = 0; i < grid.nCellsX + 1; i++)
        {
            double x = hx * i;
            double y = hy * j;
            int ind = i + j * (grid.nCellsX + 1);

            rhsUnfolded[ind] = x * x + y * y;
        }
    }

    // Solving the equation
    PoissonSolver solver(grid, bcTypes, bc);
    vector<double> solutionUnfolded = solver.Solve(rhsUnfolded);

    // Exporting the solution
    std::ofstream output;
    output.open("output (case 1).txt");

    for (int j = 0; j < grid.nCellsY + 1; j++)
    {

        for (int i = 0; i < grid.nCellsX + 1; i++)
        {
            double x = hx * i;
            double y = hy * j;
            int ind = i + j * (grid.nCellsX + 1);

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

    // Setting the boundary conditions in the nodes
    BCTypes bcTypes;
    bcTypes.alongX = BCType::Dirichlet;
    bcTypes.alongY = BCType::Periodic;

    BoundaryConditions bc;

    bc.left = vector<double>(grid.nCellsY + 1, 0);
    bc.right = vector<double>(grid.nCellsY + 1, 0);

    // Setting the charge density in the nodes
    vector<double> rhsUnfolded((grid.nCellsX + 1) * (grid.nCellsY + 1));
    for (int j = 0; j < grid.nCellsY + 1; j++)
    {
        for (int i = 0; i < grid.nCellsX + 1; i++)
        {
            double x = hx * i;
            double y = hy * j;
            int ind = i + j * (grid.nCellsX + 1);

            rhsUnfolded[ind] = y * sin(5 * pi * x) + exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 0.02);
        }
    }

    // Solving the equation
    PoissonSolver solver(grid, bcTypes, bc);
    vector<double> solutionUnfolded = solver.Solve(rhsUnfolded);

    // Exporting the solution
    std::ofstream output;
    output.open("output (case 2).txt");

    for (int j = 0; j < grid.nCellsY + 1; j++)
    {

        for (int i = 0; i < grid.nCellsX + 1; i++)
        {
            double x = hx * i;
            double y = hy * j;
            int ind = i + j * (grid.nCellsX + 1);

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

    // Setting the boundary conditions in the nodes
    BCTypes bcTypes;
    bcTypes.alongX = BCType::Dirichlet;
    bcTypes.alongY = BCType::Dirichlet;

    BoundaryConditions bc;

    bc.left = vector<double>(grid.nCellsY + 1);
    bc.right = vector<double>(grid.nCellsY + 1);
    for (int i = 0; i < grid.nCellsY + 1; i++)
    {
        double y = hy * i;

        bc.left[i] = sin(3 * pi * y);
        bc.right[i] = sin(5 * pi + 3 * pi * y);
    }

    bc.bottom = vector<double>(grid.nCellsX + 1);
    bc.top = vector<double>(grid.nCellsX + 1);
    for (int i = 0; i < grid.nCellsX + 1; i++)
    {
        double x = hx * i;

        bc.bottom[i] = sin(5 * pi * x);
        bc.top[i] = sin(5 * pi * x + 3 * pi);
    }

    // Setting the charge density in the nodes
    vector<double> rhsUnfolded((grid.nCellsX + 1) * (grid.nCellsY + 1));
    for (int j = 0; j < grid.nCellsY + 1; j++)
    {
        for (int i = 0; i < grid.nCellsX + 1; i++)
        {
            double x = hx * i;
            double y = hy * j;
            int ind = i + j * (grid.nCellsX + 1);

            rhsUnfolded[ind] = -34 * pi * pi * sin(5 * pi * x + 3 * pi * y);
        }
    }

    // Solving the equation
    PoissonSolver solver(grid, bcTypes, bc);
    vector<double> solutionUnfolded = solver.Solve(rhsUnfolded);

    // Exporting the solution
    std::ofstream output;
    output.open("output (case 3).txt");

    for (int j = 0; j < grid.nCellsY + 1; j++)
    {

        for (int i = 0; i < grid.nCellsX + 1; i++)
        {
            double x = hx * i;
            double y = hy * j;
            int ind = i + j * (grid.nCellsX + 1);

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