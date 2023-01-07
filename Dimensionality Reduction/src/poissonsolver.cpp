#include "poissonsolver.h"

#include <stdexcept>
#include <iostream>

PoissonSolver::PoissonSolver(UniformGrid grid, BCTypes bcTypes, BoundaryConditions bc) :
	_grid(grid), _bcTypes(bcTypes)
{
    // Check the validity of arguments
    bool validGrid = grid.LX > 0 &&
        grid.LY > 0 &&
        grid.nCellsX > 0 &&
        grid.nCellsY > 0 &&
        grid.nRegsX > 0 &&
        grid.nRegsY > 0;

    if (!validGrid)
        throw std::invalid_argument("Some grid parameters have not been initialized or are not positive numbers.");

    if ((int)(grid.nCellsX / grid.nRegsX) - (double)grid.nCellsX / grid.nRegsX != 0 ||
        (int)(grid.nCellsY / grid.nRegsY) - (double)grid.nCellsY / grid.nRegsY != 0)
        throw std::invalid_argument("Each region must contain an integer number of cells. Make sure that nCells can be divided by nRegs without remainder.");

    if (bcTypes.alongX == BCType::None ||
        bcTypes.alongY == BCType::None)
        throw std::invalid_argument("Missing boundary condition types.");

    bool validBCX = !((bcTypes.alongX == BCType::Dirichlet) &&
        !(bc.left.size() == _grid.nCellsY + 1 && bc.right.size() == _grid.nCellsY + 1));

    bool validBCY = !((bcTypes.alongY == BCType::Dirichlet) &&
        !(bc.bottom.size() == _grid.nCellsX + 1 && bc.top.size() == _grid.nCellsX + 1));

    if (!validBCX || !validBCY)
        throw std::invalid_argument("Invalid size of the boundary condition vectors.");

    // Calculate the number of cell per region
    _nx = grid.nCellsX / grid.nRegsX;
    _ny = grid.nCellsY / grid.nRegsY;

    if (_nx == 1 || _ny == 1)
        throw std::invalid_argument("The case when the number of cells per region is 1 is currently not supported.");

    // Calculate the number of equations
    if (bcTypes.alongX == BCType::Dirichlet && 
        bcTypes.alongY == BCType::Dirichlet)
        _nEq = grid.nRegsX * grid.nRegsY * (_nx + _ny) + 2 * (grid.nRegsX + grid.nRegsY);

    if (bcTypes.alongX == BCType::Dirichlet && 
        bcTypes.alongY == BCType::Periodic)
        _nEq = grid.nRegsX * grid.nRegsY * (_nx + _ny) + 2 * grid.nRegsY;

    if (bcTypes.alongX == BCType::Periodic && 
        bcTypes.alongY == BCType::Dirichlet)
        _nEq = grid.nRegsX * grid.nRegsY * (_nx + _ny) + 2 * grid.nRegsX;

    if (bcTypes.alongX == BCType::Periodic && 
        bcTypes.alongY == BCType::Periodic)
        _nEq = grid.nRegsX * grid.nRegsY * (_nx + _ny);

    // BC from the nodes to the centers
    if (bcTypes.alongX == BCType::Dirichlet)
    {
        _bc.left = std::vector<double>(grid.nCellsY);
        for (int i = 0; i < grid.nCellsY; i++)
            _bc.left[i] = 0.5 * (bc.left[i] + bc.left[i + 1]);
        
        _bc.right = std::vector<double>(grid.nCellsY);
        for (int i = 0; i < grid.nCellsY; i++)
            _bc.right[i] = 0.5 * (bc.right[i] + bc.right[i + 1]);
    }

    if (bcTypes.alongY == BCType::Dirichlet)
    {
        _bc.bottom = std::vector<double>(grid.nCellsX);
        for (int i = 0; i < grid.nCellsX; i++)
            _bc.bottom[i] = 0.5 * (bc.bottom[i] + bc.bottom[i + 1]);

        _bc.top = std::vector<double>(grid.nCellsX);
        for (int i = 0; i < grid.nCellsX; i++)
            _bc.top[i] = 0.5 * (bc.top[i] + bc.top[i + 1]);
    }

    // Assemble the system
    Assemble();

    // Allocate the memory
    _rhsUnfoldedC = std::vector<double>(_grid.nCellsX * _grid.nCellsY);
    _solution = Eigen::VectorXd(_nEq);
    _solutionUnfolded = std::vector<double>(_grid.nCellsX * _grid.nCellsY);
    _solutionUnfoldedN = std::vector<double>((_grid.nCellsX + 1) * (_grid.nCellsY + 1));
}

int PoissonSolver::Pos(int m, int n, Flag flag, int cell)
{
    switch (flag)
    {
    case Flag::X:
        return (_nx + _ny) * (m + n * _grid.nRegsX) + cell;
        break;
    case Flag::Y:
        return (_nx + _ny) * (m + n * _grid.nRegsX) + _nx + cell;
        break;
    case Flag::Left:
        return _grid.nRegsX * _grid.nRegsY * (_nx + _ny) + n;
        break;
    case Flag::Right:
        return _grid.nRegsX * _grid.nRegsY * (_nx + _ny) + _grid.nRegsY + n;
        break;
    case Flag::Bottom:
        if (_bcTypes.alongX == BCType::Dirichlet)
        {
            return _grid.nRegsX * _grid.nRegsY * (_nx + _ny) + 2 * _grid.nRegsY + m;
        }
        if (_bcTypes.alongX == BCType::Periodic)
        {
            return _grid.nRegsX * _grid.nRegsY * (_nx + _ny) + m;
        }
        break;
    case Flag::Top:
        if (_bcTypes.alongX == BCType::Dirichlet)
        {
            return _grid.nRegsX * _grid.nRegsY * (_nx + _ny) + 2 * _grid.nRegsY + _grid.nRegsX + m;
        }
        if (_bcTypes.alongX == BCType::Periodic)
        {
            return _grid.nRegsX * _grid.nRegsY * (_nx + _ny) + _grid.nRegsX + m;

        }
        break;
    }
}

void PoissonSolver::Assemble()
{
    std::vector<T> coeffs;
    _rhs = Eigen::VectorXd::Constant(_nEq, 0.0);

    int M = _grid.nRegsX;
    int N = _grid.nRegsY;
    int nx = _nx;
    int ny = _ny;
    double Lx = _grid.LX / M;
    double Ly = _grid.LY / N;
    double hx = (double)Lx / nx;
    double hy = (double)Ly / ny;
    BCType bcTypeX = _bcTypes.alongX;
    BCType bcTypeY = _bcTypes.alongY;

    // Internal parameters
    int regA, nRegsA;
    int regB, nRegsB;
    int nCellsA, nCellsB;
    int mPrevA, nPrevA;
    int mNextA, nNextA;
    int mPrevB, nPrevB;
    int mNextB, nNextB;
    int mFirstA, nFirstA;
    int mLastA, nLastA;
    int mFirstB, nFirstB;
    int mLastB, nLastB;
    double LA, LB;
    double hA, hB;
    Flag varA, varB;
    Flag BA1, BA2;
    Flag BB1, BB2;
    Flag bA1, bA2;
    Flag bB1, bB2;
    BCType bcTypeA, bcTypeB;

    auto setInternalParams = [&](int m, int n, Flag var)
    {
        switch (var)
        {
        case Flag::X:
            regA = m, nRegsA = M;
            regB = n, nRegsB = N;
            nCellsA = nx, nCellsB = ny;
            mPrevA = m - 1, nPrevA = n;
            mNextA = m + 1, nNextA = n;
            mPrevB = m, nPrevB = n - 1;
            mNextB = m, nNextB = n + 1;
            mFirstA = 0, nFirstA = n;
            mLastA = M - 1, nLastA = n;
            mFirstB = m, nFirstB = 0;
            mLastB = m, nLastB = N - 1;
            LA = Lx, LB = Ly;
            hA = hx, hB = hy;
            varA = X, varB = Y;
            BA1 = Left, BA2 = Right;
            BB1 = Bottom, BB2 = Top;
            bA1 = Left, bA2 = Right;
            bB1 = Bottom, bB2 = Top;
            bcTypeA = bcTypeX, bcTypeB = bcTypeY;
            break;

        case Y:
            regA = n, nRegsA = N;
            regB = m, nRegsB = M;
            nCellsA = ny, nCellsB = nx;
            mPrevA = m, nPrevA = n - 1;
            mNextA = m, nNextA = n + 1;
            mPrevB = m - 1, nPrevB = n;
            mNextB = m + 1, nNextB = n;
            mFirstA = m, nFirstA = 0;
            mLastA = m, nLastA = N - 1;
            mFirstB = 0, nFirstB = n;
            mLastB = M - 1, nLastB = n;
            LA = Ly, LB = Lx;
            hA = hy, hB = hx;
            varA = Y, varB = X;
            BA1 = Bottom, BA2 = Top;
            BB1 = Left, BB2 = Right;
            bA1 = Bottom, bA2 = Top;
            bB1 = Left, bB2 = Right;
            bcTypeA = bcTypeY, bcTypeB = bcTypeX;
            break;
        }
    };

    auto BoundaryValue = [&](Flag bFlag, double z)
    {
        int ind;
        switch (bFlag)
        {
        case Left:
            ind = (z - 0.5 * hy) / hy;
            return _bc.left[ind];
            break;
        case Right:
            ind = (z - 0.5 * hy) / hy;
            return _bc.right[ind];
            break;
        case Bottom:
            ind = (z - 0.5 * hx) / hx;
            return _bc.bottom[ind];
            break;
        case Top:
            ind = (z - 0.5 * hx) / hx;
            return _bc.top[ind];
            break;
        }
    };

    // d^2 A / da^2
    auto Term1 = [&](int l, int m, int n, int cell, Flag var, std::vector<T>& coeffs)
    {
        setInternalParams(m, n, var);

        if (cell == 0)
        {
            if (regA == 0)
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -4 * LB / (hA * hA)));
                coeffs.push_back(T(l, Pos(m, n, varA, cell + 1), 4.0 / 3 * LB / (hA * hA)));
                if (bcTypeA == BCType::Dirichlet)
                    coeffs.push_back(T(l, Pos(m, n, BA1), 8.0 / 3 * LB / (hA * hA)));
            }
            else
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -2 * LB / (hA * hA)));
                coeffs.push_back(T(l, Pos(m, n, varA, cell + 1), LB / (hA * hA)));
                coeffs.push_back(T(l, Pos(mPrevA, nPrevA, varA, nCellsA - 1), LB / (hA * hA)));
            }
        }
        else if (cell == nCellsA - 1)
        {
            if (regA == nRegsA - 1)
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -4 * LB / (hA * hA)));
                coeffs.push_back(T(l, Pos(m, n, varA, cell - 1), 4.0 / 3 * LB / (hA * hA)));
                if (bcTypeA == BCType::Dirichlet)
                    coeffs.push_back(T(l, Pos(m, n, BA2), 8.0 / 3 * LB / (hA * hA)));
            }
            else
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -2 * LB / (hA * hA)));
                coeffs.push_back(T(l, Pos(m, n, varA, cell - 1), LB / (hA * hA)));
                coeffs.push_back(T(l, Pos(mNextA, nNextA, varA, 0), LB / (hA * hA)));
            }
        }
        else
        {
            coeffs.push_back(T(l, Pos(m, n, varA, cell), -2 * LB / (hA * hA)));
            coeffs.push_back(T(l, Pos(m, n, varA, cell + 1), LB / (hA * hA)));
            coeffs.push_back(T(l, Pos(m, n, varA, cell - 1), LB / (hA * hA)));
        }
        return;
    };

    // \int (dB/da) db
    auto Term2 = [&](int l, int m, int n, int cell, Flag var, std::vector<T>& coeffs, Eigen::VectorXd& RHS)
    {
        setInternalParams(m, n, var);

        if (cell == 0)
        {
            for (int i = 0; i < nCellsB; i++)
            {
                if (regA == 0)
                {
                    if (bcTypeA == BCType::Dirichlet)
                    {
                        coeffs.push_back(T(l, Pos(m, n, varB, i), -2 * hB / (hA * hA)));
                        RHS(l) += -BoundaryValue(bA1, (regB * nCellsB + i + 0.5) * hB) * 2 * hB / (hA * hA);
                    }
                    if (bcTypeA == BCType::Periodic)
                    {
                        coeffs.push_back(T(l, Pos(m, n, varB, i), -hB / (hA * hA)));
                        coeffs.push_back(T(l, Pos(mLastA, nLastA, varB, i), hB / (hA * hA)));
                    }
                }
                else
                {
                    coeffs.push_back(T(l, Pos(m, n, varB, i), -hB / (hA * hA)));
                    coeffs.push_back(T(l, Pos(mPrevA, nPrevA, varB, i), hB / (hA * hA)));
                }
            }
        }
        if (cell == nCellsA - 1)
        {
            for (int i = 0; i < nCellsB; i++)
            {
                if (regA == nRegsA - 1)
                {
                    if (bcTypeA == BCType::Dirichlet)
                    {
                        coeffs.push_back(T(l, Pos(m, n, varB, i), -2 * hB / (hA * hA)));
                        RHS(l) += -BoundaryValue(bA2, (regB * nCellsB + i + 0.5) * hB) * 2 * hB / (hA * hA);
                    }
                    if (bcTypeA == BCType::Periodic)
                    {
                        coeffs.push_back(T(l, Pos(m, n, varB, i), -hB / (hA * hA)));
                        coeffs.push_back(T(l, Pos(mFirstA, nFirstA, varB, i), hB / (hA * hA)));
                    }
                }
                else
                {
                    coeffs.push_back(T(l, Pos(m, n, varB, i), -hB / (hA * hA)));
                    coeffs.push_back(T(l, Pos(mNextA, nNextA, varB, i), hB / (hA * hA)));
                }
            }
        }
        return;
    };

    // (dB/db) |_bB1^bB2
    auto Term3 = [&](int l, int m, int n, Flag var, std::vector<T>& coeffs)
    {
        setInternalParams(m, n, var);

        if (regB == 0)
        {
            if (bcTypeB == BCType::Dirichlet)
            {
                coeffs.push_back(T(l, Pos(m, n, varB, nCellsB - 1), -1 / hB));
                coeffs.push_back(T(l, Pos(mNextB, nNextB, varB, 0), 1 / hB));
                coeffs.push_back(T(l, Pos(m, n, varB, 0), -2 / hB));
                coeffs.push_back(T(l, Pos(m, n, BB1), 2 / hB));
            }
            if (bcTypeB == BCType::Periodic)
            {
                coeffs.push_back(T(l, Pos(m, n, varB, nCellsB - 1), -1 / hB));
                coeffs.push_back(T(l, Pos(mNextB, nNextB, varB, 0), 1 / hB));
                coeffs.push_back(T(l, Pos(m, n, varB, 0), -1 / hB));
                coeffs.push_back(T(l, Pos(mLastB, nLastB, varB, nCellsB - 1), 1 / hB));
            }
        }
        else if (regB == nRegsB - 1)
        {
            if (bcTypeB == BCType::Dirichlet)
            {
                coeffs.push_back(T(l, Pos(m, n, varB, nCellsB - 1), -2 / hB));
                if (bcTypeB == BCType::Dirichlet)
                    coeffs.push_back(T(l, Pos(m, n, BB2), 2 / hB));
                coeffs.push_back(T(l, Pos(m, n, varB, 0), -1 / hB));
                coeffs.push_back(T(l, Pos(mPrevB, nPrevB, varB, nCellsB - 1), 1 / hB));
            }
            if (bcTypeB == BCType::Periodic)
            {
                coeffs.push_back(T(l, Pos(m, n, varB, nCellsB - 1), -1 / hB));
                coeffs.push_back(T(l, Pos(mFirstB, nFirstB, varB, 0), 1 / hB));
                coeffs.push_back(T(l, Pos(m, n, varB, 0), -1 / hB));
                coeffs.push_back(T(l, Pos(mPrevB, nPrevB, varB, nCellsB - 1), 1 / hB));
            }
        }
        else
        {
            coeffs.push_back(T(l, Pos(m, n, varB, nCellsB - 1), -1 / hB));
            coeffs.push_back(T(l, Pos(mNextB, nNextB, varB, 0), 1 / hB));
            coeffs.push_back(T(l, Pos(m, n, varB, 0), -1 / hB));
            coeffs.push_back(T(l, Pos(mPrevB, nPrevB, varB, nCellsB - 1), 1 / hB));
        }
        return;
    };

    // (dA/db) |_bB1^bB2
    auto Term4 = [&](int l, int m, int n, int cell, Flag var, std::vector<T>& coeffs, Eigen::VectorXd& RHS)
    {
        setInternalParams(m, n, var);

        if (regB == 0)
        {
            if (bcTypeB == BCType::Dirichlet)
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
                coeffs.push_back(T(l, Pos(mNextB, nNextB, varA, cell), 1 / hB));
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -2 / hB));
                RHS(l) += -BoundaryValue(bB1, (regA * nCellsA + cell + 0.5) * hA) * 2 / hB;
            }
            if (bcTypeB == BCType::Periodic)
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
                coeffs.push_back(T(l, Pos(mNextB, nNextB, varA, cell), 1 / hB));
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
                coeffs.push_back(T(l, Pos(mLastB, nLastB, varA, cell), 1 / hB));
            }
        }
        else if (regB == nRegsB - 1)
        {
            if (bcTypeB == BCType::Dirichlet)
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -2 / hB));
                if (bcTypeB == BCType::Dirichlet)
                    RHS(l) += -BoundaryValue(bB2, (regA * nCellsA + cell + 0.5) * hA) * 2 / hB;
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
                coeffs.push_back(T(l, Pos(mPrevB, nPrevB, varA, cell), 1 / hB));
            }
            if (bcTypeB == BCType::Periodic)
            {
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
                coeffs.push_back(T(l, Pos(mFirstB, nFirstB, varA, cell), 1 / hB));
                coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
                coeffs.push_back(T(l, Pos(mPrevB, nPrevB, varA, cell), 1 / hB));
            }
        }
        else
        {
            coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
            coeffs.push_back(T(l, Pos(mNextB, nNextB, varA, cell), 1 / hB));
            coeffs.push_back(T(l, Pos(m, n, varA, cell), -1 / hB));
            coeffs.push_back(T(l, Pos(mPrevB, nPrevB, varA, cell), 1 / hB));
        }
        return;
    };

    // A(boundary) + \int (B/LB) db = \int (phi(boundary)/LB) db
    auto BoundaryEqs = [&](int& l, Flag bFlag, std::vector<T>& coeffs, Eigen::VectorXd& RHS)
    {
        switch (bFlag)
        {
        case Left:
            for (int n = 0; n < N; n++)
            {
                l++;
                coeffs.push_back(T(l, Pos(0, n, Left), 1));
                for (int i = 0; i < ny; i++)
                    coeffs.push_back(T(l, Pos(0, n, Y, i), hy / Ly));

                for (int i = 0; i < ny; i++)
                {
                    RHS(l) += _bc.left[n * ny + i] * hy / Ly;
                }
            }
            break;
        case Right:
            for (int n = 0; n < N; n++)
            {
                l++;
                coeffs.push_back(T(l, Pos(0, n, Right), 1));
                for (int i = 0; i < ny; i++)
                    coeffs.push_back(T(l, Pos(M - 1, n, Y, i), hy / Ly));

                for (int i = 0; i < ny; i++)
                {
                    RHS(l) += _bc.right[n * ny + i] * hy / Ly;
                }
            }
            break;
        case Bottom:
            for (int m = 0; m < M; m++)
            {
                l++;
                coeffs.push_back(T(l, Pos(m, 0, Bottom), 1));
                for (int j = 0; j < nx; j++)
                    coeffs.push_back(T(l, Pos(m, 0, X, j), hx / Lx));

                for (int j = 0; j < nx; j++)
                {
                    RHS(l) += _bc.bottom[m * nx + j] * hx / Lx;
                }
            }
            break;
        case Top:
            for (int m = 0; m < M; m++)
            {
                l++;
                coeffs.push_back(T(l, Pos(m, 0, Top), 1));
                for (int j = 0; j < nx; j++)
                    coeffs.push_back(T(l, Pos(m, N - 1, X, j), hx / Lx));

                for (int j = 0; j < nx; j++)
                {
                    RHS(l) += _bc.top[m * nx + j] * hx / Lx;
                }
            }
            break;
        }
        return;
    };

    int l = -1;
    for (int n = 0; n < _grid.nRegsY; n++)
    {
        for (int m = 0; m < _grid.nRegsX; m++)
        {
            // A = X, B = Y
            for (int j = 0; j < _nx; j++)
            {
                l++;
                Term1(l, m, n, j, Flag::X, coeffs);
                Term2(l, m, n, j, Flag::X, coeffs, _rhs);
                Term3(l, m, n, Flag::X, coeffs);
                Term4(l, m, n, j, Flag::X, coeffs, _rhs);
            }

            // A = Y, B = X
            for (int i = 0; i < _ny; i++)
            {
                l++;
                Term1(l, m, n, i, Flag::Y, coeffs);
                Term2(l, m, n, i, Flag::Y, coeffs, _rhs);
                Term3(l, m, n, Flag::Y, coeffs);
                Term4(l, m, n, i, Flag::Y, coeffs, _rhs);
            }
        }
    }

    if (_bcTypes.alongX == BCType::Dirichlet)
    {
        BoundaryEqs(l, Flag::Left, coeffs, _rhs);
        BoundaryEqs(l, Flag::Right, coeffs, _rhs);
    }
    if (_bcTypes.alongY == BCType::Dirichlet)
    {
        BoundaryEqs(l, Flag::Bottom, coeffs, _rhs);
        BoundaryEqs(l, Flag::Top, coeffs, _rhs);
    }

    _system = SpMat(_nEq, _nEq);
    _system.setFromTriplets(coeffs.begin(), coeffs.end());
    _system.makeCompressed();
    _solver.analyzePattern(_system);
    _solver.factorize(_system);
}

std::vector<double> PoissonSolver::Solve(std::vector<double> rhsUnfolded)
{
     if (rhsUnfolded.size() != (_grid.nCellsX + 1) * (_grid.nCellsY + 1))
        throw std::invalid_argument("Invalid size of the RHS vector.");

    for (int j = 0; j < _grid.nCellsY; j++)
    {
        for (int i = 0; i < _grid.nCellsX; i++)
        {
            _rhsUnfoldedC[i + _grid.nCellsX * j] = 0.25 *
                (rhsUnfolded[i + (_grid.nCellsX + 1) * j] +
                    rhsUnfolded[i + 1 + (_grid.nCellsX + 1) * j] +
                    rhsUnfolded[i + (_grid.nCellsX + 1) * (j + 1)] +
                    rhsUnfolded[i + 1 + (_grid.nCellsX + 1) * (j + 1)]);
        }
    }

    _RHS = _rhs;

    int nx = _nx;
    int ny = _ny;
    double Lx = _grid.LX / _grid.nRegsX;
    double Ly = _grid.LY / _grid.nRegsY;
    double hx = (double)Lx / nx;
    double hy = (double)Ly / ny;

    // \int RHS(b) db
    auto rhs = [&](int l, int m, int n, int cell, Flag flag, Eigen::VectorXd& RHS)
    {
        double x, y;
        switch (flag)
        {
        case Flag::X:
            for (int i = 0; i < ny; i++)
            {
                int ind = (m * nx + cell) + _grid.nCellsX * (n * ny + i);
                RHS(l) += _rhsUnfoldedC[ind] * hy;
            }
            break;
        case Flag::Y:
            for (int j = 0; j < nx; j++)
            {
                int ind = (m * nx + j) + _grid.nCellsX * (n * ny + cell);
                RHS(l) += _rhsUnfoldedC[ind] * hx;
            }
            break;
        }
        return;
    };

    int l = -1;
    for (int n = 0; n < _grid.nRegsY; n++)
    {
        for (int m = 0; m < _grid.nRegsX; m++)
        {
            // A = X, B = Y
            for (int j = 0; j < nx; j++)
            {
                l++;
                rhs(l, m, n, j, Flag::X, _RHS);
            }

            // A = Y, B = X
            for (int i = 0; i < ny; i++)
            {
                l++;
                rhs(l, m, n, i, Flag::Y, _RHS);
            }
        }
    }

    _solution = _solver.solve(_RHS);

    for (int n = 0; n < _grid.nRegsY; n++)
    {
        for (int m = 0; m < _grid.nRegsX; m++)
        {
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    int ind = (m * nx + j) + _grid.nCellsX * (n * ny + i);
                    _solutionUnfolded[ind] = _solution(Pos(m, n, Flag::X, j)) + _solution(Pos(m, n, Flag::Y, i));
                }
            }
        }
    }

    _solutionUnfoldedN[0] = _solutionUnfolded[0];
    _solutionUnfoldedN[_grid.nCellsX] = _solutionUnfolded[_grid.nCellsX - 1];
    _solutionUnfoldedN[(_grid.nCellsX + 1) * _grid.nCellsY] = _solutionUnfolded[_grid.nCellsX * (_grid.nCellsY - 1)];
    _solutionUnfoldedN[(_grid.nCellsX + 1) * (_grid.nCellsY + 1) - 1] = _solutionUnfolded[_grid.nCellsX * _grid.nCellsY - 1];

    for (int i = 1; i < _grid.nCellsX; i++)
        _solutionUnfoldedN[i] = 0.5 * 
        (_solutionUnfolded[i - 1] + 
            _solutionUnfolded[i]);
    for (int i = 1; i < _grid.nCellsX; i++)
        _solutionUnfoldedN[i + (_grid.nCellsX + 1) * _grid.nCellsY] = 0.5 * 
        (_solutionUnfolded[i - 1 + _grid.nCellsX * (_grid.nCellsY - 1)] + 
            _solutionUnfolded[i + _grid.nCellsX * (_grid.nCellsY - 1)]);
    for (int i = 1; i < _grid.nCellsY; i++)
        _solutionUnfoldedN[(_grid.nCellsX + 1) * i] = 0.5 * 
        (_solutionUnfolded[_grid.nCellsX * (i - 1)] +
            _solutionUnfolded[_grid.nCellsX * i]);
    for (int i = 1; i < _grid.nCellsY; i++)
        _solutionUnfoldedN[_grid.nCellsX + (_grid.nCellsX + 1) * i] = 0.5 *
        (_solutionUnfolded[_grid.nCellsX - 1 + _grid.nCellsX * (i - 1)] +
            _solutionUnfolded[_grid.nCellsX - 1 + _grid.nCellsX * i]);

    for (int i = 1; i < _grid.nCellsX; i++)
    {
        for (int j = 1; j < _grid.nCellsY; j++)
        {
            _solutionUnfoldedN[i + (_grid.nCellsX + 1) * j] = 0.25 *
                (_solutionUnfolded[i - 1 + _grid.nCellsX * (j - 1)] +
                    _solutionUnfolded[i + _grid.nCellsX * (j - 1)] +
                    _solutionUnfolded[i - 1 + _grid.nCellsX * j] +
                    _solutionUnfolded[i + _grid.nCellsX * j]);
        }
    }

    return _solutionUnfoldedN;
}