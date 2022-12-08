#pragma once

#include <vector>

#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

const double pi = 3.141592;

enum class BCType
{
	None,
	Dirichlet,
	Periodic
};

class UniformGrid
{
public:
	double	LX = 0;
	double	LY = 0;
	int		nCellsX = 0;
	int		nCellsY = 0;
	int		nRegsX = 0;
	int		nRegsY = 0;
};

struct BCTypes
{
	BCType alongX = BCType::None;
	BCType alongY = BCType::None;
};

struct BoundaryConditions
{
	std::vector<double> left;
	std::vector<double> right;
	std::vector<double> bottom;
	std::vector<double> top;
};

class PoissonSolver
{
public:
	PoissonSolver(UniformGrid grid, BCTypes bcTypes, BoundaryConditions bc);

	std::vector<double> Solve(std::vector<double> rhoUnfolded);

private:
	void Assemble();

	enum Flag { X, Y, Left, Right, Bottom, Top };

	int Pos(int m, int n, Flag flag, int cell = 0);

private:
	UniformGrid _grid;
	BCTypes _bcTypes;
	BoundaryConditions _bc;

	// Cells per region
	int _nx;
	int _ny;

	int	_nEq;
	SpMat _system;
	Eigen::VectorXd _rhs;
	Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> _solver;
};