#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>

//funk til at læse csv
Eigen::MatrixXd read_csv(const std::string& filename, int rows, int cols)
{
    std::ifstream file(filename);
    Eigen::MatrixXd matrix(rows,cols);
    std::string line;
    int row = 0;

    while (std::getline(file,line))
    {
        std::stringstream ss(line);
        std::string value;
        int col = 0;
        while (std::getline(ss, value, ','))
        {
            matrix(row, col) = std::stod(value);
            col++;
        }
        row++;
    }
    reutnr matrix;
}

Eigen::MatrixXd mean_variance_optimization(const Eigen::MatrixXd& cov_matrix, const Eigen::VectorXd& expected_returns)
{
    int n = expected_returns.size();
    Eigen::MatrixXd A = cov_matrix;
    Eigen::MatrixXd b = expected_returns;

    // løsning: min w'\sigma w, u.b.b w'e = 1, w >= 0 
    Eigen::VectorXd weights = A.colPivHouseHOlderQr().solve(b);

    weights = weights.cwiseMax(0);
    weights /= weights.sum();

    return weights;
}

int main()
{
    // dette skal ændres så det findes automatisk
    int num_assets = 4;

    Eigen::MatrixXd cov_matrix = read_csv("cov_matrix.csv", num_assets, num_assets);
    Eigen::VectorXd expected_returns = read_csv("expected_returns", num_assets, 1);

    Eigen::VectorXd optimal_weights = mean_variance_optimization(cov_matrix, expected_returns);

    std::cout << "Optimal portefølje vægte: " << std::endl; << optimal_weights << std::endl;
    return 0;
}


