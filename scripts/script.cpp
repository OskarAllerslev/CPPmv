// script.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

// -----------------------------------------------------
// Læs CSV UDEN header (kun tal). Dynamisk opdagelse af
// antal rækker og kolonner.
// -----------------------------------------------------
Eigen::MatrixXd read_csv_no_header(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Fejl: Kunne ikke åbne filen " + filename);
    }

    // Først tæller vi antal rækker (rows) og kolonner (cols)
    int rows = 0, cols = 0;
    {
        std::string line;
        while (std::getline(file, line))
        {
            rows++;
            // Hvis det er første linje, find antallet af kolonner
            if (rows == 1)
            {
                std::stringstream ss(line);
                std::string cell;
                while (std::getline(ss, cell, ','))
                {
                    cols++;
                }
            }
        }
    }

    // Gå tilbage til filens begyndelse
    file.clear();
    file.seekg(0);

    // Opret en MatrixXd med de fundne dimensioner
    Eigen::MatrixXd matrix(rows, cols);

    // Læs data ind
    std::string line;
    int row_idx = 0;
    while (std::getline(file, line) && row_idx < rows)
    {
        std::stringstream ss(line);
        std::string cell;
        int col_idx = 0;

        while (std::getline(ss, cell, ',') && col_idx < cols)
        {
            matrix(row_idx, col_idx) = std::stod(cell);
            col_idx++;
        }
        row_idx++;
    }

    file.close();
    return matrix;
}

// -----------------------------------------------------
// Mean-variance optimering uden short selling
// min w^T Σ w, s.t. sum(w)=1, w >= 0
// -----------------------------------------------------
Eigen::VectorXd mean_variance_optimization(const Eigen::MatrixXd &cov_matrix,
                                           const Eigen::VectorXd &expected_returns)
{
    int n = expected_returns.size();

    // Vektor af 1'ere
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);

    // Udvidet system:
    // [ Σ    1 ] [ w ] = [ -E(R) ]
    // [ 1^T  0 ] [ λ ]   [   1    ]

    Eigen::MatrixXd A(n + 1, n + 1);
    A.block(0, 0, n, n) = cov_matrix;
    A.block(0, n, n, 1) = ones;
    A.block(n, 0, 1, n) = ones.transpose();
    A(n, n) = 0;

    Eigen::VectorXd b(n + 1);
    b.head(n) = -expected_returns; // Negativ for at minimere variansen
    b(n) = 1;

    // Løs systemet
    Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);

    // Ekstraher vægte
    Eigen::VectorXd w = solution.head(n);

    // Skær negative vægte fra (ingen short selling)
    w = w.cwiseMax(0.0);

    // Normaliser, så sum(w) = 1
    double sum_w = w.sum();
    if (sum_w > 1e-12)
        w /= sum_w;
    else
        w = Eigen::VectorXd::Zero(n);

    return w;
}

int main()
{
    try
    {
        // Læs data UDEN header
        Eigen::MatrixXd cov_matrix = read_csv_no_header("cov_matrix.csv");
        Eigen::MatrixXd exp = read_csv_no_header("expected_returns.csv");

        // Forventede afkast CSV bliver en matrix med dimension (rows x cols).
        // Vi forventer at have (N x 1). Tjek:
        if (exp.cols() != 1)
        {
            throw std::runtime_error("Fejl: expected_returns.csv har ikke præcis 1 kolonne.");
        }

        // Konverter til VectorXd
        Eigen::VectorXd expected_returns = exp.col(0);

        std::cout << "Kovariansmatrix har dimension: "
                  << cov_matrix.rows() << " x " << cov_matrix.cols() << std::endl;
        std::cout << "Forventede afkast har dimension: "
                  << expected_returns.size() << std::endl;

        if (cov_matrix.rows() != cov_matrix.cols())
        {
            throw std::runtime_error("Fejl: Kovariansmatrix er ikke kvadratisk!");
        }
        if (cov_matrix.rows() != expected_returns.size())
        {
            throw std::runtime_error("Fejl: Dimensionen af kovariansmatrix passer ikke til expected_returns!");
        }

        // Kør mean-variance optimering
        Eigen::VectorXd optimal_weights = mean_variance_optimization(cov_matrix, expected_returns);

        std::cout << "Optimal portefølje vægte:\n" << optimal_weights << std::endl;

        // Gem vægte i 'optimal_weights.csv' UDEN header
        std::ofstream wfile("optimal_weights.csv");
        if (!wfile.is_open())
        {
            throw std::runtime_error("Fejl: kunne ikke skrive til optimal_weights.csv");
        }
        for (int i = 0; i < optimal_weights.size(); i++)
        {
            wfile << optimal_weights[i];
            if (i < optimal_weights.size() - 1)
                wfile << ",";
            else
                wfile << "\n";
        }
        wfile.close();
        std::cout << "Optimal vægte gemt i optimal_weights.csv" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Undtagelse: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
