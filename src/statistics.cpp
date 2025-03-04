#include "statistics.h"
#include "approx_orto.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "json.hpp"  // ����������� ���������� nlohmann::json

using json = nlohmann::json;

//// ���������� ��������������� ������������� � �������������� QR-����������
//std::pair<Eigen::VectorXd, Eigen::VectorXd> approximate_with_non_orthogonal_basis(const Eigen::VectorXd& x, const Eigen::MatrixXd& basis) {
//    // ������ ������� ������� ���������� ���������:
//    Eigen::VectorXd coeffs = basis.transpose().colPivHouseholderQr().solve(x);
//    // ��������� ������������� (�� ������������ �����)
//    Eigen::VectorXd approximation = basis.transpose() * coeffs;
//    return { approximation, coeffs };
//}

//
// ������� calculate_statistics
// ��������� ������ �� NetCDF (WaveManager � BasisManager), ����� ��� ������� ������� �������
// ��������� ������ ������� (wave_vector) � ��������������� ����� (smoothed_basis) � ���������
// ������������ ������������� ������� ������� (non orto) � ������� � ���������������� (orto).
//
void calculate_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config,
    CoeffMatrix& statistics_orto) {
    // ������������ �����
    std::string basis_path = root_folder + "/" + bath + "/" + basis;
    std::string wave_nc_path = root_folder + "/" + bath + "/" + wave + ".nc";

    BasisManager basis_manager(basis_path);
    WaveManager wave_manager(wave_nc_path);

    int width = area_config.all[0];
    int height = area_config.all[1];
    int batch_size = 48;
    int y_start_init = 75;

    statistics_orto.clear();

    // ��� ������� ����� �� ��� y
    for (int y_start = y_start_init; y_start < height / 4; y_start += batch_size) {
        int y_end = std::min(y_start + batch_size, height);
        // ��������� ������
        auto wave_data = wave_manager.load_mariogramm_by_region(y_start, y_end);
        auto fk_data = basis_manager.get_fk_region(y_start, y_end);
        if (wave_data.empty() || fk_data.empty()) continue;
        // ������������, ��� wave_data ����� ������ [T][region_height][region_width]
        int T = wave_data.size();
        int region_height = wave_data[0].size();
        int region_width = wave_data[0][0].size();
        int n_basis = fk_data.size(); // ���������� basis ������

        // � Python ������������� �������� x �� 0 �� width/3
        int x_max = width / 4;
        std::cout << "loaded\n";
        // ��� ������ ������ �������
        for (int i = 0; i < region_height; i++) {
            std::cout << "y = " << i << "\n";
            std::vector<Eigen::VectorXd> row_orto;
            std::vector<Eigen::VectorXd> row_non_orto;
            for (int x = 0; x < x_max; x++) {
                //std::cout << "x = " << x << "\n";
                // ��������� wave_vector: ��� ������� t ����� wave_data[t][i][x]
                Eigen::VectorXd wave_vector(T);
                for (int t = 0; t < T; t++) {
                    wave_vector[t] = wave_data[t][i][x];
                }
                // ��������� smoothed_basis: ������ ������� �������� (n_basis x T)
                Eigen::MatrixXd smoothed_basis(n_basis, T);
                for (int b = 0; b < n_basis; b++) {
                    for (int t = 0; t < T; t++) {
                        smoothed_basis(b, t) = fk_data[b][t][i][x];
                    }
                }
                if (smoothed_basis.cols() != wave_vector.size()) continue;
                // ��������� ������������ ������� ������������� (least squares)
                //auto [approx_non_orto, coefs_non_orto] = approximate_with_non_orthogonal_basis(wave_vector, smoothed_basis);
                // ��������� ������������ ������������� ������������� (�� approx_orto.cpp)
                Eigen::VectorXd coefs_orto = approximate_with_non_orthogonal_basis_orto(wave_vector, smoothed_basis);
                //row_non_orto.push_back(coefs_non_orto);
                row_orto.push_back(coefs_orto);
            }
            if (!row_orto.empty()) {
                statistics_orto.push_back(row_orto);
                //statistics_non_orto.push_back(row_non_orto);
            }
            /*if (i == 1)
                return;*/
        }

    }
}

// ������� ���������� ���������� ������� ������������� � CSV-����
void save_coefficients_csv(const std::string& filename, const CoeffMatrix& coeffs) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "�� ������� ������� ���� " << filename << " ��� ������.\n";
        return;
    }
    // ������ ������ � ���� ��� ����������, ������ ������ �������� ������ ������������� (����� ��������� ���������)
    for (const auto& row : coeffs) {
        bool firstCell = true;
        for (const auto& vec : row) {
            if (!firstCell) ofs << ",";
            firstCell = false;
            std::ostringstream oss;
            for (int i = 0; i < vec.size(); i++) {
                oss << vec[i];
                if (i + 1 < vec.size()) oss << " ";
            }
            ofs << oss.str();
        }
        ofs << "\n";
    }
    ofs.close();
    std::cout << "���������: " << filename << "\n";
}

void save_coefficients_json(const std::string& filename, const CoeffMatrix& coeffs) {
    nlohmann::json j;
    // �������� �� ���� ������� � �������� ������� �������������
    for (size_t row = 0; row < coeffs.size(); ++row) {
        for (size_t col = 0; col < coeffs[row].size(); ++col) {
            // ������������ ����� ���� "[row,col]"
            std::string key = "[" + std::to_string(row) + "," + std::to_string(col) + "]";
            // ����������� Eigen::VectorXd � std::vector<double>
            std::vector<double> vec(coeffs[row][col].data(), coeffs[row][col].data() + coeffs[row][col].size());
            j[key] = vec;
        }
    }

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "�� ������� ������� ���� " << filename << " ��� ������.\n";
        return;
    }

    // ������ � ���� � ��������� ��� ����������
    ofs << j.dump(4);
    ofs.close();
    std::cout << "���������: " << filename << "\n";
}

// ������� save_and_plot_statistics: ��������� ���������� � ��������� ������������ � CSV
void save_and_plot_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config) {
    CoeffMatrix statistics_orto;
    calculate_statistics(root_folder, bath, wave, basis, area_config, statistics_orto);

    std::string filename_orto = "case_statistics_hd_y_" + basis + "_o.json";
    //std::string filename_non_orto = "case_statistics_hd_y_" + basis + "_no.csv";

    save_coefficients_json(filename_orto, statistics_orto);
    /*save_coefficients_json(filename_non_orto, statistics_non_orto);*/

    // ������������ �� ����������� � ������������ ����� ������� � Excel ��� �������� � Python ��� ���������� ��������.
}
