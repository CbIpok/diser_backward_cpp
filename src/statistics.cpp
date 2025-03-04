#include "statistics.h"
#include "approx_orto.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <future>
#include <algorithm>
#include "json.hpp"  // Подключение библиотеки nlohmann::json

using json = nlohmann::json;

//// Реализация неортогональной аппроксимации с использованием QR-разложения
//std::pair<Eigen::VectorXd, Eigen::VectorXd> approximate_with_non_orthogonal_basis(const Eigen::VectorXd& x, const Eigen::MatrixXd& basis) {
//    // Решаем систему методом наименьших квадратов:
//    Eigen::VectorXd coeffs = basis.transpose().colPivHouseholderQr().solve(x);
//    // Вычисляем аппроксимацию (не используется далее)
//    Eigen::VectorXd approximation = basis.transpose() * coeffs;
//    return { approximation, coeffs };
//}

//
// Функция calculate_statistics
// Загружает данные из NetCDF (WaveManager и BasisManager), затем для каждого пикселя региона
// извлекает вектор времени (wave_vector) и соответствующий базис (smoothed_basis) и вычисляет
// коэффициенты аппроксимации обычным методом (non orto) и методом с ортогонализацией (orto).
//
void calculate_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config,
    CoeffMatrix& statistics_orto) {
    // Формирование путей
    std::string basis_path = root_folder + "/" + bath + "/" + basis;
    std::string wave_nc_path = root_folder + "/" + bath + "/" + wave + ".nc";

    BasisManager basis_manager(basis_path);
    WaveManager wave_manager(wave_nc_path);

    int width = area_config.all[0];
    int height = area_config.all[1];
    int batch_size = 8;
    int y_start_init = 75;

    statistics_orto.clear();

    // Для каждого пачки по оси y
    for (int y_start = y_start_init; y_start < height / 4; y_start += batch_size) {
        int y_end = std::min(y_start + batch_size, height);
        // Загружаем данные
        auto wave_data = wave_manager.load_mariogramm_by_region(y_start, y_end);
        auto fk_data = basis_manager.get_fk_region(y_start, y_end);
        if (wave_data.empty() || fk_data.empty()) continue;
        // Предполагаем, что wave_data имеет размер [T][region_height][region_width]
        int T = wave_data.size();
        int region_height = wave_data[0].size();
        int region_width = wave_data[0][0].size();
        int n_basis = fk_data.size(); // количество basis файлов

        // В Python использовался диапазон x от 0 до width/3
        int x_max = width / 4;
        std::cout << "loaded\n";
        // Create a vector to hold futures for each row.
        std::vector<std::future<std::vector<Eigen::VectorXd>>> futures;
        futures.reserve(region_height);

        // Launch an async task for each row.
        for (int i = 0; i < region_height; i++) {
            // Capture i by value; capture wave_data and fk_data by reference.
            futures.push_back(std::async(std::launch::async, [i, T, x_max, n_basis, &wave_data, &fk_data]() -> std::vector<Eigen::VectorXd> {
                std::vector<Eigen::VectorXd> row_orto;
                // Process each x in the row.
                for (int x = 0; x < x_max; x++) {
                    // Build wave_vector from wave_data[t][i][x] for each t.
                    Eigen::VectorXd wave_vector(T);
                    for (int t = 0; t < T; t++) {
                        wave_vector[t] = wave_data[t][i][x];
                    }
                    // Build the smoothed_basis matrix (n_basis x T) from fk_data[b][t][i][x].
                    Eigen::MatrixXd smoothed_basis(n_basis, T);
                    for (int b = 0; b < n_basis; b++) {
                        for (int t = 0; t < T; t++) {
                            smoothed_basis(b, t) = fk_data[b][t][i][x];
                        }
                    }
                    // Ensure dimensions match.
                    if (smoothed_basis.cols() != wave_vector.size()) continue;
                    // Compute the orthogonal approximation coefficients.
                    Eigen::VectorXd coefs_orto = approximate_with_non_orthogonal_basis_orto(wave_vector, smoothed_basis);
                    row_orto.push_back(coefs_orto);
                }
                return row_orto;
                }));
        }

        // Retrieve the results in order and push them into statistics_orto.
        for (auto& future : futures) {
            auto row_orto = future.get();
            if (!row_orto.empty()) {
                statistics_orto.push_back(row_orto);
            }
        }

    }
}

// Функция сохранения двумерного массива коэффициентов в CSV-файл
void save_coefficients_csv(const std::string& filename, const CoeffMatrix& coeffs) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << " для записи.\n";
        return;
    }
    // Каждая строка – один ряд статистики, каждая ячейка содержит вектор коэффициентов (числа разделены пробелами)
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
    std::cout << "Сохранено: " << filename << "\n";
}

void save_coefficients_json(const std::string& filename, const CoeffMatrix& coeffs) {
    nlohmann::json j;
    // Проходим по всем строкам и столбцам матрицы коэффициентов
    for (size_t row = 0; row < coeffs.size(); ++row) {
        for (size_t col = 0; col < coeffs[row].size(); ++col) {
            // Формирование ключа вида "[row,col]"
            std::string key = "[" + std::to_string(row) + "," + std::to_string(col) + "]";
            // Преобразуем Eigen::VectorXd в std::vector<double>
            std::vector<double> vec(coeffs[row][col].data(), coeffs[row][col].data() + coeffs[row][col].size());
            j[key] = vec;
        }
    }

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << " для записи.\n";
        return;
    }

    // Запись в файл с отступами для читаемости
    ofs << j.dump(4);
    ofs.close();
    std::cout << "Сохранено: " << filename << "\n";
}

// Функция save_and_plot_statistics: вычисляет статистику и сохраняет коэффициенты в CSV
void save_and_plot_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config) {
    CoeffMatrix statistics_orto;
    calculate_statistics(root_folder, bath, wave, basis, area_config, statistics_orto);

    std::string filename_orto = "case_statistics_hd_y_" + basis + bath + "_o.json";
    //std::string filename_non_orto = "case_statistics_hd_y_" + basis + "_no.csv";

    save_coefficients_json(filename_orto, statistics_orto);
    /*save_coefficients_json(filename_non_orto, statistics_non_orto);*/

    // Визуализация не реализована – коэффициенты можно открыть в Excel или передать в Python для построения графиков.
}
