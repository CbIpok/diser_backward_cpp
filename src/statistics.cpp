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

int count_from_name(const std::string& name) {
    std::size_t underscorePos = name.find('_');
    if (underscorePos != std::string::npos) {
        std::string numberPart = name.substr(underscorePos + 1);
        return std::stoi(numberPart);
    }
    return 0;
}

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
    int gigabyte_size = 116;
    int memory_in_gb = 16;
    int batch_size = gigabyte_size*memory_in_gb / count_from_name(basis);
    int y_start_init = 75;

    statistics_orto.clear();

    for (int y_start = y_start_init; y_start < height / 4; y_start += batch_size) {
        auto start = std::chrono::high_resolution_clock::now();
        int y_end = std::min(y_start + batch_size, height);
        // Получаем 3D view для волновых данных
        auto wave_data = wave_manager.load_mariogramm_by_region(y_start, y_end);
        // Получаем вектор 3D view для basis-данных
        auto fk_data = basis_manager.get_fk_region(y_start, y_end);
        // Фиксируем конечное время
        auto end = std::chrono::high_resolution_clock::now();

        // Вычисляем продолжительность выполнения в миллисекундах
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "time: " << duration << " мс" << std::endl;
        // Если данные не получены, переходим к следующему блоку
        if (wave_data.data.empty() || fk_data.empty()) continue;

        int T = wave_data.T_dim;
        int region_height = wave_data.Y_dim;
        int region_width = wave_data.X_dim;
        int n_basis = static_cast<int>(fk_data.size());

        int x_max = width / 4;
        std::cout << "loaded\n";

        std::vector<std::future<std::vector<CoefficientData>>> futures;
        futures.reserve(region_height);

        for (int i = 0; i < region_height; i++) {
            futures.push_back(std::async(std::launch::async, [i, T, x_max, n_basis, &wave_data, &fk_data]() -> std::vector<CoefficientData> {
                std::vector<CoefficientData> row_data;
                for (int x = 0; x < x_max; x++) {
                    // Формируем вектор значений сигнала для текущего пикселя
                    Eigen::VectorXd wave_vector(T);
                    for (int t = 0; t < T; t++) {
                        wave_vector[t] = wave_data(t, i, x);
                    }
                    // Формируем матрицу базиса (n_basis x T)
                    Eigen::MatrixXd smoothed_basis(n_basis, T);
                    for (int b = 0; b < n_basis; b++) {
                        for (int t = 0; t < T; t++) {
                            smoothed_basis(b, t) = fk_data[b](t, i, x);
                        }
                    }
                    if (smoothed_basis.cols() != wave_vector.size()) continue;

                    // Вычисление коэффициентов аппроксимации методом с ортогонализацией
                    Eigen::VectorXd coefs_orto = approximate_with_non_orthogonal_basis_orto(wave_vector, smoothed_basis);

                    // Вычисление аппроксимированного сигнала
                    Eigen::VectorXd approximation = smoothed_basis.transpose() * coefs_orto;
                    // Вычисление среднеквадратичной ошибки (RMSE)
                    double error = std::sqrt((wave_vector - approximation).squaredNorm() / wave_vector.size());

                    CoefficientData pixelData;
                    pixelData.coefs = coefs_orto;
                    pixelData.aprox_error = error;
                    row_data.push_back(pixelData);
                }
                return row_data;
                }));
        }

        for (auto& future : futures) {
            auto row_data = future.get();
            if (!row_data.empty()) {
                statistics_orto.push_back(row_data);
            }
        }
    }
}

void save_coefficients_json(const std::string& filename, const CoeffMatrix& coeffs) {
    nlohmann::json j;
    for (size_t row = 0; row < coeffs.size(); ++row) {
        for (size_t col = 0; col < coeffs[row].size(); ++col) {
            std::string key = "[" + std::to_string(row) + "," + std::to_string(col) + "]";
            std::vector<double> vec(coeffs[row][col].coefs.data(),
                coeffs[row][col].coefs.data() + coeffs[row][col].coefs.size());
            double error = coeffs[row][col].aprox_error;
            j[key] = { {"coefs", vec}, {"aprox_error", error} };
        }
    }
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << " для записи.\n";
        return;
    }
    ofs << j.dump(4);
    ofs.close();
    std::cout << "Сохранено: " << filename << "\n";
}

void save_and_plot_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config) {
    CoeffMatrix statistics_orto;
    calculate_statistics(root_folder, bath, wave, basis, area_config, statistics_orto);

    std::string filename_orto = "case_statistics_hd_y_" + basis + bath + "_o.json";

    save_coefficients_json(filename_orto, statistics_orto);
    // Визуализация не реализована – коэффициенты можно открыть в Excel или передать в Python для построения графиков.
}
