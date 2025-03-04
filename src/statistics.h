#ifndef STATISTICS_H
#define STATISTICS_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "stable_data_structs.h"
#include "managers.h"

// Определяем тип для хранения коэффициентов для каждого пикселя:
// двумерный массив (размер region_height x region_width) элементов, где каждый элемент – вектор коэффициентов.
using CoeffMatrix = std::vector<std::vector<Eigen::VectorXd>>;

// Функция для вычисления статистики аппроксимации по всему региону
void calculate_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config,
    CoeffMatrix& statistics_orto);

// Функция для сохранения статистики в CSV-файлы
void save_and_plot_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config);

#endif // STATISTICS_H
