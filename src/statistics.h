#ifndef STATISTICS_H
#define STATISTICS_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "stable_data_structs.h"
#include "managers.h"

// Новая структура для хранения коэффициентов и ошибки аппроксимации
struct CoefficientData {
    Eigen::VectorXd coefs;
    double aprox_error;
};

// Тип для хранения данных по всем пикселям: двумерный массив объектов CoefficientData
using CoeffMatrix = std::vector<std::vector<CoefficientData>>;

// Функция для вычисления статистики аппроксимации по всему региону
void calculate_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config,
    CoeffMatrix& statistics_orto);

// Функция для сохранения статистики в JSON-файл
void save_and_plot_statistics(const std::string& root_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config);

#endif // STATISTICS_H
