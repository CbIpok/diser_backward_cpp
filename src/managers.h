#ifndef MANAGERS_H
#define MANAGERS_H

#include <string>
#include <vector>
#include "stable_data_structs.h"

// Класс для работы с данными basis, содержащимися в NetCDF-файлах
class BasisManager {
public:
    std::string folder; // путь к каталогу с basis-файлами (NetCDF-файлы)

    explicit BasisManager(const std::string& folder_) : folder(folder_) {}

    // Функция чтения данных basis для региона [y_start, y_end)
    // Возвращает 4D массив: [num_files][T][region_height][X]
    std::vector<std::vector<std::vector<std::vector<double>>>> get_fk_region(int y_start, int y_end);
};

// Класс для работы с мариограммами (Wave data)
class WaveManager {
public:
    std::string nc_file; // путь к NetCDF-файлу с мариограммами

    explicit WaveManager(const std::string& nc_file_) : nc_file(nc_file_) {}

    // Функция загрузки данных переменной "height" для региона [y_start, y_end)
    // Возвращает 3D массив: [T][region_height][X]
    std::vector<std::vector<std::vector<double>>> load_mariogramm_by_region(int y_start, int y_end);
};

#endif // MANAGERS_H
