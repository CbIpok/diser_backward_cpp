#ifndef MANAGERS_H
#define MANAGERS_H

#include <string>
#include <vector>
#include "stable_data_structs.h"
#include "hdf5.h"
#include <filesystem>
#include <stdexcept>
#include <regex>
#include <limits>
#include <future>
#include <algorithm>

namespace fs = std::filesystem;

// Шаблонная обёртка для представления 3D-массива в виде непрерывного блока памяти
template <typename T>
struct Array3DView {
    std::vector<T> data; // непрерывный буфер
    hsize_t T_dim, Y_dim, X_dim;

    Array3DView(hsize_t T, hsize_t Y, hsize_t X)
        : data(T* Y* X), T_dim(T), Y_dim(Y), X_dim(X) {
    }

    // Доступ к элементу (t, y, x)
    T& operator()(hsize_t t, hsize_t y, hsize_t x) {
        return data[t * Y_dim * X_dim + y * X_dim + x];
    }

    const T& operator()(hsize_t t, hsize_t y, hsize_t x) const {
        return data[t * Y_dim * X_dim + y * X_dim + x];
    }
};

// Функция для открытия HDF5-файла с проверкой ошибок.
int open_nc_file(const std::string& filename, hid_t& file);

// Функция чтения данных из HDF5-файла, возвращающая 3D view.
// Считываются данные из набора данных "height" для региона [y_start, y_end)
Array3DView<double> read_nc_file(const fs::path& filePath, int y_start, int y_end);

// Класс для работы с данными basis, содержащимися в NetCDF-файлах
class BasisManager {
public:
    std::string folder; // путь к каталогу с basis-файлами (NetCDF-файлы)

    explicit BasisManager(const std::string& folder_) : folder(folder_) {}

    // Функция чтения данных basis для региона [y_start, y_end)
    // Возвращает вектор 3D view: [num_files]{[T][region_height][X]}
    std::vector<Array3DView<double>> get_fk_region(int y_start, int y_end);
};

// Класс для работы с мариограммами (Wave data)
class WaveManager {
public:
    std::string nc_file; // путь к NetCDF-файлу с мариограммами

    explicit WaveManager(const std::string& nc_file_) : nc_file(nc_file_) {}

    // Функция загрузки данных переменной "height" для региона [y_start, y_end)
    // Возвращает 3D view: [T][region_height][X]
    Array3DView<double> load_mariogramm_by_region(int y_start, int y_end);
};

#endif // MANAGERS_H
