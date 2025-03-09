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

// ��������� ������ ��� ������������� 3D-������� � ���� ������������ ����� ������
template <typename T>
struct Array3DView {
    std::vector<T> data; // ����������� �����
    hsize_t T_dim, Y_dim, X_dim;

    Array3DView(hsize_t T, hsize_t Y, hsize_t X)
        : data(T* Y* X), T_dim(T), Y_dim(Y), X_dim(X) {
    }

    // ������ � �������� (t, y, x)
    T& operator()(hsize_t t, hsize_t y, hsize_t x) {
        return data[t * Y_dim * X_dim + y * X_dim + x];
    }

    const T& operator()(hsize_t t, hsize_t y, hsize_t x) const {
        return data[t * Y_dim * X_dim + y * X_dim + x];
    }
};

// ������� ��� �������� HDF5-����� � ��������� ������.
int open_nc_file(const std::string& filename, hid_t& file);

// ������� ������ ������ �� HDF5-�����, ������������ 3D view.
// ����������� ������ �� ������ ������ "height" ��� ������� [y_start, y_end)
Array3DView<double> read_nc_file(const fs::path& filePath, int y_start, int y_end);

// ����� ��� ������ � ������� basis, ������������� � NetCDF-������
class BasisManager {
public:
    std::string folder; // ���� � �������� � basis-������� (NetCDF-�����)

    explicit BasisManager(const std::string& folder_) : folder(folder_) {}

    // ������� ������ ������ basis ��� ������� [y_start, y_end)
    // ���������� ������ 3D view: [num_files]{[T][region_height][X]}
    std::vector<Array3DView<double>> get_fk_region(int y_start, int y_end);
};

// ����� ��� ������ � ������������� (Wave data)
class WaveManager {
public:
    std::string nc_file; // ���� � NetCDF-����� � �������������

    explicit WaveManager(const std::string& nc_file_) : nc_file(nc_file_) {}

    // ������� �������� ������ ���������� "height" ��� ������� [y_start, y_end)
    // ���������� 3D view: [T][region_height][X]
    Array3DView<double> load_mariogramm_by_region(int y_start, int y_end);
};

#endif // MANAGERS_H
