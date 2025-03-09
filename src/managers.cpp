#include "managers.h"
#include <iostream>

//
// ������� ��� �������� HDF5-����� � ��������� ������.
//
int open_nc_file(const std::string& filename, hid_t& file) {
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "������ �������� ����� " << filename << std::endl;
        return -1;
    }
    return 0;
}

//
// ��������� ������� ������ ������ �� HDF5-�����: ������ ����������� ��������
// � ����������� �����, �������������� Array3DView, ��� �������������� �����������.
//
Array3DView<double> read_nc_file(const fs::path& filePath, int y_start, int y_end) {
    std::cout << "loading: " << filePath << std::endl;

    hid_t file;
    if (open_nc_file(filePath.string(), file) != 0)
        throw std::runtime_error("������ �������� �����");

    // ��������� ����� ������ "height"
    hid_t dataset = H5Dopen(file, "height", H5P_DEFAULT);
    if (dataset < 0) {
        std::cerr << "������ �������� ������ ������ 'height' � ����� " << filePath.string() << std::endl;
        H5Fclose(file);
        throw std::runtime_error("������ �������� ������ ������");
    }

    // �������� dataspace ������ ������
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0) {
        std::cerr << "������ ��������� dataspace ������ ������ 'height' � ����� " << filePath.string() << std::endl;
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("������ ��������� dataspace");
    }

    // ��������� ����� ���������
    int ndims = H5Sget_simple_extent_ndims(dataspace);
    if (ndims != 3) {
        std::cerr << "��������� 3 ��������� � ����� " << filePath.string() << std::endl;
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("�������� ����� ���������");
    }

    // �������� ������� ���������
    hsize_t dims[3];
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);
    hsize_t T = dims[0];
    hsize_t Y = dims[1];
    hsize_t X = dims[2];

    // ������������� y_end, ���� �� ������� �� ������� ������
    int local_y_end = y_end > static_cast<int>(Y) ? static_cast<int>(Y) : y_end;
    hsize_t region_height = local_y_end - y_start;

    // ���������� hyperslab � �������� dataspace
    hsize_t offset[3] = { 0, static_cast<hsize_t>(y_start), 0 };
    hsize_t count[3] = { T, region_height, X };
    herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    if (status < 0) {
        std::cerr << "������ ������ hyperslab � ����� " << filePath.string() << std::endl;
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("������ ������ hyperslab");
    }

    // ������ dataspace ��� ������ � ���� �� ���������
    hid_t memspace = H5Screate_simple(3, count, nullptr);
    if (memspace < 0) {
        std::cerr << "������ �������� memory dataspace ��� hyperslab � ����� " << filePath.string() << std::endl;
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("������ �������� memory dataspace");
    }

    // ������ Array3DView ��� �������� ������ ��� ��������������� �����������
    Array3DView<double> view(T, region_height, X);

    // ��������� ������ ��������������� � ����������� ����� view.data
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, view.data.data());
    if (status < 0) {
        std::cerr << "������ ������ ������ �� ����� " << filePath << std::endl;
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("������ ������ ������");
    }

    // ����������� ������� HDF5
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);

    return view;
}

//
// ���������� ������ WaveManager::load_mariogramm_by_region � �������������� Array3DView.
//
Array3DView<double> WaveManager::load_mariogramm_by_region(int y_start, int y_end) {
    return read_nc_file(nc_file, y_start, y_end);
}

//
// ������� ��� ���������� ������� �� ����� �����
//
int extractIndex(const fs::path& filePath) {
    std::regex regexPattern("_(\\d+)\\.nc");
    std::smatch match;
    std::string filename = filePath.filename().string();
    if (std::regex_search(filename, match, regexPattern)) {
        return std::stoi(match[1].str());
    }
    return std::numeric_limits<int>::max();
}

//
// ������� ��� ��������� ���������������� ������ ������
//
std::vector<fs::path> getSortedFileList(const std::string& folder) {
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.is_regular_file()) {
            fs::path filePath = entry.path();
            if (filePath.extension() == ".nc" &&
                filePath.filename().string().find('_') != std::string::npos) {
                files.push_back(filePath);
            }
        }
    }
    std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
        return extractIndex(a) < extractIndex(b);
        });
    return files;
}

//
// ���������� ������ BasisManager::get_fk_region � �������������� Array3DView.
// ��� ������� ����� �� �������� �������� 3D view, ������������ � �������.
//
std::vector<Array3DView<double>> BasisManager::get_fk_region(int y_start, int y_end) {
    std::vector<Array3DView<double>> fk;
    std::vector<fs::path> files = getSortedFileList(folder);

    // ����������� ������ ��� ������� �����
    std::vector<std::future<Array3DView<double>>> futures;
    for (const auto& file : files) {
        futures.push_back(std::async(std::launch::async, read_nc_file, file, y_start, y_end));
    }

    // �������� ����������
    for (auto& fut : futures) {
        fk.push_back(fut.get());
    }
    return fk;
}
