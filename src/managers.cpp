#include "managers.h"
#include <netcdf.h>
#include <iostream>
<<<<<<< Updated upstream
#include <filesystem>
#include <algorithm>
#include <vector>
#include <future>
#include <regex>
namespace fs = std::filesystem;

// Функция для открытия NetCDF-файла с проверкой ошибок
int open_nc_file(const std::string& filename, int& ncid) {
    int retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Ошибка открытия файла " << filename << " : " << nc_strerror(retval) << std::endl;
=======
#include <adios2.h>
//
// Функция для открытия HDF5-файла с проверкой ошибок.
//
int open_nc_file(const std::string& filename, hid_t& file) {
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "Ошибка открытия файла " << filename << std::endl;
        return -1;
>>>>>>> Stashed changes
    }
    return retval;
}
std::vector<std::vector<std::vector<double>>> read_nc_file(const fs::path& file, int y_start, int y_end) {
    std::vector<std::vector<std::vector<double>>> data;
    std::cout << "loading: " << file << std::endl;
    int ncid;

<<<<<<< Updated upstream
    if (open_nc_file(file.string(), ncid) != NC_NOERR)
        return data;

    int varid;
    int retval = nc_inq_varid(ncid, "height", &varid);
    if (retval != NC_NOERR) {
        std::cerr << "Переменная 'height' не найдена в " << file.string() << std::endl;
        nc_close(ncid);
        return data;
    }

    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);
    if (ndims != 3) {
        std::cerr << "Ожидалось 3 измерения в файле " << file.string() << std::endl;
        nc_close(ncid);
        return data;
=======
//
// Изменённая функция чтения данных из HDF5-файла: данные считываются напрямую
// в непрерывный буфер, представленный Array3DView, без промежуточного копирования.
//
Array3DView<double> read_nc_file(const fs::path& filePath, int y_start, int y_end)
{
    std::cout << "loading: " << filePath << std::endl;

    // Инициализация ADIOS2 (используем настройки по умолчанию)
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("ReadIO");
    io.SetEngine("BPFile");
    // Открываем файл в режиме чтения
    adios2::Engine reader = io.Open(filePath.string(), adios2::Mode::Read);
    if (!reader)
    {
        throw std::runtime_error("Ошибка открытия файла " + filePath.string());
    }

    // Начинаем шаг чтения (для файлового режима BeginStep сразу после открытия)
    reader.BeginStep();

    // Получаем переменную "height"
    adios2::Variable<double> var = io.InquireVariable<double>("height");
    if (!var)
    {
        throw std::runtime_error("Переменная 'height' не найдена в файле " + filePath.string());
    }

    // Получаем размеры глобального массива (ожидается 3 измерения: [T, Y, X])
    std::vector<size_t> shape = var.Shape();
    if (shape.size() != 3)
    {
        throw std::runtime_error("Ожидалось 3 измерения в файле " + filePath.string());
>>>>>>> Stashed changes
    }
    size_t T = shape[0];
    size_t Y = shape[1];
    size_t X = shape[2];

<<<<<<< Updated upstream
    int dimids[3];
    size_t T, Y, X;
    nc_inq_vardimid(ncid, varid, dimids);
    nc_inq_dimlen(ncid, dimids[0], &T);
    nc_inq_dimlen(ncid, dimids[1], &Y);
    nc_inq_dimlen(ncid, dimids[2], &X);

    int local_y_end = y_end;
    if (static_cast<size_t>(local_y_end) > Y)
        local_y_end = static_cast<int>(Y);
    size_t region_height = local_y_end - y_start;
    data.resize(T, std::vector<std::vector<double>>(region_height, std::vector<double>(X, 0.0)));

    size_t start[3] = { 0, static_cast<size_t>(y_start), 0 };
    size_t count[3] = { T, region_height, X };
    std::vector<double> buffer(T * region_height * X, 0.0);
    std::cout << "start\n";
    retval = nc_get_vara_double(ncid, varid, start, count, buffer.data());
    std::cout << "end\n";
    if (retval != NC_NOERR) {
        std::cerr << "Ошибка чтения файла " << file.string() << " : " << nc_strerror(retval) << std::endl;
        nc_close(ncid);
        return data;
    }
    nc_close(ncid);

    // Копирование данных из буфера в 3D-вектор
    for (size_t t = 0; t < T; t++) {
        for (size_t i = 0; i < region_height; i++) {
            for (size_t x = 0; x < X; x++) {
                size_t idx = t * region_height * X + i * X + x;
                data[t][i][x] = buffer[idx];
            }
        }
    }
    return data;
=======
    // Корректировка y_end, если выходит за пределы данных
    int local_y_end = y_end > static_cast<int>(Y) ? static_cast<int>(Y) : y_end;
    size_t region_height = local_y_end - y_start;

    // Определяем область (selection) для чтения: [0, y_start, 0] размером [T, region_height, X]
    std::vector<size_t> start{ 0, static_cast<size_t>(y_start), 0 };
    std::vector<size_t> count{ T, region_height, X };
    var.SetSelection({ start, count });

    // Создаем объект для хранения прочитанных данных
    Array3DView<double> view(T, region_height, X);

    // Запрос на чтение данных в буфер view.data
    reader.Get<double>(var, view.data.data());
    reader.EndStep();
    reader.Close();

    return view;
>>>>>>> Stashed changes
}

// Реализация метода WaveManager::load_mariogramm_by_region с использованием netcdf.h
std::vector<std::vector<std::vector<double>>> WaveManager::load_mariogramm_by_region(int y_start, int y_end) {
    
    return read_nc_file(nc_file,y_start,y_end);
}


// Функция для извлечения индекса из имени файла
int extractIndex(const fs::path& filePath) {
    // Регулярное выражение для поиска шаблона _<число>.nc
    std::regex regexPattern("_(\\d+)\\.nc");
    std::smatch match;
    std::string filename = filePath.filename().string();
    if (std::regex_search(filename, match, regexPattern)) {
        return std::stoi(match[1].str());
    }
    // Если индекс не найден, возвращаем максимальное значение,
    // чтобы файл оказался в конце отсортированного списка.
    return std::numeric_limits<int>::max();
}

// Функция для получения отсортированного списка файлов
<<<<<<< Updated upstream
std::vector<fs::path> getSortedFileList(const std::string& folder) {
    std::vector<fs::path> files;

    // Перебираем все файлы в заданной директории
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.is_regular_file()) {
            fs::path filePath = entry.path();
            // Проверяем, что файл имеет расширение ".nc" и содержит символ '_'
            if (filePath.extension() == ".nc" &&
                filePath.filename().string().find('_') != std::string::npos) {
                files.push_back(filePath);
            }
        }
    }

    // Сортируем файлы по числовому значению, извлечённому из имени файла
    std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
        return extractIndex(a) < extractIndex(b);
        });

    return files;
}


std::vector<std::vector<std::vector<std::vector<double>>>> BasisManager::get_fk_region(int y_start, int y_end) {
    std::vector<std::vector<std::vector<std::vector<double>>>> fk;
    std::vector<fs::path> files = getSortedFileList(folder);
=======
//
std::vector<fs::path> getSortedFolderList(const std::string& folder) {
    std::vector<fs::path> directories;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.is_directory()) {
            fs::path dirPath = entry.path();
            if (dirPath.extension() == ".bp" &&
                dirPath.filename().string().find('_') != std::string::npos) {
                directories.push_back(dirPath);
            }
        }
    }
    std::sort(directories.begin(), directories.end(), [](const fs::path& a, const fs::path& b) {
        return extractIndex(a) < extractIndex(b);
        });
    return directories;
}

//
// Реализация метода BasisManager::get_fk_region с использованием Array3DView.
// Для каждого файла из каталога создаётся 3D view, возвращаемый в векторе.
//
std::vector<Array3DView<double>> BasisManager::get_fk_region(int y_start, int y_end) {
    std::vector<Array3DView<double>> fk;
    std::vector<fs::path> files = getSortedFolderList(folder);
>>>>>>> Stashed changes

    // Последовательная обработка файлов
    for (const auto& file : files) {
        auto file_data = read_nc_file(file, y_start, y_end);
        if (!file_data.empty()) {
            fk.push_back(file_data);
        }
    }
    return fk;
}