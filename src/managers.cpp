#include "managers.h"
#include <netcdf.h>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <vector>
#include <future>
namespace fs = std::filesystem;

// Функция для открытия NetCDF-файла с проверкой ошибок
int open_nc_file(const std::string& filename, int& ncid) {
    int retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Ошибка открытия файла " << filename << " : " << nc_strerror(retval) << std::endl;
    }
    return retval;
}

// Реализация метода WaveManager::load_mariogramm_by_region с использованием netcdf.h
std::vector<std::vector<std::vector<double>>> WaveManager::load_mariogramm_by_region(int y_start, int y_end) {
    std::vector<std::vector<std::vector<double>>> data;
    int ncid;
    if (open_nc_file(nc_file, ncid) != NC_NOERR) return data;

    int varid;
    int retval = nc_inq_varid(ncid, "height", &varid);
    if (retval != NC_NOERR) {
        std::cerr << "Переменная 'height' не найдена в файле " << nc_file << std::endl;
        nc_close(ncid);
        return data;
    }

    // Проверяем, что переменная имеет 3 измерения
    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);
    if (ndims != 3) {
        std::cerr << "Ожидалось 3 измерения, получено " << ndims << std::endl;
        nc_close(ncid);
        return data;
    }

    // Получаем размеры измерений
    int dimids[3];
    nc_inq_vardimid(ncid, varid, dimids);
    size_t T, Y, X;
    nc_inq_dimlen(ncid, dimids[0], &T);
    nc_inq_dimlen(ncid, dimids[1], &Y);
    nc_inq_dimlen(ncid, dimids[2], &X);

    if (static_cast<size_t>(y_end) > Y) y_end = Y;
    size_t region_height = y_end - y_start;

    // Выделяем память для данных: data[T][region_height][X]
    data.resize(T, std::vector<std::vector<double>>(region_height, std::vector<double>(X, 0.0)));

    // Задаем массивы start и count для чтения
    size_t start[3] = { 0, static_cast<size_t>(y_start), 0 };
    size_t count[3] = { T, region_height, X };

    // Буфер для хранения всех данных
    std::vector<double> buffer(T * region_height * X, 0.0);
    retval = nc_get_vara_double(ncid, varid, start, count, buffer.data());
    if (retval != NC_NOERR) {
        std::cerr << "Ошибка чтения переменной 'height': " << nc_strerror(retval) << std::endl;
        nc_close(ncid);
        return data;
    }
    nc_close(ncid);

    // Заполняем data из буфера
    for (size_t t = 0; t < T; t++) {
        for (size_t i = 0; i < region_height; i++) {
            for (size_t x = 0; x < X; x++) {
                size_t idx = t * region_height * X + i * X + x;
                data[t][i][x] = buffer[idx];
            }
        }
    }
    return data;
}



std::vector<std::vector<std::vector<std::vector<double>>>> BasisManager::get_fk_region(int y_start, int y_end) {
    std::vector<std::vector<std::vector<std::vector<double>>>> fk;
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.path().extension() == ".nc")
            files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());

    std::vector<std::future<std::vector<std::vector<std::vector<double>>>>> futures;

    for (const auto& file : files) {
        futures.push_back(std::async(std::launch::deferred, [file, y_start, y_end]() -> std::vector<std::vector<std::vector<double>>> {
            std::vector<std::vector<std::vector<double>>> data;
            std::cout << file << std::endl;
            int ncid;

            // Синхронизация открытия файла и других вызовов netCDF
            {
          
                if (open_nc_file(file.string(), ncid) != NC_NOERR)
                    return data;
            }

            int varid;
            {
             
                int retval = nc_inq_varid(ncid, "height", &varid);
                if (retval != NC_NOERR) {
                    std::cerr << "Переменная 'height' не найдена в " << file.string() << std::endl;
                    nc_close(ncid);
                    return data;
                }
            }

            int ndims;
            {
        
                nc_inq_varndims(ncid, varid, &ndims);
            }
            if (ndims != 3) {
                std::cerr << "Ожидалось 3 измерения в файле " << file.string() << std::endl;
       
                nc_close(ncid);
                return data;
            }

            int dimids[3];
            size_t T, Y, X;
            {
        
                nc_inq_vardimid(ncid, varid, dimids);
                nc_inq_dimlen(ncid, dimids[0], &T);
                nc_inq_dimlen(ncid, dimids[1], &Y);
                nc_inq_dimlen(ncid, dimids[2], &X);
            }

            int local_y_end = y_end;
            if (static_cast<size_t>(local_y_end) > Y) local_y_end = Y;
            size_t region_height = local_y_end - y_start;
            data.resize(T, std::vector<std::vector<double>>(region_height, std::vector<double>(X, 0.0)));

            size_t start[3] = { 0, static_cast<size_t>(y_start), 0 };
            size_t count[3] = { T, region_height, X };
            std::vector<double> buffer(T * region_height * X, 0.0);

            {
                //std::lock_guard<std::mutex> lock(netcdf_mutex);
                int retval = nc_get_vara_double(ncid, varid, start, count, buffer.data());
                if (retval != NC_NOERR) {
                    std::cerr << "Ошибка чтения файла " << file.string() << " : " << nc_strerror(retval) << std::endl;
                    nc_close(ncid);
                    return data;
                }
                nc_close(ncid);
            }

            for (size_t t = 0; t < T; t++) {
                for (size_t i = 0; i < region_height; i++) {
                    for (size_t x = 0; x < X; x++) {
                        size_t idx = t * region_height * X + i * X + x;
                        data[t][i][x] = buffer[idx];
                    }
                }
            }
            return data;
            }));
    }

    for (auto& fut : futures) {
        auto file_data = fut.get();
        if (!file_data.empty()) {
            fk.push_back(file_data);
        }
    }
    return fk;
}