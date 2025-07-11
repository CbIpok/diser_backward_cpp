// nc_reader_child.cpp
#define NOMINMAX 
#include <netcdf.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <windows.h>
#include <algorithm>
#include <stdexcept>
#include <string>

// Структура для хранения размерностей переменной из NetCDF-файла
struct NcDimensions {
    size_t T;
    size_t Y;
    size_t X;
};

// Функция чтения размерностей переменной "height" из NetCDF-файла
NcDimensions read_nc_dimensions(int ncid, int varid) {
    NcDimensions dims;
    int dimids[3];
    int retval = nc_inq_vardimid(ncid, varid, dimids);
    if (retval != NC_NOERR) {
        throw std::runtime_error("Ошибка получения идентификаторов измерений: " + std::string(nc_strerror(retval)));
    }
    retval = nc_inq_dimlen(ncid, dimids[0], &dims.T);
    if (retval != NC_NOERR) {
        throw std::runtime_error("Ошибка чтения размера измерения T: " + std::string(nc_strerror(retval)));
    }
    retval = nc_inq_dimlen(ncid, dimids[1], &dims.Y);
    if (retval != NC_NOERR) {
        throw std::runtime_error("Ошибка чтения размера измерения Y: " + std::string(nc_strerror(retval)));
    }
    retval = nc_inq_dimlen(ncid, dimids[2], &dims.X);
    if (retval != NC_NOERR) {
        throw std::runtime_error("Ошибка чтения размера измерения X: " + std::string(nc_strerror(retval)));
    }
    return dims;
}

int main(int argc, char* argv[]) {
    // Новое использование: ожидается <filename> <y_start> <region_height> [<shmName>]
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <filename> <y_start> <region_height> [<shmName>]" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    int y_start = std::atoi(argv[2]);
    int region_height_arg = std::atoi(argv[3]);

    bool useSharedMemory = (argc >= 5);
    std::string shmName;
    if (useSharedMemory) {
        shmName = argv[4];
    }

    int ncid;
    int retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Ошибка открытия файла " << filename << ": " << nc_strerror(retval) << std::endl;
        return 1;
    }

    int varid;
    retval = nc_inq_varid(ncid, "height", &varid);
    if (retval != NC_NOERR) {
        std::cerr << "Переменная 'height' не найдена в файле " << filename << std::endl;
        nc_close(ncid);
        return 1;
    }

    // Получаем размерности из файла
    NcDimensions dims;
    try {
        dims = read_nc_dimensions(ncid, varid);
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        nc_close(ncid);
        return 1;
    }

    // Корректировка region_height: не превышать оставшееся количество строк в файле
    int region_height = std::min(region_height_arg, static_cast<int>(dims.Y) - y_start);
    size_t T = dims.T;
    size_t X = dims.X;
    size_t totalElements = T * region_height * X;

    // Создаём общую память или выполняем вывод в stdout
    HANDLE hMapFile = NULL;
    LPVOID pBuf = NULL;
    if (useSharedMemory) {
        // Формируем уникальное имя для общей памяти (можно оставить прежним способом)
        DWORD pid = GetCurrentProcessId();
        DWORD tick = GetTickCount();
        shmName += "_" + std::to_string(pid) + "_" + std::to_string(tick);
        hMapFile = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
            0, static_cast<DWORD>(totalElements * sizeof(double)), shmName.c_str());
        if (hMapFile == NULL) {
            std::cerr << "Не удалось создать объект отображения файла (" << GetLastError() << ")\n";
            nc_close(ncid);
            return 1;
        }
        pBuf = MapViewOfFile(hMapFile, FILE_MAP_WRITE, 0, 0, totalElements * sizeof(double));
        if (pBuf == NULL) {
            std::cerr << "Не удалось отобразить файл (" << GetLastError() << ")\n";
            CloseHandle(hMapFile);
            nc_close(ncid);
            return 1;
        }
    }

    // Считываем данные из NetCDF в буфер
    size_t start[3] = { 0, static_cast<size_t>(y_start), 0 };
    size_t count[3] = { T, static_cast<size_t>(region_height), X };
    std::vector<double> buffer(totalElements);
    retval = nc_get_vara_double(ncid, varid, start, count, buffer.data());
    if (retval != NC_NOERR) {
        std::cerr << "Ошибка чтения данных: " << nc_strerror(retval) << std::endl;
        nc_close(ncid);
        if (useSharedMemory) {
            UnmapViewOfFile(pBuf);
            CloseHandle(hMapFile);
        }
        return 1;
    }
    nc_close(ncid);

    // Переупорядочивание данных:
    // Исходный порядок (из файла): [T][region_height][X]
    // Желаемый порядок: [region_height][X][T]
    if (useSharedMemory) {
        double* dest = reinterpret_cast<double*>(pBuf);
        for (size_t t = 0; t < T; t++) {
            for (int y = 0; y < region_height; y++) {
                for (size_t x = 0; x < X; x++) {
                    size_t src_idx = t * (region_height * X) + y * X + x;
                    size_t dst_idx = y * (X * T) + x * T + t;
                    dest[dst_idx] = buffer[src_idx];
                }
            }
        }
        UnmapViewOfFile(pBuf);
        CloseHandle(hMapFile);
    }
    else {
        // Фолбэк: переупорядочивание в дополнительный буфер и вывод в stdout
        std::vector<double> reordered(totalElements);
        for (size_t t = 0; t < T; t++) {
            for (int y = 0; y < region_height; y++) {
                for (size_t x = 0; x < X; x++) {
                    size_t src_idx = t * (region_height * X) + y * X + x;
                    size_t dst_idx = y * (X * T) + x * T + t;
                    reordered[dst_idx] = buffer[src_idx];
                }
            }
        }
        std::cout.write(reinterpret_cast<const char*>(reordered.data()), totalElements * sizeof(double));
        std::cout.flush();
    }

    return 0;
}
