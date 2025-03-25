// nc_reader_child.cpp
#include <netcdf.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <windows.h>

int main(int argc, char* argv[]) {
    // Ожидаем: <filename> <y_start> <region_height> <T> <X> [<shmName>]
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <filename> <y_start> <region_height> <T> <X> [<shmName>]" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    int y_start = std::atoi(argv[2]);
    int region_height = std::atoi(argv[3]);
    size_t T = static_cast<size_t>(std::atoi(argv[4]));
    size_t X = static_cast<size_t>(std::atoi(argv[5]));
    size_t total_elements = T * region_height * X;

    bool useSharedMemory = (argc >= 7);
    HANDLE hMapFile = NULL;
    LPVOID pBuf = NULL;
    if (useSharedMemory) {
        std::string shmName = argv[6];
        hMapFile = OpenFileMapping(FILE_MAP_WRITE, FALSE, shmName.c_str());
        if (hMapFile == NULL) {
            std::cerr << "Could not open shared memory (" << GetLastError() << ")\n";
            return 1;
        }
        pBuf = MapViewOfFile(hMapFile, FILE_MAP_WRITE, 0, 0, total_elements * sizeof(double));
        if (pBuf == NULL) {
            std::cerr << "Could not map view of file (" << GetLastError() << ")\n";
            CloseHandle(hMapFile);
            return 1;
        }
    }

    int ncid;
    int retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Error opening file " << filename << ": " << nc_strerror(retval) << std::endl;
        if (useSharedMemory) {
            UnmapViewOfFile(pBuf);
            CloseHandle(hMapFile);
        }
        return 1;
    }

    int varid;
    retval = nc_inq_varid(ncid, "height", &varid);
    if (retval != NC_NOERR) {
        std::cerr << "Variable 'height' not found in file " << filename << std::endl;
        nc_close(ncid);
        if (useSharedMemory) {
            UnmapViewOfFile(pBuf);
            CloseHandle(hMapFile);
        }
        return 1;
    }

    // Считываем данные из NetCDF в буфер в порядке [T][region_height][X]
    size_t start[3] = { 0, static_cast<size_t>(y_start), 0 };
    size_t count[3] = { T, static_cast<size_t>(region_height), X };
    std::vector<double> buffer(total_elements);
    retval = nc_get_vara_double(ncid, varid, start, count, buffer.data());
    if (retval != NC_NOERR) {
        std::cerr << "Error reading data: " << nc_strerror(retval) << std::endl;
        nc_close(ncid);
        if (useSharedMemory) {
            UnmapViewOfFile(pBuf);
            CloseHandle(hMapFile);
        }
        return 1;
    }
    nc_close(ncid);

    if (useSharedMemory) {
        // Переупорядочиваем данные напрямую в общей памяти.
        // Исходный порядок (из файла): [T][region_height][X]
        // Индекс: src_idx = t * (region_height * X) + y * X + x
        // Желаемый порядок: [region_height][X][T]
        // Индекс: dst_idx = y * (X * T) + x * T + t
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

        // Завершаем работу с общей памятью
        UnmapViewOfFile(pBuf);
        CloseHandle(hMapFile);
    }
    else {
        // Фолбэк: если общая память не используется, выполняем переупорядочивание в дополнительный буфер и выводим в stdout.
        std::vector<double> reordered(total_elements);
        for (size_t t = 0; t < T; t++) {
            for (int y = 0; y < region_height; y++) {
                for (size_t x = 0; x < X; x++) {
                    size_t src_idx = t * (region_height * X) + y * X + x;
                    size_t dst_idx = y * (X * T) + x * T + t;
                    reordered[dst_idx] = buffer[src_idx];
                }
            }
        }
        std::cout.write(reinterpret_cast<const char*>(reordered.data()), total_elements * sizeof(double));
        std::cout.flush();
    }

    return 0;
}
