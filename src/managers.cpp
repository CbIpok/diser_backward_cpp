#define NOMINMAX 
#include "managers.h"
#include <windows.h>
#include <string>
#include <sstream>
#include <netcdf.h>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <vector>
#include <future>
#include <regex>
namespace fs = std::filesystem;
#define GET_TIME(code) \
    do { \
       /* auto start = std::chrono::high_resolution_clock::now(); */\
        code; \
        /*auto end = std::chrono::high_resolution_clock::now(); \
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); \
        { \
            static std::mutex cout_mutex; \
            std::lock_guard<std::mutex> lock(cout_mutex); \
            std::cout << __FILE__ << " line " << __LINE__ \
                      << " time: " << duration << "ms" \
                      << " note: " << file << std::endl; \
        }*/ \
    } while(0)
// Ôóíêöèÿ äëÿ îòêðûòèÿ NetCDF-ôàéëà ñ ïðîâåðêîé îøèáîê
int open_nc_file(const std::string& filename, int& ncid) {
    int retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Îøèáêà îòêðûòèÿ ôàéëà " << filename << " : " << nc_strerror(retval) << std::endl;
    }
    return retval;
}
std::vector<std::vector<std::vector<double>>> read_nc_file(const fs::path& file, int y_start, int y_end) {
    std::vector<std::vector<std::vector<double>>> data;
    std::cout << "Loading file: " << file << std::endl;
    int ncid;

    GET_TIME(
        // Открываем файл для получения размеров.
        if (open_nc_file(file.string(), ncid) != NC_NOERR)
            return data;
    );

    int varid;
    int retval = nc_inq_varid(ncid, "height", &varid);
    if (retval != NC_NOERR) {
        std::cerr << "Variable 'height' not found in file " << file.string() << std::endl;
        nc_close(ncid);
        return data;
    }

    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);
    if (ndims != 3) {
        std::cerr << "Expected 3 dimensions in file " << file.string() << std::endl;
        nc_close(ncid);
        return data;
    }

    int dimids[3];
    size_t T, Y, X;
    nc_inq_vardimid(ncid, varid, dimids);
    nc_inq_dimlen(ncid, dimids[0], &T);
    nc_inq_dimlen(ncid, dimids[1], &Y);
    nc_inq_dimlen(ncid, dimids[2], &X);

    int local_y_end = y_end;
    if (static_cast<size_t>(local_y_end) > Y)
        local_y_end = static_cast<int>(Y);
    int region_height = local_y_end - y_start;

    GET_TIME(
        // Предвыделяем 3D-структуру для данных.
        data.resize(T);
    for (size_t t = 0; t < T; t++) {
        data[t].resize(region_height);
        for (int i = 0; i < region_height; i++) {
            data[t][i].resize(X);
        }
    });

    nc_close(ncid); // Получение размеров завершено.

    // Вычисляем общий размер в байтах.
    size_t totalBytes = T * region_height * X * sizeof(double);

    // Создаём уникальное имя общей памяти.
    std::ostringstream shmNameStream;
    shmNameStream << "Local\\MySharedMemory_" << GetCurrentProcessId() << "_" << GetTickCount();
    std::string shmName = shmNameStream.str();

    HANDLE hMapFile = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
        0, static_cast<DWORD>(totalBytes), shmName.c_str());
    if (hMapFile == NULL) {
        std::cerr << "Could not create file mapping object (" << GetLastError() << ")\n";
        return data;
    }

    // Формируем командную строку для запуска дочернего процесса с передачей имени общей памяти.
    std::string childExe = "nc_reader_child.exe"; // Убедитесь, что путь указан верно.
    std::ostringstream oss;
    oss << "\"" << childExe << "\" "
        << "\"" << file.string() << "\" "
        << y_start << " "
        << region_height << " "
        << T << " "
        << X << " "
        << "\"" << shmName << "\"";
    std::string commandLine = oss.str();

    // Запускаем дочерний процесс.
    PROCESS_INFORMATION piProcInfo;
    GET_TIME(
        ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));
    );
    STARTUPINFO siStartInfo;
    GET_TIME(
        ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
    );
    siStartInfo.cb = sizeof(STARTUPINFO);

    GET_TIME(
        BOOL bSuccess = CreateProcess(NULL,
            const_cast<LPSTR>(commandLine.c_str()),
            NULL,
            NULL,
            TRUE,
            0,
            NULL,
            NULL,
            &siStartInfo,
            &piProcInfo);
    if (!bSuccess) {
        std::cerr << "CreateProcess failed (" << GetLastError() << ")\n";
        CloseHandle(hMapFile);
        return data;
    });

    GET_TIME(
        // Ожидаем завершения дочернего процесса.
        WaitForSingleObject(piProcInfo.hProcess, INFINITE);
    CloseHandle(piProcInfo.hProcess);
    CloseHandle(piProcInfo.hThread););
    LPVOID pBuf;
    GET_TIME(
        // Отображаем общую память для чтения.
        pBuf = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, totalBytes);
    if (pBuf == NULL) {
        std::cerr << "Could not map view of file (" << GetLastError() << ")\n";
        CloseHandle(hMapFile);
        return data;
    });

    // Копируем данные из общей памяти в 3D-вектор.
    double* dptr = reinterpret_cast<double*>(pBuf);
    GET_TIME(
        for (size_t t = 0; t < T; t++) {
            for (int i = 0; i < region_height; i++) {
                size_t idx = t * region_height * X + i * X;
                std::copy(dptr + idx, dptr + idx + X, data[t][i].begin());
            }
        });
    GET_TIME(
        UnmapViewOfFile(pBuf);
    CloseHandle(hMapFile);
        );
    return data;
}


// Ðåàëèçàöèÿ ìåòîäà WaveManager::load_mariogramm_by_region ñ èñïîëüçîâàíèåì netcdf.h
std::vector<std::vector<std::vector<double>>> WaveManager::load_mariogramm_by_region(int y_start, int y_end) {

    return read_nc_file(nc_file, y_start, y_end);
}


// Ôóíêöèÿ äëÿ èçâëå÷åíèÿ èíäåêñà èç èìåíè ôàéëà
int extractIndex(const fs::path& filePath) {
    // Ðåãóëÿðíîå âûðàæåíèå äëÿ ïîèñêà øàáëîíà _<÷èñëî>.nc
    std::regex regexPattern("_(\\d+)\\.nc");
    std::smatch match;
    std::string filename = filePath.filename().string();
    if (std::regex_search(filename, match, regexPattern)) {
        return std::stoi(match[1].str());
    }
    // Åñëè èíäåêñ íå íàéäåí, âîçâðàùàåì ìàêñèìàëüíîå çíà÷åíèå,
    // ÷òîáû ôàéë îêàçàëñÿ â êîíöå îòñîðòèðîâàííîãî ñïèñêà.
    return std::numeric_limits<int>::max();
}

// Ôóíêöèÿ äëÿ ïîëó÷åíèÿ îòñîðòèðîâàííîãî ñïèñêà ôàéëîâ
std::vector<fs::path> getSortedFileList(const std::string& folder) {
    std::vector<fs::path> files;

    // Ïåðåáèðàåì âñå ôàéëû â çàäàííîé äèðåêòîðèè
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.is_regular_file()) {
            fs::path filePath = entry.path();
            // Ïðîâåðÿåì, ÷òî ôàéë èìååò ðàñøèðåíèå ".nc" è ñîäåðæèò ñèìâîë '_'
            if (filePath.extension() == ".nc" &&
                filePath.filename().string().find('_') != std::string::npos) {
                files.push_back(filePath);
            }
        }
    }

    // Ñîðòèðóåì ôàéëû ïî ÷èñëîâîìó çíà÷åíèþ, èçâëå÷¸ííîìó èç èìåíè ôàéëà
    std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
        return extractIndex(a) < extractIndex(b);
        });

    return files;
}


std::vector<std::vector<std::vector<std::vector<double>>>> BasisManager::get_fk_region(int y_start, int y_end) {
    std::vector<std::vector<std::vector<std::vector<double>>>> fk;
    std::vector<fs::path> files = getSortedFileList(folder);

    // Создаем вектор будущих результатов (фьючерсов)
    std::vector<std::future<std::vector<std::vector<std::vector<double>>>>> futures;

    // Запускаем read_nc_file параллельно для каждого файла
    for (const auto& file : files) {
        futures.push_back(std::async(std::launch::async, [this, file, y_start, y_end]() {
            return read_nc_file(file, y_start, y_end);
            }));
    }

    // Синхронно собираем результаты
    for (auto& fut : futures) {
        auto file_data = fut.get();
        if (!file_data.empty()) {
            fk.push_back(file_data);
        }
    }

    return fk;
}