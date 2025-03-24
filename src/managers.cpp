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
#define GET_TIME(code) do { code; } while(0)

int open_nc_file(const std::string& filename, int& ncid) {
    int retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Error opening file " << filename << " : " << nc_strerror(retval) << std::endl;
    }
    return retval;
}

std::vector<double> read_nc_file(const std::string& filename, int y_start, int y_end) {
    std::vector<double> result;
    fs::path file(filename);
    std::cout << "Loading file (mapped): " << file << std::endl;

    // Получаем размеры из файла.
    int ncid;
    if (open_nc_file(file.string(), ncid) != NC_NOERR)
        return result;

    int varid;
    int retval = nc_inq_varid(ncid, "height", &varid);
    if (retval != NC_NOERR) {
        std::cerr << "Variable 'height' not found in file " << file.string() << std::endl;
        nc_close(ncid);
        return result;
    }

    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);
    if (ndims != 3) {
        std::cerr << "Expected 3 dimensions in file " << file.string() << std::endl;
        nc_close(ncid);
        return result;
    }

    int dimids[3];
    size_t T, Y, X;
    nc_inq_vardimid(ncid, varid, dimids);
    nc_inq_dimlen(ncid, dimids[0], &T);
    nc_inq_dimlen(ncid, dimids[1], &Y);
    nc_inq_dimlen(ncid, dimids[2], &X);
    nc_close(ncid);

    int local_y_end = y_end;
    if (static_cast<size_t>(local_y_end) > Y)
        local_y_end = static_cast<int>(Y);
    int region_height = local_y_end - y_start;

    size_t totalElements = T * region_height * X;
    size_t totalBytes = totalElements * sizeof(double);

    // Создаем уникальное имя для общей памяти.
    std::ostringstream shmNameStream;
    shmNameStream << "Local\\MySharedMemory_" << GetCurrentProcessId() << "_" << GetTickCount();
    std::string shmName = shmNameStream.str();

    HANDLE hMapFile = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
        0, static_cast<DWORD>(totalBytes), shmName.c_str());
    if (hMapFile == NULL) {
        std::cerr << "Could not create file mapping object (" << GetLastError() << ")\n";
        return result;
    }

    // Формируем командную строку для запуска дочернего процесса.
    std::string childExe = "nc_reader_child.exe"; // путь должен быть корректным
    std::ostringstream oss;
    oss << "\"" << childExe << "\" "
        << "\"" << file.string() << "\" "
        << y_start << " "
        << region_height << " "
        << T << " "
        << X << " "
        << "\"" << shmName << "\"";
    std::string commandLine = oss.str();

    // Запускаем дочерний процесс с повторными попытками в случае неудачи.
    const int maxAttempts = 3;
    int attempt = 0;
    DWORD exitCode = 1; // по умолчанию ненулевой
    while (attempt < maxAttempts && exitCode != 0) {
        PROCESS_INFORMATION piProcInfo;
        ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));
        STARTUPINFO siStartInfo;
        ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
        siStartInfo.cb = sizeof(STARTUPINFO);

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
            return result;
        }

        WaitForSingleObject(piProcInfo.hProcess, INFINITE);
        if (!GetExitCodeProcess(piProcInfo.hProcess, &exitCode)) {
            std::cerr << "Failed to get child process exit code (" << GetLastError() << ")\n";
            exitCode = 1;
        }
        CloseHandle(piProcInfo.hProcess);
        CloseHandle(piProcInfo.hThread);

        if (exitCode != 0) {
            std::cerr << "Child process returned " << exitCode << ", retrying ("
                << (attempt + 1) << "/" << maxAttempts << ")...\n";
        }
        attempt++;
    }
    if (exitCode != 0) {
        std::cerr << "Child process failed after " << maxAttempts << " attempts.\n";
        CloseHandle(hMapFile);
        return result;
    }

    // Мапим общую память для чтения.
    LPVOID pBuf = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, totalBytes);
    if (pBuf == NULL) {
        std::cerr << "Could not map view of file (" << GetLastError() << ")\n";
        CloseHandle(hMapFile);
        return result;
    }

    // Копируем данные из общей памяти в вектор.
    result.resize(totalElements);
    memcpy(result.data(), pBuf, totalBytes);

    UnmapViewOfFile(pBuf);
    CloseHandle(hMapFile);

    return result;
}

std::vector<std::vector<double>> BasisManager::get_fk_region(int y_start, int y_end) {
    std::vector<std::vector<double>> fk;
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

    auto extractIndex = [](const fs::path& filePath) -> int {
        std::regex regexPattern("_(\\d+)\\.nc");
        std::smatch match;
        std::string filename = filePath.filename().string();
        if (std::regex_search(filename, match, regexPattern)) {
            return std::stoi(match[1].str());
        }
        return std::numeric_limits<int>::max();
        };

    std::sort(files.begin(), files.end(), [&](const fs::path& a, const fs::path& b) {
        return extractIndex(a) < extractIndex(b);
        });

    std::vector<std::future<std::vector<double>>> futures;
    // Запускаем асинхронное выполнение для каждого файла
    for (const auto& file : files) {
        futures.emplace_back(std::async(std::launch::async, [file, y_start, y_end]() -> std::vector<double> {
            return read_nc_file(file.string(), y_start, y_end);
            }));
    }

    // Сбор результатов
    for (auto& fut : futures) {
        std::vector<double> file_data = fut.get();
        if (!file_data.empty()) {
            fk.push_back(std::move(file_data));
        }
    }
    return fk;
}

std::vector<double> WaveManager::load_mariogramm_by_region(int y_start, int y_end) {
    return read_nc_file(nc_file, y_start, y_end);
}
