cmake_minimum_required(VERSION 3.14)
project(TsunamiCoefficientsCalculator VERSION 1.0 LANGUAGES CXX)

# Используем стандарт C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Подключаем Eigen через FetchContent
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

# Подключаем nlohmann/json (можно добавить через FetchContent или vcpkg)
include(FetchContent)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.2
)
FetchContent_MakeAvailable(json)

# Указываем путь до установленного NetCDF (NetCDF‑C, а также NetCDF‑C++ API если требуется)
set(NetCDF_DIR "C:/Program Files/netCDF_4.9.3")
include_directories("${NetCDF_DIR}/include")
link_directories("${NetCDF_DIR}/lib")

# Добавляем исходники проекта
add_executable(TsunamiCoefficientsCalculator
    src/main.cpp
    src/approx_orto.cpp
    src/approx_orto.h
    src/stable_data_structs.h
    src/managers.h
    src/json.hpp
    src/managers.cpp
    src/statistics.h
    src/statistics.cpp
)

# Линкуем с библиотекой netcdf.lib (убедитесь, что она находится в NetCDF_DIR/lib)
target_link_libraries(TsunamiCoefficientsCalculator PRIVATE netcdf.lib)

# Если требуется, можно добавить явное указание для Eigen (хотя FetchContent уже подключает заголовки):
target_include_directories(TsunamiCoefficientsCalculator PRIVATE ${eigen_SOURCE_DIR})