cmake_minimum_required(VERSION 3.14)
project(TsunamiCoefficientsCalculator VERSION 1.0 LANGUAGES CXX)

# Используем стандарт C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
<<<<<<< Updated upstream

=======
set(DH5_HAVE_THREADSAFE)
# Найти HDF5 (vcpkg предоставляет импортированную цель)
find_package(HDF5 REQUIRED)
message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
find_package(adios2 CONFIG REQUIRED)
# Проверка наличия цели hdf5::hdf5_cpp-shared
if(TARGET hdf5::hdf5_cpp-shared)
    message(STATUS "Импортированная цель hdf5::hdf5_cpp-shared найдена")
else()
    message(WARNING "Импортированная цель hdf5::hdf5_cpp-shared не найдена. Доступные цели HDF5:")
    get_cmake_property(_targets TARGETS)
    foreach(target ${_targets})
        if(target MATCHES "hdf5")
            message(STATUS "  ${target}")
        endif()
    endforeach()
endif()
>>>>>>> Stashed changes

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
<<<<<<< Updated upstream
    src/statistics.cpp
)

# Линкуем с библиотекой netcdf.lib (убедитесь, что она находится в NetCDF_DIR/lib)
target_link_libraries(TsunamiCoefficientsCalculator PRIVATE netcdf.lib)

# Если требуется, можно добавить явное указание для Eigen (хотя FetchContent уже подключает заголовки):
target_include_directories(TsunamiCoefficientsCalculator PRIVATE ${eigen_SOURCE_DIR})
=======
    src/statistics.cpp)
message(STATUS "Добавлена исполняемая цель TsunamiCoefficientsCalculator со списком исходников")

# Линкуем с HDF5 C++ API
# ВАЖНО: поскольку цель hdf5::hdf5_cpp-shared не найдена, используем hdf5::hdf5-shared
target_link_libraries(TsunamiCoefficientsCalculator PRIVATE hdf5::hdf5-shared adios2::adios2)
message(STATUS "Ссылка на hdf5::hdf5-shared добавлена для TsunamiCoefficientsCalculator")

# Добавляем пути к заголовочным файлам Eigen
target_include_directories(TsunamiCoefficientsCalculator PRIVATE ${eigen_SOURCE_DIR})
message(STATUS "Пути к заголовочным файлам Eigen добавлены: ${eigen_SOURCE_DIR}")
add_executable(convert_h5_to_bp src/convector.cpp)
target_link_libraries(convert_h5_to_bp PRIVATE adios2::adios2  hdf5::hdf5_cpp-shared)
>>>>>>> Stashed changes
