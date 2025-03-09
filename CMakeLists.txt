cmake_minimum_required(VERSION 3.14)
project(TsunamiCoefficientsCalculator VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(DH5_HAVE_THREADSAFE)
# Найти HDF5 (vcpkg предоставляет импортированную цель)
find_package(HDF5 REQUIRED)
message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")

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

# Подключение Eigen через FetchContent
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)
message(STATUS "Eigen успешно получен: ${eigen_SOURCE_DIR}")

# Подключение nlohmann/json через FetchContent
include(FetchContent)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.2
)
FetchContent_MakeAvailable(json)
message(STATUS "nlohmann/json успешно получен: ${json_SOURCE_DIR}")

# Добавляем исходники проекта
add_executable(TsunamiCoefficientsCalculator
    src/main.cpp
    src/approx_orto.cpp
    src/approx_orto.h
    src/stable_data_structs.h
    src/managers.h
    src/managers.cpp
    src/statistics.h
    src/statistics.cpp
)
message(STATUS "Добавлена исполняемая цель TsunamiCoefficientsCalculator со списком исходников")

# Линкуем с HDF5 C++ API
# ВАЖНО: поскольку цель hdf5::hdf5_cpp-shared не найдена, используем hdf5::hdf5-shared
target_link_libraries(TsunamiCoefficientsCalculator PRIVATE hdf5::hdf5-shared)
message(STATUS "Ссылка на hdf5::hdf5-shared добавлена для TsunamiCoefficientsCalculator")

# Добавляем пути к заголовочным файлам Eigen
target_include_directories(TsunamiCoefficientsCalculator PRIVATE ${eigen_SOURCE_DIR})
message(STATUS "Пути к заголовочным файлам Eigen добавлены: ${eigen_SOURCE_DIR}")
