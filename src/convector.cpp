#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include "H5Cpp.h"
#include <adios2.h>

namespace fs = std::filesystem;

// Функция, выполняющая преобразование одного файла
bool convertFile(const fs::path& inputFile, const fs::path& outputFile)
{
    try
    {
        // Открываем HDF5-файл для чтения
        H5::H5File h5file(inputFile.string(), H5F_ACC_RDONLY);
        H5::DataSet dataset = h5file.openDataSet("height");

        // Получаем dataspace и проверяем число измерений
        H5::DataSpace dataspace = dataset.getSpace();
        int ndims = dataspace.getSimpleExtentNdims();
        if (ndims != 3)
        {
            std::cerr << "Ожидалось 3 измерения, а получено " << ndims
                << " в файле: " << inputFile << std::endl;
            return false;
        }

        // Получаем размеры: [T, Y, X]
        hsize_t dims[3];
        dataspace.getSimpleExtentDims(dims, nullptr);
        size_t T = dims[0], Y = dims[1], X = dims[2];

        // Читаем все данные в буфер
        std::vector<double> data(T * Y * X);
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
        std::cout << "HDF5-файл успешно прочитан: " << inputFile << std::endl;

        // Инициализируем ADIOS2 для записи BP-файла
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("WriteIO");
        io.SetEngine("BPFile");

        // Определяем оператор для сжатия ZFP
        auto zfpOp = adios.DefineOperator("zfpOp", "zfp");

        // Открываем BP-файл для записи
        adios2::Engine writer = io.Open(outputFile.string(), adios2::Mode::Write);

        // Определяем переменную "height" с глобальными размерами и полной выборкой
        adios2::Variable<double> var = io.DefineVariable<double>(
            "height", { T, Y, X }, { 0, 0, 0 }, { T, Y, X });

        // Добавляем операцию сжатия: параметр "Accuracy" задаёт требуемую точность
        var.AddOperation(zfpOp, { {"Accuracy", "0.001"} });

        // Записываем данные в BP-файл
        writer.BeginStep();
        writer.Put(var, data.data());
        writer.EndStep();
        writer.Close();

        std::cout << "Конвертация завершена: " << outputFile << std::endl;
    }
    catch (H5::Exception& err)
    {
        std::cerr << "Ошибка HDF5: " << err.getDetailMsg()
            << " в файле: " << inputFile << std::endl;
        return false;
    }
    catch (std::exception& ex)
    {
        std::cerr << "Ошибка: " << ex.what()
            << " в файле: " << inputFile << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char* argv[])
{
 /*   if (argc < 3)
    {
        std::cerr << "Использование: " << argv[0]
            << " <input_directory> <output_directory>" << std::endl;
        return 1;
    }*/

    //fs::path inputDir = "T:/tsunami_res_folder";
    //fs::path outputDir = "T:/tsunami_adios_res";

    //if (!fs::exists(inputDir) || !fs::is_directory(inputDir))
    //{
    //    std::cerr << "Входная директория не существует или не является директорией: "
    //        << inputDir << std::endl;
    //    return 1;
    //}

    //// Рекурсивный обход входной директории
    //for (const auto& entry : fs::recursive_directory_iterator(inputDir))
    //{
    //    if (entry.is_regular_file())
    //    {
    //        fs::path inputFilePath = entry.path();
    //        // Обрабатываем только файлы с расширением ".nc"
    //        if (inputFilePath.extension() == ".nc")
    //        {
    //            // Получаем путь относительно входной директории
    //            fs::path relativePath = fs::relative(inputFilePath, inputDir);
    //            // Заменяем расширение на .bp
    //            relativePath.replace_extension(".bp");
    //            fs::path outputFilePath = outputDir / relativePath;

    //            // Создаём каталог для выходного файла, если он отсутствует
    //            fs::create_directories(outputFilePath.parent_path());

    //            std::cout << "Конвертация: " << inputFilePath
    //                << " -> " << outputFilePath << std::endl;
    //            if (!convertFile(inputFilePath, outputFilePath))
    //            {
    //                std::cerr << "Не удалось конвертировать файл: "
    //                    << inputFilePath << std::endl;
    //            }
    //        }
    //    }
    //}

    std::string bpFile = "T:/tsunami_adios_res/parabola_200_2000/basis_6/basis_0.bp";

    try
    {
        // Инициализация ADIOS2 с использованием движка BPFile для чтения
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("BPFile");

        // Открываем BP-файл для чтения
        adios2::Engine reader = io.Open(bpFile, adios2::Mode::Read);

        // Начинаем шаг чтения для получения метаданных
        reader.BeginStep();

        std::cout << "BP File: " << bpFile << std::endl;
        std::cout << "Available Variables:" << std::endl;

        // Получаем список переменных, доступных в файле
        const auto& variables = io.AvailableVariables();
        for (const auto& varPair : variables)
        {
            const std::string& varName = varPair.first;
            const auto& params = varPair.second; // std::map<std::string, std::string>

            std::cout << "----------------------------------------" << std::endl;
            std::cout << "Variable: " << varName << std::endl;

            // Вывод параметров переменной (например, тип, размеры, и т.д.)
            for (const auto& p : params)
            {
                std::cout << "  " << p.first << " : " << p.second << std::endl;
            }
        }

        // Вывод атрибутов IO (если имеются)
        const auto& ioAttributes = io.AvailableAttributes();
        if (!ioAttributes.empty())
        {
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "IO Attributes:" << std::endl;
            for (const auto& attrPair : ioAttributes)
            {
                const std::string& attrName = attrPair.first;
                const auto& params = attrPair.second; // std::map<std::string, std::string>

                std::cout << "Attribute: " << attrName << std::endl;
                for (const auto& p : params)
                {
                    std::cout << "  " << p.first << " : " << p.second << std::endl;
                }
            }
        }

        reader.EndStep();
        reader.Close();
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
