#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "approx_orto.h"
#include "stable_data_structs.h"
#include "statistics.h"
// Функция для проведения тестов
void run_tests() {
    // Размерности для тестовых случаев
    std::vector<int> dimensions = { 3, 4, 5, 6, 8 };
    // Порог допуска (например, 1e-6)
    double tol = 1e-6;
    bool all_passed = true;

    std::cout << "tests (approximate_with_non_orthogonal_basis_orto):\n";

    for (int n : dimensions) {
        // Генерируем случайную квадратную матрицу n x n,
        // представляющую базис (каждая строка – базисный вектор).
        Eigen::MatrixXd M;
        // Обеспечиваем обратимость (регенерируем, если определитель слишком мал)
        do {
            M = Eigen::MatrixXd::Random(n, n);
        } while (std::abs(M.determinant()) < 1e-3);

        // Генерируем случайный вектор коэффициентов c длины n.
        Eigen::VectorXd c = Eigen::VectorXd::Random(n);
        // Вычисляем вектор x как линейную комбинацию базисных векторов:
        // x = c[0]*M.row(0) + ... + c[n-1]*M.row(n-1).
        // Чтобы получить столбцовый вектор x, вычисляем:
        Eigen::VectorXd x = M.transpose() * c;

        // Вычисляем коэффициенты с использованием функции аппроксимации.
        Eigen::VectorXd b = approximate_with_non_orthogonal_basis_orto(x, M);

        // Считаем ошибку (норма разности)
        double error = (b - c).norm();
        std::cout << "dim " << n << ": err = " << error;
        if (error < tol) {
            std::cout << " [PASSED]\n";
        }
        else {
            std::cout << " [FAILED]\n";
            all_passed = false;
        }
    }

    if (all_passed) {
        std::cout << "OK.\n";
    }
    else {
        std::cout << "FAIL.\n";
    }
}

#ifdef UNIT_TESTS
int main() {
    run_tests();
    return 0;
}
#else
int main() {

    // Запускаем тесты, если необходимо
    run_tests();
    // Параметры проекта
    std::string root_folder = "C:/dmitrienkomy/cache/";
    std::string bath = "y_200_2000";
    std::string wave = "gaus_single_2_h";
    std::string basis = "basis_6";

    // Инициализация конфигурации области (файл zones.json должен быть корректным)
    AreaConfigurationInfo area_config("T:/tsunami_res_folder/info/zones.json");

    // Вычисляем и сохраняем статистику аппроксимации
    save_and_plot_statistics(root_folder, bath, wave, basis, area_config);
    return 0;
}
#endif
