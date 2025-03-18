#include <iostream>
#include <vector>
#include <filesystem>
#include <Eigen/Dense>
#include "approx_orto.h"
#include "stable_data_structs.h"
#include "statistics.h"

// For convenience
namespace fs = std::filesystem;

// Function to copy a folder (recursively)
bool copyFolder(const std::string& source, const std::string& destination) {
    try {
        fs::create_directories(destination);
        for (const auto& entry : fs::recursive_directory_iterator(source)) {
            const auto& path = entry.path();
            auto relativePath = fs::relative(path, source);
            fs::copy(path, fs::path(destination) / relativePath, fs::copy_options::recursive | fs::copy_options::overwrite_existing);
        }
    }
    catch (fs::filesystem_error& e) {
        std::cerr << "Folder copy error: " << e.what() << std::endl;
        return false;
    }
    return true;
}

// Function to delete a folder
bool deleteFolder(const std::string& folder) {
    try {
        fs::remove_all(folder);
    }
    catch (fs::filesystem_error& e) {
        std::cerr << "Folder deletion error: " << e.what() << std::endl;
        return false;
    }
    return true;
}

// Function to copy a file
bool copyFile(const std::string& source, const std::string& destination) {
    try {
        fs::create_directories(fs::path(destination).parent_path());
        fs::copy_file(source, destination, fs::copy_options::overwrite_existing);
    }
    catch (fs::filesystem_error& e) {
        std::cerr << "File copy error: " << e.what() << std::endl;
        return false;
    }
    return true;
}

// Function to delete a file
bool deleteFile(const std::string& file) {
    try {
        fs::remove(file);
    }
    catch (fs::filesystem_error& e) {
        std::cerr << "File deletion error: " << e.what() << std::endl;
        return false;
    }
    return true;
}

// Function to check if a file or folder exists
bool fileExists(const std::string& file) {
    return fs::exists(file);
}

// It is assumed that AreaConfigurationInfo and the function save_and_plot_statistics are already defined
// For example:
// class AreaConfigurationInfo { /* ... */ };
// void save_and_plot_statistics(const std::string&, const std::string&, const std::string&, const std::string&, const AreaConfigurationInfo&);

int runWithPrePost(const std::string& root_folder,
    const std::string& cache_folder,
    const std::string& bath,
    const std::string& wave,
    const std::string& basis,
    const AreaConfigurationInfo& area_config) {
    // Build paths for copying the basis folder
    std::string sourceBasisFolder = root_folder + "/" + bath + "/" + basis;
    std::string destBathFolder = cache_folder + "/" + bath;
    std::string destBasisFolder = destBathFolder + "/" + basis;

    // Copy the basis folder
    bool copiedFolder = copyFolder(sourceBasisFolder, destBasisFolder);
    if (!copiedFolder) {
        std::cerr << "Failed to copy folder: " << sourceBasisFolder << std::endl;
        return -1;
    }

    // Build paths for the wave folder
    std::string sourceWaveFolder = root_folder + "/" + bath + "/" + wave + ".bp";
    std::string destWaveFolder = destBathFolder + "/" + wave + ".bp";

    bool folderAlreadyExists = fileExists(destWaveFolder);
    bool copiedWaveFolder = false;

    // If the folder does not exist in the destination and the source folder exists, copy it
    if (!folderAlreadyExists && fs::exists(sourceWaveFolder)) {
        copiedWaveFolder = copyFolder(sourceWaveFolder, destWaveFolder);
        if (!copiedWaveFolder) {
            std::cerr << "Failed to copy folder: " << sourceWaveFolder << std::endl;
            // If copying the folder failed, delete the previously copied basis folder
            deleteFolder(destBasisFolder);
            return -1;
        }
    }

    // Execute the main function
    save_and_plot_statistics(cache_folder, bath, wave, basis, area_config);

    // Delete the copied basis folder from the cache
    if (!deleteFolder(destBasisFolder)) {
        std::cerr << "Failed to delete folder: " << destBasisFolder << std::endl;
    }

    // If the wave folder was copied (i.e. it did not exist beforehand) – delete it
    if (copiedWaveFolder) {
        if (!deleteFolder(destWaveFolder)) {
            std::cerr << "Failed to delete folder: " << destWaveFolder << std::endl;
        }
    }

    return 0;
}

// Function to run tests
void run_tests() {
    // Dimensions for test cases
    std::vector<int> dimensions = { 3, 4, 5, 6, 8 };
    // Tolerance (for example, 1e-6)
    double tol = 1e-6;
    bool all_passed = true;

    std::cout << "tests (approximate_with_non_orthogonal_basis_orto):\n";

    for (int n : dimensions) {
        // Generate a random square matrix n x n,
        // representing the basis (each row is a basis vector).
        Eigen::MatrixXd M;
        // Ensure invertibility (regenerate if the determinant is too small)
        do {
            M = Eigen::MatrixXd::Random(n, n);
        } while (std::abs(M.determinant()) < 1e-3);

        // Generate a random coefficient vector c of length n.
        Eigen::VectorXd c = Eigen::VectorXd::Random(n);
        // Compute the vector x as a linear combination of basis vectors:
        // x = c[0]*M.row(0) + ... + c[n-1]*M.row(n-1).
        // To get a column vector x, compute:
        Eigen::VectorXd x = M.transpose() * c;

        // Compute the coefficients using the approximation function.
        Eigen::VectorXd b = approximate_with_non_orthogonal_basis_orto(x, M);

        // Calculate the error (norm of the difference)
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
int main(int argc, char* argv[]) {
    // Check for required arguments: bath, wave, basis
    if (argc < 4) {
        std::cerr << "usage: " << argv[0] << " bath wave basis" << std::endl;
        return 1;
    }

    // Read command line arguments
    std::string bath = argv[1];
    std::string wave = argv[2];
    std::string basis = argv[3];

    // Other parameters can be fixed or also obtained from arguments
    std::string root_folder = "T:/tsunami_adios_res";
    std::string cache_folder = "C:/dmitrienkomy/cache/";
<<<<<<< Updated upstream
    std::string bath = "parabola_200_2000";
    std::string wave = "gaus_single_2_h";
    std::string basis = "basis_48";
    std::vector<std::string> folderNames = {
        /*"basis_6",
        "basis_8",
        "basis_9",
        "basis_10",
        "basis_12",
        "basis_15",
        "basis_16",
        "basis_18",
        "basis_20",*/
        "basis_24",
        "basis_25",
        "basis_30",
        "basis_36",
        "basis_40",
        "basis_48"
    };
=======
>>>>>>> Stashed changes

    // Initialize area configuration (zones.json must be correct)
    AreaConfigurationInfo area_config("T:/tsunami_res_folder/info/zones.json");

<<<<<<< Updated upstream
    // Вычисляем и сохраняем статистику аппроксимации
    for (auto& basis : folderNames)
    {
        runWithPrePost(root_folder, cache_folder, bath, wave, basis, area_config);
=======
    // If needed, run run_tests(), for example for debugging:
    // run_tests();

    // Execute processing for given parameters
    // You can use either runWithPrePost or directly save_and_plot_statistics
    // Example usage of runWithPrePost:
    if (runWithPrePost(root_folder, cache_folder, bath, wave, basis, area_config) != 0) {
        std::cerr << "run error in runWithPrePost.\n";
        return 1;
>>>>>>> Stashed changes
    }
    //save_and_plot_statistics(cache_folder, "x_200_2000", "gaus_double_1_2", "basis_6", area_config);
    return 0;
}
#endif
