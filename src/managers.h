#ifndef MANAGERS_H
#define MANAGERS_H

#include <string>
#include <vector>
#include "stable_data_structs.h"

// Êëàññ äëÿ ðàáîòû ñ äàííûìè basis, ñîäåðæàùèìèñÿ â NetCDF-ôàéëàõ
class BasisManager {
public:
    std::string folder; // ïóòü ê êàòàëîãó ñ basis-ôàéëàìè (NetCDF-ôàéëû)

    explicit BasisManager(const std::string& folder_) : folder(folder_) {}

    // Ôóíêöèÿ ÷òåíèÿ äàííûõ basis äëÿ ðåãèîíà [y_start, y_end)
    // Âîçâðàùàåò 4D ìàññèâ: [num_files][T][region_height][X]
    std::vector<std::vector<std::vector<std::vector<double>>>> get_fk_region(int y_start, int y_end);
};

// Êëàññ äëÿ ðàáîòû ñ ìàðèîãðàììàìè (Wave data)
class WaveManager {
public:
    std::string nc_file; // ïóòü ê NetCDF-ôàéëó ñ ìàðèîãðàììàìè

    explicit WaveManager(const std::string& nc_file_) : nc_file(nc_file_) {}

    // Ôóíêöèÿ çàãðóçêè äàííûõ ïåðåìåííîé "height" äëÿ ðåãèîíà [y_start, y_end)
    // Âîçâðàùàåò 3D ìàññèâ: [T][region_height][X]
    std::vector<std::vector<std::vector<double>>> load_mariogramm_by_region(int y_start, int y_end);
};

#endif // MANAGERS_H