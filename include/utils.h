#ifndef UTILS_H
#define UTILS_H
#include<vector>
#include<string>
#include "defs.h"

std::vector<double> std_soc_proj(const std::vector<double>& x);

void std_soc_proj_inplace(std::vector<double>& x);

std::vector<double> box_proj(const std::vector<double>&x, 
                             const std::vector<double>& ub, 
                             const std::vector<double>& lb);

void box_proj_inplace(std::vector<double>&x, 
                      const std::vector<double>& ub,
                      const std::vector<double>& lb);

std::string generate_variable_name(std::string x);

#endif