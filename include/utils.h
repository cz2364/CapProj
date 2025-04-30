#ifndef UTILS_H
#define UTILS_H
#include<vector>
#include<string>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "defs.h"


std::vector<double> std_soc_proj(const std::vector<double>& x);

void std_soc_proj_inplace(std::vector<double>& x);

Eigen::VectorXd std_soc_proj(Eigen::VectorXd & x);

std::vector<double> box_proj(const std::vector<double>&x, 
                             const std::vector<double>& ub, 
                             const std::vector<double>& lb);

void box_proj_inplace(std::vector<double>&x, 
                      const std::vector<double>& ub,
                      const std::vector<double>& lb);

Eigen::VectorXd box_proj(const Eigen::VectorXd x, 
                         const std::unordered_map<long, std::pair<double, double> > bounds);



std::string generate_variable_name(std::string x);


double weighted_norm(Eigen::VectorXd x, Eigen::VectorXd y, double w);


#endif