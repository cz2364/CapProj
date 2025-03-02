#include "defs.h"
#include<string>
#include<vector>

StdSoc::StdSoc(std::string name, int dim, std::vector<std::string>& var_names, std::vector<long>& var_index){
    cone_name = name;
    cone_dim = dim;
    cone_var_names = std::vector<std::string>(dim);
    std::copy(var_names.begin(), var_names.end(), cone_var_names.begin());
    cone_var_index = std::vector<long>(dim);
    std::copy(var_index.begin(), var_index.end(), cone_var_index.begin());
}
std::string StdSoc::get_name(){
    return cone_name;
}
int StdSoc::get_dim(){
    return cone_dim;
}
