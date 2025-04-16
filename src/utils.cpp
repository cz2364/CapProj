#include "defs.h"
#include "utils.h"
#include<string>
#include<vector>


void box_proj_inplace(std::vector<double>&x, std::vector<double>& ub, std::vector<double>& lb){
    int n = x.size();
    for(int i = 0; i < n; i++){
        double t = (x[i] > ub[i]) ? ub[i]:x[i];
        x[i] = (t < lb[i]) ? lb[i]:t;
    }
}

std::vector<double> box_proj(const std::vector<double>& x, const std::vector<double>& ub, const std::vector<double>& lb){
    std::vector<double> px(x);
    int n = x.size();
    for(int i = 0; i < n; i++){
        double t = (x[i] > ub[i]) ? ub[i]:x[i];
        px[i] = (t < lb[i]) ? lb[i]:t;
    }
    return px;
}

std::string generate_variable_name(std::string x) {
    std::string pattern = "__COPY_";
    std::size_t found = x.rfind(pattern);
    if (found == std::string::npos){
        return x + pattern + "1";
    }else{
        std::string valid_name = x.substr(0, found);
        if (found + 7 == x.size()){
            return x + pattern + "1";
        }else{
            std::string suffix = x.substr(found + 7);
            try{
                int suffix_num = std::stoi(suffix);
                return valid_name + pattern + std::to_string(suffix_num + 1);
            } catch (const std::invalid_argument& e){
                return x + pattern + "1";
            }
        }
    }
}

std::vector<double> std_soc_proj(const std::vector<double>& x){
    int n = x.size();
    std::vector<double> result(n, 0.0);
    
    double norm = 0.0;
    for(int i = 1; i<n; i++){
        norm += x[i] * x[i];
    }
    norm = sqrt(norm);

    if (norm < -x[0]){
        return result;
    }else if(norm < x[0]) {
        for(int i = 0; i < n; i++){
            result[i] = x[i];
        }
        return result;
    }else{
        double multiplier = (1 + x[0] / norm) / 2.0;
        result[0] = multiplier * norm;
        for(int i = 1; i < n; i++){
            result[i] = multiplier * x[i];
        }
        return result;
    }

}