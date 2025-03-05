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

std::vector<double> box_proj(const std::vector<double>&x, std::vector<double>& ub, std::vector<double>& lb){
    std::vector<double> px(x);
    int n = x.size();
    for(int i = 0; i < n; i++){
        double t = (x[i] > ub[i]) ? ub[i]:x[i];
        px[i] = (t < lb[i]) ? lb[i]:t;
    }
    return px;
}