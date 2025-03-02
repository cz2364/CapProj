#ifndef DEFS_H
#define DEFS_H

#include <string>
#include <vector>

class StdSoc{
    public:
        StdSoc();
        StdSoc(std::string name, 
               int dim, 
               std::vector<std::string>& var_name, 
               std::vector<long>& var_index);
        
        std::string get_name();
        int get_dim();

    private:
        std::string cone_name;
        int cone_dim;
        std::vector<std::string> cone_var_names;
        std::vector<long> cone_var_index;

};
#endif