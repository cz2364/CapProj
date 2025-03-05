#ifndef DEFS_H
#define DEFS_H

#include <string>
#include <vector>
#include <Eigen/Core>

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

enum class Domain{
    leq,
    eq,
    geq,
    incone,
};

class Variable{
    public:
        Variable();
        Variable(std::string name, double lb, double ub);
        double get_lower();
        double get_upper();
    private:
        std::string variable_name;
        double variable_lower_bound;
        double variable_upper_bound;
        bool is_valid();
};

class Constraint {
    public:
        Constraint();
        Constraint(std::string name, 
                   std::vector<std::string>& var_names, 
                   const Eigen::MatrixXf& A, 
                   const Eigen::MatrixXf& b, 
                   Domain sense);
    private:
        std::string constr_name;

};


class Task{
    class VariablesManager{
        // every variable must be registered at the variables manager
        // an index will be assigned to a variable when registering
    };
    class ConstraintsManager{
        // every constraint must be registered at the constraints manager
        // an index will be assigned to a constraint when registering
        // track dual variables
    };
};

#endif