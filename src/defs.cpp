#include "defs.h"
#include "utils.h"
#include<string>
#include<vector>
#include <cmath>
#include "absl/log/log.h"


Task::VariablesManager::VariablesManager(){
    primal_variable_names = std::vector<std::string>();
    primal_variable_indices = std::vector<long>();
    primal_variable_values = std::vector<double>();
    primal_variable_name_index_value = std::unordered_map<std::string, std::pair<long, double> >();
    dual_variable_names = std::vector<std::string>();
    dual_variable_indices = std::vector<long>();
    dual_variable_values = std::vector<double>();
    dual_variable_name_index_value = std::unordered_map<std::string, std::pair<long, double> >();
    num_primal_variables = 0;
    num_dual_variables = 0;
}

void Task::VariablesManager::add_primal(std::string var_name, double lb, double ub){
    primal_variable_names.push_back(var_name);
    primal_variable_indices.push_back(num_primal_variables);
    primal_variable_values.push_back(0.0);
    primal_variable_name_index_value[var_name] = std::make_pair(num_primal_variables, 0.0);

    primal_variable_lower_bounds.push_back(lb);
    primal_variable_upper_bounds.push_back(ub);
    primal_variable_name_index_bounds[var_name] = std::make_pair(num_primal_variables, std::make_pair(lb, ub));
    num_primal_variables += 1;
}

void Task::VariablesManager::add_dual(std::string var_name, Domain sense){
    // need to add a logger for record the exception.
    dual_variable_names.push_back(var_name);
    dual_variable_indices.push_back(num_dual_variables);
    dual_variable_values.push_back(0.0);
    dual_variable_name_index_value[var_name] = std::make_pair(num_dual_variables, 0.0);
    if(sense == Domain::leq){
        dual_variable_name_index_bounds[var_name] = std::make_pair(num_dual_variables, 
                                                                   std::make_pair(0.0, 
                                                                                  std::numeric_limits<double>::infinity()));
        num_dual_variables += 1;                            
    }else if(sense == Domain::geq) {
        dual_variable_name_index_bounds[var_name] = std::make_pair(num_dual_variables, 
                                                                   std::make_pair(-std::numeric_limits<double>::infinity(), 
                                                                                  0.0));
        num_dual_variables += 1;                                                                          
    }else if(sense == Domain::eq){
        dual_variable_name_index_bounds[var_name] = std::make_pair(num_dual_variables, 
                                                                std::make_pair(-std::numeric_limits<double>::infinity(), 
                                                                               std::numeric_limits<double>::infinity()));
        num_dual_variables += 1;
    }else{
        // add exceptions log here
    }

    
}

void Task::VariablesManager::add_duals(std::vector<std::string> var_names, Domain sense) {
    for(int i = 0; i < var_names.size(); i++){
        add_dual(var_names[i], sense);
    }
}

void Task::VariablesManager::add_duals(std::vector<std::string> var_names, std::vector<Domain> senses){
    if (var_names.size() != senses.size()){
        // terminate and long
    }else{
        for(int i = 0; i < var_names.size(); i++){
            add_dual(var_names[i], senses[i]);
        }
    }
}

void Task::VariablesManager::add_cone(std::string var_name, double lb, double ub){
    cone_variable_names.push_back(var_name);
    cone_variable_indices.push_back(num_cone_variables);
    cone_variable_values.push_back(0.0);
    cone_variable_name_index_value[var_name] = std::make_pair(num_cone_variables, 0.0);
    cone_variable_name_index_bounds[var_name] = std::make_pair(num_cone_variables, std::make_pair(lb, ub));
    num_cone_variables += 1;

}

const std::unordered_map<std::string, std::pair<long, double> > & Task::VariablesManager:: get_primal_variable_name_index_value() const{
    return primal_variable_name_index_value;
}

const std::unordered_map<std::string, std::pair<long, std::pair<double, double> > > & Task::VariablesManager::get_primal_variable_name_index_bounds() const{
    return primal_variable_name_index_bounds;
}


const std::unordered_map<std::string, std::pair<long, double> > & Task::VariablesManager::get_cone_variable_name_index_value() const {
    return cone_variable_name_index_value;
}

const std::unordered_map<std::string, std::pair<long, std::pair<double, double> > > & Task::VariablesManager:: get_cone_variable_name_index_bounds() const {
    return cone_variable_name_index_bounds;
}

const std::unordered_map<std::string, std::pair<long, std::pair<double, double> > > & Task::VariablesManager::get_dual_variable_name_index_bounds() const {
    return dual_variable_name_index_bounds;
}

long Task::VariablesManager::get_num_primal_variables(){
    return num_primal_variables;
}

long Task::VariablesManager::get_num_dual_variables(){
    return num_dual_variables;
}

long Task::VariablesManager::get_num_cone_variables(){
    return num_cone_variables;
}


Task::ConstraintsManager::ConstraintsManager(){
    constraints_names = std::vector<std::string>();
    constraints_indices = std::vector<long>();
    constraints_A = std::vector<std::vector<double> >();
    constraints_b = std::vector<double>();
    constraints_vars_lists = std::vector<std::vector<std::string> >();
    constraints_senses = std::vector<Domain>();
    constraints_linear_var_names_set = std::set<std::string>();
    num_constraints = 0;
    num_nnz = 0;
}



std::set<std::string> Task::ConstraintsManager::get_constraints_linear_var_names_set(){
    return constraints_linear_var_names_set;
}

void Task::ConstraintsManager::add_linear_constraint(std::string constr_name, 
                                                     std::vector<std::string> constr_vars,
                                                     std::vector<double> a, 
                                                     double b, 
                                                     Domain sense,
                                                     VariablesManager& vars){
        if (constr_vars.size() != a.size()){
            // terminate and log
        }
        for(int i = 0; i < constr_vars.size(); i++){
            if(vars.get_primal_variable_name_index_value().find(constr_vars[i]) == vars.get_primal_variable_name_index_value().end()){
                // A variable in the constraint is not registered, need to register before using it.
                // terminate and log
            }
            constraints_linear_var_names_set.insert(constr_vars[i]);
        }
        constraints_names.push_back(constr_name);
        constraints_vars_lists.push_back(constr_vars);
        constraints_A.push_back(a);
        constraints_b.push_back(b);
        constraints_senses.push_back(sense);
        constraints_indices.push_back(num_constraints);
        num_constraints += 1;
        num_nnz += a.size();

        // add the corresponding dual variable
        vars.add_dual("dual_" + constr_name, sense);

    }

void Task::ConstraintsManager::add_linear_constraints(std::vector<std::string> constr_names, 
                                                      std::vector<std::vector<std::string> > constr_vars_lists, 
                                                      std::vector<std::vector<double>>& A, 
                                                      std::vector<double> b,
                                                      std::vector<Domain> senses, 
                                                      VariablesManager& vars){
        // check if the size of vectors agree
        if(constr_names.size()!=constr_vars_lists.size() || constr_names.size()!=A.size() ){
            // terminate and log
        }
        if (constr_names.size() != A.size() || A.size() != b.size() || b.size() != senses.size()){
            // terminate and log
        }
        for(int i = 0; i < constr_names.size(); i++){
            add_linear_constraint(constr_names[i], constr_vars_lists[i], A[i], b[i], senses[i], vars);
        }

    }

void Task::ConstraintsManager::add_coupling_constraint(std::string constr_name,
                                                       std::vector<std::string> constr_linear_vars,
                                                       std::vector<double> a,
                                                       std::string cone_var,
                                                       VariablesManager& vars) {

    // register the variables
    if (constr_linear_vars.size() !=  a.size()){
        // log and terminate
    }
    if(vars.get_cone_variable_name_index_value().find(cone_var) == vars.get_cone_variable_name_index_value().end()){
        // register the cone variable before use
        // log and terminate
    }
    for(int i = 0; i < constr_linear_vars.size(); i++){
        if(vars.get_primal_variable_name_index_value().find(constr_linear_vars[i]) == vars.get_primal_variable_name_index_value().end()){
            // regiter the linear variable before use
            // log and terminaate
        }
    }
    std::vector<std::string> temp_vars_list(constr_linear_vars);
    temp_vars_list.push_back(cone_var);
    std::vector<double> temp_a(a);
    temp_a.push_back(-1.0);
    double temp_b = 0.0;
    Domain sense = Domain::eq;

    constraints_names.push_back(constr_name);
    constraints_vars_lists.push_back(temp_vars_list);
    constraints_A.push_back(temp_a);
    constraints_b.push_back(temp_b);
    constraints_senses.push_back(sense);
    constraints_indices.push_back(num_constraints);
    num_constraints += 1;
    num_nnz += temp_a.size();

    // add the corresponding dual variable
    vars.add_dual("dual_" + constr_name, sense);
}

void Task::ConstraintsManager::add_coupling_constraints(std::vector<std::string> constr_names, 
                                                        std::vector<std::vector<std::string> > constr_linear_vars, 
                                                        std::vector<std::vector<double> > A, 
                                                        std::vector<std::string> cone_vars,
                                                        VariablesManager& vars) {
    // check sizes of inputs
    if(constr_names.size() != constr_linear_vars.size() || constr_names.size() != A.size()){
        // terminate and log
    }
    if(constr_names.size() != cone_vars.size()){
        // terminate and log
    }

    for(int i = 0; i < constr_names.size(); i++){
        add_coupling_constraint(constr_names[i], constr_linear_vars[i], A[i], cone_vars[i], vars);
    }

}



long Task::ConstraintsManager::get_num_constraints(){
    return num_constraints;
}

long Task::ConstraintsManager::get_num_nnz(){
    return num_nnz;
}
std::vector<double> const & Task::ConstraintsManager::get_constraints_b() const{
    return constraints_b;
}
std::vector<std::vector<std::string> > const & Task::ConstraintsManager::get_constraints_vars_lists() const{
    return constraints_vars_lists;
}

const std::vector<std::vector<double> > & Task::ConstraintsManager::get_constraints_A() const {
    return constraints_A;
}


Task::ConesManager::ConesManager(){
    cones_names = std::vector<std::string>();
    cones_indices = std::vector<long>();
    cones_var_names = std::vector<std::vector<std::string> >();
    cones_valid_var_names = std::vector<std::vector<std::string> >();
    cones_var_names_set = std::set<std::string>();
    num_cones = 0;
}


void Task::ConesManager::add_cone(std::string cone_name, 
                                  std::vector<std::string> var_names, 
                                  VariablesManager& vars,
                                  ConstraintsManager& cons){
    if (var_names.size() < 2) {
        // terminate and log
    }

    // check if the variables are registered at primal, if so, terminate and log
    for(int i = 0; i < var_names.size(); i++){
        if(vars.get_primal_variable_name_index_value().find(var_names[i]) != vars.get_primal_variable_name_index_value().end()){
            // terminate and log
        }
    }

    
    // check if the variables are registered, if not, register at the variables manager  
    for(int i = 0; i < var_names.size(); i++){
        if(vars.get_cone_variable_name_index_value().find(var_names[i]) == vars.get_cone_variable_name_index_value().end()){
            if(i == 0){
                vars.add_cone(var_names[i], 0.0, std::numeric_limits<double>::infinity());
            }else{
                vars.add_cone(var_names[i], 
                                -std::numeric_limits<double>::infinity(), 
                                std::numeric_limits<double>::infinity());
            }
        }
    }

    // check if there are existing variables in other cones, if so create a new variable and add a constraint
    std::vector<std::string> new_valid_names;
    for(int i = 0; i < var_names.size(); i++){
        if(cones_var_names_set.find(var_names[i]) != cones_var_names_set.end()){
            std::string new_var_name = generate_variable_name(var_names[i]);
            new_valid_names.push_back(new_var_name);
            cones_var_names_set.insert(new_var_name);
            // register the new variable at the variables manager
            std::pair<double, double> temp_bound = vars.get_cone_variable_name_index_bounds().at(var_names[i]).second;
            vars.add_cone(new_var_name, temp_bound.first, temp_bound.second);

            // add a constraint to the constraints manager
            std::vector<std::string> temp_vars_list = {new_var_name, var_names[i]};
            std::vector<double> temp_a = {1.0 ,-1.0};
            double temp_b = 0.0;
            cons.add_linear_constraint("arti_eq_constraint_" + new_var_name + "_" + var_names[i], 
                                       temp_vars_list, 
                                       temp_a, 
                                       temp_b, 
                                       Domain::eq, 
                                       vars);
        
        }else{
            new_valid_names.push_back(var_names[i]);
            cones_var_names_set.insert(var_names[i]);
        }
    }
    cones_names.push_back(cone_name);
    cones_indices.push_back(num_cones);
    cones_var_names.push_back(var_names);
    cones_valid_var_names.push_back(new_valid_names);
    num_cones += 1;
}

long Task::ConesManager::get_num_cones(){
    return num_cones;
}
const std::vector<std::vector<std::string> >& Task::ConesManager::get_cones_var_names() const{
    return cones_var_names;
}

void::Task::ConesManager::add_rotated_cone(std::string cone_name, 
    std::vector<std::string> var_names, 
    VariablesManager& vars,
    ConstraintsManager& cons){
    // check if the variable one is registered
    if(vars.get_cone_variable_name_index_value().find("one") == vars.get_cone_variable_name_index_value().end()){
        vars.add_cone("one", 1.0, 1.0);
        cons.add_linear_constraint("one_eq", {"one"}, {1.0}, 1.0, Domain::eq, vars);
    }

    if (var_names.size() < 3) {
        // terminate and log
    }
    // check if the variables are registered, if not, register at the variables manager  
    for(int i = 0; i < var_names.size(); i++){
        if(vars.get_cone_variable_name_index_value().find(var_names[i]) == vars.get_cone_variable_name_index_value().end()){
            if(i == 0 || i == 1){
                vars.add_cone(var_names[i], 0.0, std::numeric_limits<double>::infinity());
            }else{
                    vars.add_cone(var_names[i], 
                                  -std::numeric_limits<double>::infinity(), 
                                  std::numeric_limits<double>::infinity());
            }
        }
    }

    // convert the rotated cone to a second order cone
    std::string temp_add_name = "converted_" + var_names[0] + "_add_" + var_names[1];
    std::string temp_sub_name = "converted_" + var_names[0] + "_sub_" + var_names[1];

    std::vector<std::string> converted_cone_var_names;
    converted_cone_var_names.push_back(temp_add_name);
    converted_cone_var_names.push_back(temp_sub_name);

    for(int i = 2; i < var_names.size(); i++){
        converted_cone_var_names.push_back(var_names[i]);
    }

    cons.add_linear_constraint("rotated_cone_reform_converted_" + var_names[0] + "_add_" + var_names[1],
                               {converted_cone_var_names[0], var_names[0], var_names[1]}, 
                               {1.0, -1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)}, 
                               0.0, Domain::eq, vars);
    cons.add_linear_constraint("rotated_cone_reform_converted_" + var_names[0] + "_sub_" + var_names[1],
                               {converted_cone_var_names[1], var_names[0], var_names[1]},
                               {1.0, -1.0 / std::sqrt(2.0), 1.0/ std::sqrt(2.0)},
                               0.0, Domain::eq, vars);


    // check if there are existing variables in other cones, if so create a new variable and add a constraint
    std::vector<std::string> new_valid_names;
    for(int i = 0; i < converted_cone_var_names.size(); i++){
        if(cones_var_names_set.find(converted_cone_var_names[i]) != cones_var_names_set.end()){
            std::string new_var_name = generate_variable_name(converted_cone_var_names[i]);
            new_valid_names.push_back(new_var_name);
            cones_var_names_set.insert(new_var_name);
            // register the new variable at the variables manager
            std::pair<double, double> temp_bound = vars.get_cone_variable_name_index_bounds().at(converted_cone_var_names[i]).second;
            vars.add_cone(new_var_name, temp_bound.first, temp_bound.second);

            // add a constraint to the constraints manager
            std::vector<std::string> temp_vars_list = {new_var_name, converted_cone_var_names[i]};
            std::vector<double> temp_a = {1.0 ,-1.0};
            double temp_b = 0.0;
            cons.add_linear_constraint("arti_eq_constraint_" + new_var_name + "_" + converted_cone_var_names[i], 
                                       temp_vars_list, 
                                       temp_a, 
                                       temp_b, 
                                       Domain::eq, 
                                       vars);
        
        }else{
            new_valid_names.push_back(converted_cone_var_names[i]);
            cones_var_names_set.insert(converted_cone_var_names[i]);
        }
    }
    cones_names.push_back(cone_name);
    cones_indices.push_back(num_cones);
    cones_var_names.push_back(var_names);
    cones_valid_var_names.push_back(new_valid_names);
    num_cones += 1;


}

Task::ObjManager::ObjManager(){
    obj_sense = Sense::minimize;
    obj_coeffs = std::unordered_map<std::string, double>();
}

void Task::ObjManager::set_obj_sense(Sense sense){
    obj_sense = sense;
}
void Task::ObjManager::set_obj_coeff(std::string var_name, double coeff, VariablesManager& vars){
    if(vars.get_primal_variable_name_index_value().find(var_name) == vars.get_primal_variable_name_index_value().end()){
        // terminate and log
    }
    obj_coeffs[var_name] = coeff;
}

void Task::ObjManager::set_obj_coeffs(std::unordered_map<std::string, double> coeffs, VariablesManager& vars){
    for(auto it = coeffs.begin(); it != coeffs.end(); ++it){
        if(vars.get_primal_variable_name_index_value().find(it->first) == vars.get_primal_variable_name_index_value().end()){
            // terminate and log
        }
        obj_coeffs[it->first] = it->second;
    }
}

const std::unordered_map<std::string, double> & Task::ObjManager::get_obj_coeffs() const{
    return obj_coeffs;
}


long Task::find_variable_sol_index(std::string vn){
    //first primal then cone last dual
    long num_primal_var =vars.get_num_primal_variables();
    long num_dual_var = vars.get_num_dual_variables();
    long num_cone_var = vars.get_num_cone_variables();
    if(vars.get_primal_variable_name_index_value().find(vn) != vars.get_primal_variable_name_index_value().end()){
        return vars.get_primal_variable_name_index_value().at(vn).first;
    }else if(vars.get_cone_variable_name_index_value().find(vn) != vars.get_cone_variable_name_index_value().end()){
        return num_primal_var + vars.get_cone_variable_name_index_value().at(vn).first;
    }else if(vars.get_dual_variable_name_index_value().find(vn) != vars.get_dual_variable_name_index_value().end()){
        return num_primal_var + num_cone_var + vars.get_dual_variable_name_index_value().at(vn).first;
    }else{
        // log and terminate
        return -1;
    }
    
}


Solver::Solver(SpMat c, SpMat A, Eigen::VectorXd b, 
    std::vector<std::vector<long> > cones_var_indices,
    std::unordered_map<long, std::pair<double, double> > primal_proj_bounds,
    std::unordered_map<long, std::pair<double, double> > dual_proj_bounds,
    long num_primal_var,
    long num_cone_var,
    long num_dual_var,
    SolverParams params){
        c = c;
        A = A;
        b = b;
        cones_var_indices = cones_var_indices;
        primal_proj_bounds = primal_proj_bounds;
        num_primal_var = num_primal_var;
        num_cone_var = num_cone_var;
        num_dual_var = num_dual_var;
        params = params;
    }

std::pair<double, double> Solver::onePDHGStep(Eigen::VectorXd& p_n_c, 
                                              Eigen::VectorXd& p_n_c_new,
                                              Eigen::VectorXd& d,
                                              Eigen::VectorXd& d_new,
                                              double primal_weight, 
                                              double step_size, 
                                              long n_outer_iter){
    double cur_ss = step_size;

    while(true){
        p_n_c_new = p_n_c - (cur_ss / primal_weight) * (c.transpose() - A.transpose() * d.transpose());
        Eigen::VectorXd p_temp = p_n_c_new.head(num_primal_var);
        Eigen::VectorXd p_new = box_proj(p_temp, primal_proj_bounds);
        p_n_c_new.segment(0, num_primal_var) = p_new;

        // cone projection

        for(int i = 0; i < cones_var_indices.size(); i++){
            Eigen::VectorXd temp_cone_var(cones_var_indices[i].size());
            for(int j = 0; j < cones_var_indices[i].size(); j++){
                temp_cone_var[j] = p_n_c_new[cones_var_indices[i][j]];
            }
            Eigen::VectorXd temp_cone_var_proj = std_soc_proj(temp_cone_var);
            for (int j = 0; j < cones_var_indices[i].size(); j++){
                p_n_c_new[cones_var_indices[i][j]] = temp_cone_var_proj[j];
            }
        }
        //
        Eigen::VectorXd d_temp_1 = 1.0 / (primal_weight * cur_ss) * d - A * (2.0 * p_n_c_new - p_n_c);
        Eigen::VectorXd d_temp_2 = box_proj(d_temp_1, dual_proj_bounds);
        d_new = d - A * (2.0 * p_n_c_new - p_n_c)  - (cur_ss * primal_weight) * d_temp_2;

        double test_ss_numerator = std::pow(weighted_norm(p_n_c_new - p_n_c, d_new - d, primal_weight), 2);
        double test_ss_denorm = 2 * (d_new.transpose() - d.transpose()) * A * (p_n_c_new - p_n_c);

        double test_ss = test_ss_numerator / test_ss_denorm;

        double nxt_ss = std::min((1 - std::pow(n_outer_iter + 1, -0.3)) * test_ss, 
                                 (1 + std::pow(n_outer_iter + 1, -0.6)) * cur_ss);
        if(cur_ss < test_ss){
            return std::make_pair(cur_ss, nxt_ss);
            
        }else{
            cur_ss = nxt_ss;
        }
    }
    

}




void Task::solve(){
    LOG(INFO) << "Preparing to solve the task";
    long num_primal_var =vars.get_num_primal_variables();
    long num_dual_var = vars.get_num_dual_variables();
    long num_cone_var = vars.get_num_cone_variables();

    std::vector<double> all_variables(num_primal_var + num_cone_var + num_dual_var, 0.0);
    // form the matrix A in sparse matrix
    std::vector<T> tripletList;
    tripletList.reserve(cons.get_num_nnz());
    
    const std::vector<std::vector<std::string> > & constr_vars_lists = cons.get_constraints_vars_lists();
    const std::vector<std::vector<double> > & constr_A = cons.get_constraints_A();

    for(int i = 0; i < cons.get_num_constraints(); i++){
        for (int j = 0; j < constr_vars_lists[i].size(); j++){
            long var_index = find_variable_sol_index(constr_vars_lists[i][j]);
            tripletList.push_back(T(i, var_index, constr_A[i][j]));
        }
    }
    SpMat A(cons.get_num_constraints(), num_primal_var + num_cone_var);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // form the b vector

    std::vector<T> tripletList_b;
    tripletList_b.reserve(cons.get_num_constraints());
    for(int i = 0; i < cons.get_constraints_b().size(); i++){
        tripletList_b.push_back(T(0, i, cons.get_constraints_b()[i]));
    }
    SpMat b(1, cons.get_num_constraints());

    // form the objective vector
    std::vector<T> tripletList_obj;
    tripletList.reserve(objs.get_obj_coeffs().size());
    for(auto it = objs.get_obj_coeffs().begin(); it != objs.get_obj_coeffs().end(); ++it){
        long var_index = find_variable_sol_index(it->first);
        tripletList_obj.push_back(T(0, var_index, it->second));
    }
    SpMat c(1, num_primal_var + num_cone_var + num_dual_var);
    c.setFromTriplets(tripletList_obj.begin(), tripletList_obj.end());

    // form list of cone variables indices
    std::vector<std::vector<long> > cones_var_indices;
    for(int i = 0; i < cones.get_num_cones(); i++){
        std::vector<long> temp_cone_var_indices;
        for(int j = 0; j < cones.get_cones_var_names()[i].size(); j++){
            long var_index = find_variable_sol_index(cones.get_cones_var_names().at(i)[j]);
            temp_cone_var_indices.push_back(var_index);
        }
        cones_var_indices.push_back(temp_cone_var_indices);
    }
    
    // form list of primal variable indices and bounds
    std::unordered_map<long, std::pair<double, double> > primal_proj_bounds;
    for(auto it = vars.get_primal_variable_name_index_bounds().begin(); it != vars.get_primal_variable_name_index_bounds().end(); ++it){
        long var_index = it->second.first;
        std::pair<double, double> bounds = it->second.second;
        primal_proj_bounds[var_index] = bounds;
    }

    // form list of dual variable indices and bounds
    std::unordered_map<long, std::pair<double, double> > dual_proj_bounds;
    for(auto it = vars.get_dual_variable_name_index_bounds().begin(); it != vars.get_dual_variable_name_index_bounds().end(); ++it){
        long var_index = it->second.first;
        std::pair<double, double> bounds = it->second.second;
        dual_proj_bounds[var_index] = bounds;
    }


    // initialize solution

    while(true){
        // Adpative step of PDHG
        // Get restart candidate

        // Primal weight update
        // check stopping criteria
    }


    

}