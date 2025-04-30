#ifndef DEFS_H
#define DEFS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <limits>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat; // declares a row-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


enum class Domain{
    leq,
    eq,
    geq,
    incone,
};

enum class Sense{
    minimize,
    maximize,
};



struct SolverParams{
    int max_iter;
    double tol;
    double rho;
    double alpha;
    double beta;
    double gamma;
};

class Solver{
    // assume the problem has been cleaned up to the following form 
    // min c^\top x + t_1 + \ldots + t_p
    // s.t. l_v \leq x \leq u_v   # l_v can be -inf and u_v can be inf
    //      A1x <= b1   # A1 is assume to be a sparse matrix (eigen triplet) 
    //      A2x >= b2   # A2 is assume to be a sparse matrix (eigen triplet) 
    //      A3x = b3   # the equaility constraints
    //      Fx = u   # coupling the x linear variables and the u conic variables 
    //      (t, u) \in K_i # the standard cones 

    // the problem can be reformulated as the following saddle point problem
    // min_{l_v \leq x \leq u_v, (t, u) \in K_i} \max_{p1, p2 \geq 0, p3, q free} 
    //      c^\top x + t_1 + \ldots + t_p + p1^\top(A1x - b1) - p2^\top (A2x - b2) + p3^\top (A3x - b3) + q^\top (Fx - u)

    // need to input the dim(p1), dim(p2), dim(p3), and dim(q) 
      
    
    public:
        Solver(SpMat c, SpMat A, Eigen::VectorXd b, 
               std::vector<std::vector<long> > cones_var_indices,
               std::unordered_map<long, std::pair<double, double> > primal_proj_bounds,
               std::unordered_map<long, std::pair<double, double> > dual_proj_bounds,
               long num_primal_var,
               long num_cone_var,
               long num_dual_var,
               SolverParams params);
        Solver(const Solver& solver) = delete;
        Solver& operator=(const Solver& solver) = delete;
        void solve();
        std::pair<double, double> onePDHGStep(Eigen::VectorXd& p_n_c, 
                                              Eigen::VectorXd& p_n_c_new,
                                              Eigen::VectorXd& d,
                                              Eigen::VectorXd& d_new, 
                                              double primal_weight, 
                                              double step_size, 
                                              long n_outer_iter);
        void computeRestartCandidate();
        void primalWeightUpdate();
    private:
        SpMat c; // 1 * (np + nc)
        SpMat A; // nd * (np + nc)
        Eigen::VectorXd b; // nd * 1
        std::vector<std::vector<long> > cones_var_indices; 
        std::unordered_map<long, std::pair<double, double> > primal_proj_bounds; // np
        std::unordered_map<long, std::pair<double, double> > dual_proj_bounds;   //nd
        long num_primal_var, num_cone_var, num_dual_var; 
        SolverParams params;
        Eigen::VectorXd var_p_n_c;
        Eigen::VectorXd var_d;

};


class Task{
    class VariablesManager{
        // every variable must be registered at the variables manager
        // an index will be assigned to a variable when registering
        public:
            VariablesManager();
            VariablesManager(const VariablesManager& vars) = delete;
            VariablesManager& operator=(const VariablesManager& vars) = delete;
            void add_primal(std::string var_name, double lb, double ub);
            void add_dual(std::string var_name, Domain sense);
            void add_duals(std::vector<std::string> var_names, Domain sense);
            void add_duals(std::vector<std::string> var_names, std::vector<Domain> senses);
            void add_cone(std::string var_name, double lb, double ub);

            const std::vector<std::string>&  get_primal_variable_names() const;
            const std::vector<long>&  get_primal_variable_indices() const;
            const std::vector<double>& get_primal_variable_values() const;
            const std::unordered_map<std::string, std::pair<long, double> >& get_primal_variable_name_index_value() const;
            const std::vector<double>& get_primal_variable_lower_bounds() const;
            const std::vector<double>& get_primal_variable_upper_bounds() const;
            const std::unordered_map<std::string, std::pair<long, std::pair<double, double> > >& get_primal_variable_name_index_bounds() const;
            const std::vector<std::string>& get_dual_variable_names() const;
            const std::vector<long>& get_dual_variable_indices() const;
            const std::vector<double>& get_dual_variable_values() const;
            const std::unordered_map<std::string, std::pair<long, double> >& get_dual_variable_name_index_value() const;
            const std::unordered_map<std::string, std::pair<long, std::pair<double, double> > >& get_dual_variable_name_index_bounds() const;
            const std::unordered_map<std::string, std::pair<long, double> >& get_cone_variable_name_index_value() const;
            const std::unordered_map<std::string, std::pair<long, std::pair<double, double> > >& get_cone_variable_name_index_bounds() const;

            long get_num_primal_variables();
            long get_num_dual_variables();
            long get_num_cone_variables();

        private:
            std::vector<std::string> primal_variable_names;
            std::vector<long> primal_variable_indices;
            std::vector<double> primal_variable_values;
            std::unordered_map<std::string, std::pair<long, double> > primal_variable_name_index_value;
            
            std::vector<double> primal_variable_lower_bounds;
            std::vector<double> primal_variable_upper_bounds;
            std::unordered_map<std::string, std::pair<long, std::pair<double, double> > > primal_variable_name_index_bounds;


            // register when a constraint is added to the task
            std::vector<std::string> dual_variable_names;
            std::vector<long> dual_variable_indices;
            std::vector<double> dual_variable_values;
            std::unordered_map<std::string, std::pair<long, double> > dual_variable_name_index_value;
            std::unordered_map<std::string, std::pair<long, std::pair<double, double> > > dual_variable_name_index_bounds;

            // register cone variables when a cone is added to the task
            std::vector<std::string> cone_variable_names;
            std::vector<long> cone_variable_indices;
            std::vector<double> cone_variable_values;
            std::unordered_map<std::string, std::pair<long, double> > cone_variable_name_index_value;
            std::unordered_map<std::string, std::pair<long, std::pair<double, double> > > cone_variable_name_index_bounds;


            long num_primal_variables;
            long num_dual_variables;
            long num_cone_variables;

    };
    class ConstraintsManager{
        // every constraint must be registered at the constraints manager
        // an index will be assigned to a constraint when registering
        // track dual variables
        
        public:
            ConstraintsManager();
            void add_linear_constraint(std::string constr_name, 
                                std::vector<std::string> constr_vars, 
                                std::vector<double> a, 
                                double b,
                                Domain sense, 
                                VariablesManager& vars);
            void add_linear_constraints(std::vector<std::string> constr_names, 
                                        std::vector<std::vector<std::string> > constr_vars_lists, 
                                        std::vector<std::vector<double> >& A, 
                                        std::vector<double> b,
                                        std::vector<Domain> senses,
                                        VariablesManager& vars);
            void add_coupling_constraint(std::string constr_name,
                                       std::vector<std::string> constr_linear_vars,
                                       std::vector<double> a,
                                       std::string cone_var,
                                       VariablesManager& vars);
            void add_coupling_constraints(std::vector<std::string> constr_names,
                                          std::vector<std::vector<std::string> > constr_linear_vars,
                                          std::vector<std::vector<double> > A,
                                          std::vector<std::string> cone_vars,
                                          VariablesManager& vars); 
            std::set<std::string> get_constraints_linear_var_names_set();
            long get_num_constraints();
            long get_num_nnz();
            const std::vector<std::vector<std::string> >& get_constraints_vars_lists() const;
            const std::vector<double>& get_constraints_b() const;
            const std::vector<std::vector<double> > & get_constraints_A() const;
        private:
            std::vector<std::string> constraints_names;
            std::vector<long> constraints_indices;
            std::vector<std::vector<double> > constraints_A;
            std::vector<double> constraints_b;
            std::vector<Domain> constraints_senses;
            std::vector<std::vector<std::string> > constraints_vars_lists;
            std::vector<long> constr_indices;
            std::set<std::string> constraints_linear_var_names_set;
            std::set<std::string> constraints_cone_var_names_set;
            long num_constraints;
            long num_nnz;
    };
    class ConesManager{
        // every second order cone must be registered at the constraints manager
        // an index will be assigned to a constraint when registering
        public:
            ConesManager();
            void add_cone(std::string cone_name, 
                          std::vector<std::string> var_names, 
                          VariablesManager& vars,
                          ConstraintsManager& cons);
            void add_rotated_cone(std::string cone_name, 
                                  std::vector<std::string> var_names, 
                                  VariablesManager& vars,
                                  ConstraintsManager& cons);
            long get_num_cones();
            const std::vector<std::vector<std::string> >& get_cones_var_names() const;

        private:
            std::vector<std::string> cones_names;
            std::vector<long> cones_indices;
            std::vector<std::vector<std::string> > cones_var_names;
            std::vector<std::vector<std::string> > cones_valid_var_names;
            std::set<std::string> cones_var_names_set;
            long num_cones;
    };
    class ObjManager{
        public:
            ObjManager();
            void set_obj_sense(Sense sense);
            void set_obj_coeff(std::string var_name, double coeff, VariablesManager& vars);
            void set_obj_coeffs(std::unordered_map<std::string, double> coeffs, VariablesManager& vars);
            void set_obj_coeffs(std::vector<std::string> var_names, std::vector<double> coeffs, VariablesManager& vars);
            const std::unordered_map<std::string, double> & get_obj_coeffs() const;
        private:
            Sense obj_sense;
            std::unordered_map<std::string, double> obj_coeffs;
    };
    

    public:
        Task();
        void solve();
        long find_variable_sol_index(std::string vn);

    private:
        VariablesManager vars;
        ConstraintsManager cons;
        ConesManager cones;
        ObjManager objs;
    
};

#endif