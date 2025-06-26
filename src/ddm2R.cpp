#include "ggdmcHeaders/ddm_simulation.h"
#include <RcppArmadillo.h>
#include <ggdmcHeaders/simulation_type_casting.h>

//' @rdname fastdm
//' @export
// [[Rcpp::export]]
std::vector<double> dfastdm(Rcpp::NumericVector rt_r,
                            Rcpp::NumericVector parameters_r,
                            bool is_lower = true)
{
    auto parameters = Rcpp::as<std::vector<double>>(parameters_r);
    auto ddm = ddm::ddm_class(parameters, is_lower);
    auto rt = Rcpp::as<std::vector<double>>(rt_r);
    auto out = ddm.dddm(rt);
    return out;
}

//' @rdname fastdm
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pfastdm(Rcpp::NumericVector rt_r,
                            Rcpp::NumericVector parameters_r,
                            bool is_lower = true)
{
    const size_t nrt = rt_r.size();
    auto parameters = Rcpp::as<std::vector<double>>(parameters_r);

    /* p_vector, passing also RT and precision
     * 0  1   2   3   4   5   6  7   8   9
     * a  v  zr   d  sz  sv  t0 st0  s   precision
     */
    // is_lower use g- density function.
    // ddm::ddm_class ddm;
    // ddm.set_parameters(parameters, is_lower);
    // ddm.print("Parameters: ");

    auto ddm = ddm::ddm_class(parameters, is_lower);
    auto f_ptr = select_f_calculator(ddm);
    const double scaled = ddm.m_zr * ddm.m_a;

    f_ptr->start(true); // Force set upper to calculate offset.
    const double offset = get_cdf(f_ptr, 0, scaled);

    // Determine direction and compute CDF values
    const bool use_upper = ddm.m_plus;
    f_ptr->start(use_upper);

    Rcpp::NumericVector out(nrt);
    for (size_t i = 0; i < nrt; ++i)
    {
        out[i] = use_upper ? (get_cdf(f_ptr, rt_r[i], scaled) - offset)
                           : (offset - get_cdf(f_ptr, rt_r[i], scaled));
    }

    return out;
}

//' @rdname fastdm
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame
rfastdm(unsigned int n = 1,
        Rcpp::NumericVector parameters_r = Rcpp::NumericVector::create(),
        Rcpp::NumericVector time_parameters_r = Rcpp::NumericVector::create(),
        bool debug = false)
{
    // Set default parameters if not provided
    if (parameters_r.size() == 0)
    {
        parameters_r = Rcpp::NumericVector::create(
            Rcpp::Named("a") = 1.0, Rcpp::Named("v") = 1.5,
            Rcpp::Named("zr") = 0.5, Rcpp::Named("d") = 0.0,
            Rcpp::Named("szr") = 0.0, Rcpp::Named("sv") = 0.0,
            Rcpp::Named("t0") = 0.2, Rcpp::Named("st0") = 0.0);
    }

    if (time_parameters_r.size() == 0)
    {
        time_parameters_r = Rcpp::NumericVector::create(
            Rcpp::Named("t_min") = -0.5, Rcpp::Named("t_max") = 0.5,
            Rcpp::Named("dt") = 0.01);
    }

    auto parameters = Rcpp::as<std::vector<double>>(parameters_r);
    auto time_parameters = Rcpp::as<std::vector<double>>(time_parameters_r);
    ddm::ddm_class ddm_obj;
    ddm_obj.set_parameters(parameters);
    ddm_obj.set_times(time_parameters);
    bool is_valid = ddm_obj.validate_parameters(debug);

    if (!is_valid && debug)
    {
        ddm_obj.print("simulate_ddm_trials");
        throw std::runtime_error("Invalid paraemters");
    }

    auto trials = ddm::simulate_trials(ddm_obj, n, debug);

    SimulationResults results(n);
    for (const auto &trial : trials)
    {
        results.add_trial(trial.first, trial.second);
    }

    return Rcpp::DataFrame::create(Rcpp::Named("RT") = results.reaction_times,
                                   Rcpp::Named("R") = results.responses);
}

//' @rdname simulate_ddm_trials
//' @export
// [[Rcpp::export]]
bool validate_ddm_parameters(const Rcpp::S4 &rt_model_r,
                             const Rcpp::NumericVector &parameters_r,
                             bool debug = false)
{
    auto d_ptr = new_design_light_rt_model(rt_model_r);
    size_t nparameter = parameters_r.size();
    if (d_ptr->m_n_free_parameter != nparameter)
    {
        std::string nparameter = std::to_string(d_ptr->m_n_free_parameter);
        Rcpp::stop("Invalid parameter vector: Expected length " + nparameter +
                   ", but got " + nparameter);
    }

    // auto l_ptr = new_likelihood_for_simulation(rt_model_r);
    auto parameters = Rcpp::as<std::vector<double>>(parameters_r);

    bool is_valid = false;
    ddm::ddm_class ddm_obj;

    for (size_t cell_idx = 0; cell_idx < d_ptr->m_n_cell; ++cell_idx)
    {
        d_ptr->set_parameter_values(cell_idx, parameters);
        auto cell_parameters = d_ptr->m_parameter_matrix[cell_idx];

        ddm_obj.set_parameters(cell_parameters);

        is_valid = ddm_obj.validate_parameters(debug);
        if (!is_valid)
        {
            Rcpp::Rcout << "Invalid DDM parameters at the condition, "
                        << d_ptr->m_cell_names[cell_idx] << ": \n";
            ddm_obj.print();
            break;
        }
    }
    return is_valid;
}

void simulate_conditions(const std::shared_ptr<design::design_class> &design,
                         ddm::ddm_class &ddm_obj,
                         const std::vector<double> &parameters,
                         unsigned int n_trial_per_cell,
                         SimulationResults &results, bool debug)
{
    std::string prev_condition;
    std::string node_1 = design->m_accumulator_names[0];

    for (size_t cell_idx = 0; cell_idx < design->m_n_cell; ++cell_idx)
    {
        const auto &cell_name = design->m_cell_names[cell_idx];
        const auto [condition, accumulator, valid] = parse_cell_name(cell_name);

        if (accumulator != node_1)
        {
            if (debug)
            {
                Rcpp::Rcout
                    << "Use the parameter at the first response to simulate\n";
                Rcpp::Rcout << "Current and previous conditions: " << condition
                            << ", " << prev_condition << std::endl;
                Rcpp::Rcout << "Current node/accumulator is: " << accumulator
                            << ", which is not the node 1, " << node_1
                            << ". skip\n"
                            << std::endl;
            }
            continue;
        }
        if (!valid || (cell_idx > 0 && condition == prev_condition))
        {
            if (debug)
            {
                Rcpp::Rcout << "This condition has done before:\n";
                Rcpp::Rcout << "Current and previous conditions:" << condition
                            << ", " << prev_condition << ". skip \n";
            }

            continue;
        }

        prev_condition = condition;

        design->set_parameter_values(cell_idx, parameters);
        auto cell_parameters = design->m_parameter_matrix[cell_idx];

        // is_lower is not critical, I use the default set in ddm.h
        ddm_obj.set_parameters(cell_parameters);

        bool is_valid = ddm_obj.validate_parameters(debug);

        if (!is_valid && debug)
        {
            Rcpp::Rcout << "DDM model (current condition = " << cell_name
                        << ") to simulate all " << design->m_n_accumulator
                        << " possible responses\n";
            ddm_obj.print("simulate_ddm_trials");
            throw std::runtime_error("Invalid paraemters");
        }

        // auto cell_data = ddm::simulate_trials(ddm_obj, n_trial_per_cell,
        // debug);
        auto cell_data = ddm::simulate_trials(ddm_obj, n_trial_per_cell, debug);

        for (const auto &trial : cell_data)
        {
            results.add_trial(trial.first, trial.second, condition);
        }
    }
}

//' @rdname simulate_ddm_trials
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame simulate_ddm_trials(
    const Rcpp::S4 &rt_model_r, const Rcpp::NumericVector &parameters_r,
    Rcpp::NumericVector time_parameters_r = Rcpp::NumericVector::create(),
    unsigned int n_trial = 1L, bool debug = false)
{
    if (time_parameters_r.size() == 0)
    {
        time_parameters_r = Rcpp::NumericVector::create(
            Rcpp::Named("t_min") = -0.5, Rcpp::Named("t_max") = 0.5,
            Rcpp::Named("dt") = 0.01);
    }

    auto time_parameters = Rcpp::as<std::vector<double>>(time_parameters_r);

    // 1. Set pu a new design instance
    auto d_ptr = new_design_light_rt_model(rt_model_r);

    // 2. Trial Distribution Setup
    unsigned int n_condition = d_ptr->m_n_cell / d_ptr->m_n_accumulator;
    unsigned int n_trial_per_cell = n_trial / n_condition;

    // 3. Simulation
    SimulationResults results(n_trial);
    static ddm::ddm_class ddm_obj;
    ddm_obj.set_times(time_parameters);

    auto parameters = Rcpp::as<std::vector<double>>(parameters_r);

    simulate_conditions(d_ptr, ddm_obj, parameters, n_trial_per_cell, results,
                        debug);

    // new_DataFrame is in lba_type_casting.
    return new_DataFrame(results);
}