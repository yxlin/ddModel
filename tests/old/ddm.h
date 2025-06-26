#pragma once
#include <RcppArmadillo.h>

namespace ddm
{
enum ParamIndex
{
    PARAM_a = 0,
    PARAM_d,
    PARAM_precision,
    PARAM_s,
    PARAM_st0,
    PARAM_sv,
    PARAM_szr,
    PARAM_t0,
    PARAM_v,
    PARAM_zr,
    NUM_PARAMS // Optional: total number of parameters
};
// Core parameters: a      d       precision       s       st0     sv      sz t0
// v      z
const double EPSILON = 1e-6; // used in ddiffusion

class ddm_class
{
  public:
    double m_a;  // Boundary separation
    double m_v;  // Mean of the drift
    double m_t0; // Non-decision time
    double m_d;  // Difference between boundaries of non-decision time

    // the user enters the reltive szr
    double m_szr; // width of zr distribution
    double m_sv;  // standard deviation of v distribution
    double m_st0; // width of t0 distribution
    double m_zr;  // Mean of diffusion starting point relative to boundaries
    double m_s;   // standard deviation; sqrt(diffusion constant)
    double m_precision; // precision
    // m_plus = true uses g_plus (upper); other use g_minus (lower)
    double m_plus;
    double m_sz; // absolute sz
    // m_sz = m_szr * m_a;
    double m_t_offset;

    double TUNE_DZ;
    double TUNE_DV;
    double TUNE_DT0;

    // If std=c++11 we can use C++ defaults to set as = 1e-6;
    double TUNE_PDE_DT_MIN;
    double TUNE_PDE_DT_MAX;   // ... we can default to = 1e-6;
    double TUNE_PDE_DT_SCALE; // ... we can default to = 0.0;
    double TUNE_INT_T0;
    double TUNE_INT_Z;
    double TUNE_SV_EPSILON;
    double TUNE_SZ_EPSILON;
    double TUNE_ST0_EPSILON;

    ddm_class(const std::vector<double> &parameters = {1.0, 1.5, 0.5, 0.0, 0.0,
                                                       0.0, 0.1, 0.0, 1.0, 2.5},
              bool is_minus = true);

    void set_parameters(const std::vector<double> &parameters,
                        bool is_minus = true);

    void set_parameters(const std::vector<std::vector<double>> &parameters,
                        bool is_minus = true);

    bool validate_parameters(bool is_debug = false);
    void print(std::string header = "") const;

    // Density function
    double compute_g_series(double ta2, double zr, bool small_time, int N);
    double compute_g_factor(double dt, double a, double zr, double v,
                            bool no_var_case);

    std::pair<int, int> get_N(double dt, double ta2, double eps);

    // double g_no_var(double dt, double zr);
    double g_no_var(double dt, double a, double zr, double v);
    double integral_v(double dt, double a, double zr, double v);

    double integrate_v_over_zr(double dt, double zr);
    double integral_z(double dt);

    double integrate_z_over_t(double dt);
    double integral_t0(double dt);
    double g(double rt);

  private:
    void set_precision();
};
class f_calculator
{
  public:
    virtual ~f_calculator() = default;
    virtual const std::vector<double> get_f(double t) = 0;
    virtual double get_z(int i) const = 0;
    virtual void start(bool plus) = 0;

    int m_N;
    bool m_plus;
};
class f_plain : public f_calculator
{
  public:
    struct Data
    {
        double a = 0.0; // Initialize with default values
        double v = 0.0;
        double d = 0.0;
        double t0 = 0.0;
        double dt_min = 0.01;
        double dt_max = 0.1;
        double dt_scale = 1.5;

        double dz = 0.1;       /* z step-size */
        double t_offset = 0.0; /* time adjustment, resulting from t0 and d */
        double curr_t = 0.0; /* adjusted time, corresponding to the vector F */
        std::vector<double> f; // state at time t + t_offset; ie CDFs
        double f_limit(double z);

        Data() = default;
    };

    // f_plain();
    f_plain(const ddm_class &p);
    ~f_plain();

    // void start(int plus) override;
    const std::vector<double> get_f(double t) override;
    double get_z(int i) const override;
    void start(bool plus) override;
    // void print(std::string header = "") const;
    void advance_to(double t1);
    void make_step(double dt);

    void solve(double left, double mid, double right);

  public:
    Data m_data;
    std::vector<double> m_tmp_buffer;
    std::vector<double> m_tmp_vector;
    int m_n_solver;
};
class f_sz : public f_calculator
{
  public:
    struct Data
    {
        // Now can be default-constructed as nullptr
        std::unique_ptr<f_plain> f_plain_ptr;
        std::vector<double> avg; // the computed averages
        int k;                   // the average involves 2*k+1 cells
        double q;                // unused part of the outermost cells
        double dzoversz;         // scale factor for the integration

        // Default constructor
        Data() = default;
    };

    f_sz(const ddm_class &p);
    ~f_sz();

    const std::vector<double> get_f(double t) override;
    double get_z(int i) const override;
    void start(bool plus) override;

  public:
    Data m_data;
};
class f_sv : public f_calculator
{
  public:
    struct Data
    {
        int nv = 0; // number of points in integration
        std::vector<std::unique_ptr<f_sz>> f_sz_ptrs;
        std::vector<double> avg;

        // Default constructor
        Data() = default;
    };
    f_sv(const ddm_class &p);
    ~f_sv();

    const std::vector<double> get_f(double t) override;
    double get_z(int i) const override;
    void start(bool plus) override;

  public:
    Data m_data;
};
class f_st0 : public f_calculator
{
  public:
    struct Data
    {
        std::unique_ptr<f_sv> f_sv_ptr;
        double st0;        // variability of t0
        int M;             // number of stored grid lines
        double start_time; // t-value of first stored grid line
        double dt;         // t-spacing of stored grid lines

        // stored grid lines (length M*(N+1))
        std::vector<std::vector<double>> value_matrix;
        std::vector<bool> valid; // which lines in 'values' are valid

        // first grid line starts at pos. base**(N + 1)
        int base;
        std::vector<double> avg; // the computed average (size N+1)

        // Default constructor
        Data() = default;
    };
    f_st0(const ddm_class &p);
    ~f_st0();

    const std::vector<double> get_f(double t) override;
    double get_z(int i) const override;
    void start(bool plus) override;

    void get_row(int j, double a);

  public:
    Data m_data;
    int m_ncol;
};

std::unique_ptr<f_calculator> select_f_calculator(const ddm_class &p,
                                                  bool debug = false);
double get_density(std::unique_ptr<f_calculator> &f_ptr, double t, double z);
int find_slot(double target, const std::vector<double> &values, int l, int r);

std::vector<double> generate_uniform_samples(unsigned int n, double min = 0.0,
                                             double max = 1.0);
void find_time_range(std::unique_ptr<f_calculator> &f_ptr, double z,
                     double &t_min, double &t_max, double F_min, double F_max);
std::vector<double> build_F_table(std::unique_ptr<f_calculator> &f_ptr,
                                  double z, double t_min, double t_max, int N);

// std::pair<std::vector<double>, std::vector<unsigned int>>
void calculate_results(const std::vector<double> &Fs,
                       const std::vector<double> &F, double t_min, double dt,
                       std::vector<std::pair<unsigned int, double>> &out);

std::vector<std::pair<unsigned int, double>>
simulate_trials(const ddm_class &ddm_obj, unsigned int n, bool debug = false);
} // namespace ddm