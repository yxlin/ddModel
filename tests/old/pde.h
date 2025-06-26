#pragma once
#include "ddm.h"
#include <RcppArmadillo.h>

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
    f_plain(const ddm::ddm_class &p);
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

    f_sz(const ddm::ddm_class &p);
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
    f_sv(const ddm::ddm_class &p);
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
    f_st0(const ddm::ddm_class &p);
    ~f_st0();

    const std::vector<double> get_f(double t) override;
    double get_z(int i) const override;
    void start(bool plus) override;

    void get_row(int j, double a);

  public:
    Data m_data;
    int m_ncol;
};

// void ddm_class::r(unsigned int n,
//                   std::vector<std::pair<unsigned int, double>> &out)
// {
//     // Initialize min/max values
//     double Fs_min = 1.0;
//     double Fs_max = 0.0;
//     double t_min = -0.5;
//     double t_max = 0.5;

//     // Generate F values using vector instead of raw array
//     std::vector<double> Fs(n);
//     for (auto &val : Fs)
//     {
//         val = Rf_runif(0.0, 1.0);
//         Fs_min = std::min(Fs_min, val);
//         Fs_max = std::min(Fs_max, val);
//     }
// }
