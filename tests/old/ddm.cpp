#include "ddm.h"

namespace ddm
{
void ddm_class::set_precision()
{
    // Try to achieve an accuracy of approximately 10^{-p} for the CDF.
    TUNE_PDE_DT_MIN = std::pow(10, -0.400825 * m_precision - 1.422813);
    TUNE_PDE_DT_MAX = std::pow(10, -0.627224 * m_precision + 0.492689);
    TUNE_PDE_DT_SCALE = std::pow(10, -1.012677 * m_precision + 2.261668);

    TUNE_DZ = std::pow(10, -0.5 * m_precision - 0.033403); // CDF
    TUNE_DV = std::pow(10, -1.0 * m_precision + 1.4);
    TUNE_DT0 = std::pow(10, -0.5 * m_precision - 0.323859);

    TUNE_INT_T0 = 0.089045 * std::exp(-1.037580 * m_precision); // PDF
    TUNE_INT_Z = 0.508061 * std::exp(-1.022373 * m_precision);

    TUNE_SV_EPSILON = std::pow(10, -(m_precision + 2.0)); // Used by pdiffusion

    // Used by ddiffusion and pdiffusion
    TUNE_SZ_EPSILON = std::pow(10, -(m_precision + 2.0));
    TUNE_ST0_EPSILON = std::pow(10, -(m_precision + 2.0));
}

ddm_class::ddm_class(const std::vector<double> &parameters, bool is_minus)
{
    set_parameters(parameters, is_minus);
}

void ddm_class::set_parameters(const std::vector<double> &parameters,
                               bool is_minus)
{
    // Store common expressions
    double s = parameters[PARAM_s];
    double scaling = (s == 1) ? 1.0 : (1.0 / s);

    // Used in F_calculator
    m_s = s;
    m_precision = parameters[PARAM_precision];
    m_d = parameters[PARAM_d];
    m_szr = parameters[PARAM_szr];
    m_st0 = parameters[PARAM_st0];

    m_zr = is_minus ? parameters[PARAM_zr] : 1 - parameters[PARAM_zr];

    m_a = parameters[PARAM_a] * scaling;
    m_sv = parameters[PARAM_sv] * scaling;
    m_v = is_minus ? parameters[PARAM_v] * scaling
                   : (-parameters[PARAM_v]) * scaling;

    m_t_offset = parameters[PARAM_t0] + 0.5 * parameters[PARAM_st0];
    m_t0 = (m_st0 == 0) ? parameters[PARAM_t0] : m_t_offset;

    m_plus = is_minus ? false : true;
    m_sz = m_szr * m_a;

    set_precision();
}

void ddm_class::set_parameters(
    const std::vector<std::vector<double>> &parameters, bool is_minus)
{
    // Rcpp::Rcout << input[0][0] << " " << input[1][0] << " " << input[2][0];

    // Store common expressions;
    // fastdm is a relative model and has only one accumulator.
    // Thus, both columns in the parameter matrix are the same.
    double s = parameters[PARAM_s][0];
    double scaling = (s == 1) ? 1.0 : (1.0 / s);

    // Used in F_calculator
    m_s = s;
    m_precision = parameters[PARAM_precision][0];
    m_d = parameters[PARAM_d][0];
    m_szr = parameters[PARAM_szr][0];
    m_st0 = parameters[PARAM_st0][0];

    m_zr = is_minus ? parameters[PARAM_zr][0] : 1 - parameters[PARAM_zr][0];

    m_a = parameters[PARAM_a][0] * scaling;
    m_sv = parameters[PARAM_sv][0] * scaling;
    m_v = is_minus ? parameters[PARAM_v][0] * scaling
                   : (-parameters[PARAM_v][0]) * scaling;

    m_t_offset = parameters[PARAM_t0][0] + 0.5 * parameters[PARAM_st0][0];
    m_t0 = (m_st0 == 0) ? parameters[PARAM_t0][0] : m_t_offset;

    m_plus = is_minus ? false : true;
    m_sz = m_szr * m_a;

    set_precision();
}

bool ddm_class::validate_parameters(bool is_debug)
{
    bool valid = true;

    if (m_a <= 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter a = " << m_a << std::endl;
        }
    }
    if (m_szr < 0 || m_szr > 1)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter szr = " << m_szr
                        << std::endl;
        }
    }
    if (m_st0 < 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter st0 = " << m_st0
                        << std::endl;
        }
    }
    if (m_sv < 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter sv = " << m_sv
                        << std::endl;
        }
    }
    if (m_t0 - std::fabs(0.5 * m_d) - 0.5 * m_st0 < 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter combination t0 = " << m_t0
                        << ", d = " << m_d << ", st0 =" << m_st0 << std::endl;
        }
    }
    if (m_zr - 0.5 * m_szr <= 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter combination zr = " << m_zr
                        << ", szr = " << m_szr << std::endl;
        }
    }
    if (m_zr + 0.5 * m_szr >= 1)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter combination zr = " << m_zr
                        << ", szr = " << m_szr << std::endl;
        }
    }
    if (m_s <= 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid diffusion constant " << m_s;
        }
    }
    return valid;
}

void ddm_class::print(std::string header) const
{
    // Formatting constants
    const int col_width = 10;
    const std::string separator = " | ";

    std::string boundary = m_plus ? "upper" : "lower";

    // Print header
    Rcpp::Rcout << "\n=== " << header << " ===\n\n";

    // Parameter table
    Rcpp::Rcout << std::left << std::setw(col_width) << "Parameter" << separator
                << std::setw(col_width) << "Value" << "\n";

    Rcpp::Rcout << std::string(col_width + separator.length() + col_width, '-')
                << "\n";

    // Core parameters
    Rcpp::Rcout << std::left << std::setw(col_width) << "a" << separator
                << std::setw(col_width) << m_a << "\n"
                << std::setw(col_width) << "v" << separator
                << std::setw(col_width) << m_v << "\n"
                << std::setw(col_width) << "t0" << separator
                << std::setw(col_width) << m_t0 << "\n"
                << std::setw(col_width) << "d" << separator
                << std::setw(col_width) << m_d << "\n"
                << std::setw(col_width) << "zr" << separator
                << std::setw(col_width) << m_zr << "\n";

    // Variability parameters
    // Rcpp::Rcout << "\nVariability Parameters:\n";
    // Rcpp::Rcout << std::string(col_width + separator.length() + col_width,
    // '-')
    //             << "\n";

    Rcpp::Rcout << std::left << std::setw(col_width) << "szr" << separator
                << std::setw(col_width) << m_szr << "\n"
                << std::setw(col_width) << "sv" << separator
                << std::setw(col_width) << m_sv << "\n"
                << std::setw(col_width) << "st0" << separator
                << std::setw(col_width) << m_st0
                << "\n"
                // Diffusion constant
                << std::setw(col_width) << "s (dc)" << separator
                << std::setw(col_width) << m_s << "\n"
                << std::setw(col_width) << "precision" << separator
                << std::setw(col_width) << m_precision << "\n"
                << std::setw(col_width) << "boundary" << separator
                << std::setw(col_width) << boundary << "\n\n";
}

// DDM theorteical density -----------------------------
double ddm_class::compute_g_series(double ta2, double zr, bool small_time,
                                   int N)
{
    if (ta2 <= 0)
    {
        throw std::runtime_error("dt / (a * a) be greater than 0.");
    }
    double out = 0;
    if (small_time)
    {
        // A3 Formula on p1217; Appendix Mathematical Details V, R, & V (2004)
        // See also Feller (1971, p359 & p370 Problem 22);
        // Note zr = z/a. v = 0 & a = 1
        const double ta2_3 = ta2 * ta2 * ta2;
        const double norm_factor = 1.0 / std::sqrt(M_2PI * ta2_3);

        for (int i = -N / 2; i <= N / 2; i++)
        {
            const double step = 2.0 * i + zr;
            out += std::exp(-step * step / (2.0 * ta2)) * step;
        }
        return out * norm_factor;
    }
    else
    {
        // Large-time approximation (A4 formula, p1217)
        for (int i = 1; i <= N; i++)
        {
            const double step = i * M_PI;
            out += std::exp(-0.5 * step * step * ta2) * std::sin(step * zr) * i;
        }
        return out * M_PI;
    }
}

double ddm_class::compute_g_factor(double dt, double a, double zr, double v,
                                   bool no_var_case)
{
    const double a2 = a * a;
    const double v2 = v * v;

    if (no_var_case)
    {
        // Case when sv = 0 (no variability)
        const double factor = std::exp(-a * zr * v - 0.5 * v2 * dt) / a2;
        return std::isfinite(factor) ? factor : 0.0;
    }
    else
    {
        // Case with variability (sv != 0)
        const double sv2 = m_sv * m_sv; // use the member variable
        const double zr2 = zr * zr;

        const double denominator = dt * sv2 + 1.0;
        const double exponent =
            -0.5 * (v2 * dt + 2 * v * a * zr - a2 * zr2 * sv2) / denominator;
        const double factor =
            std::exp(exponent) / (a2 * std::sqrt(denominator));

        return std::isfinite(factor) ? factor : 0.0;
    }
}

std::pair<int, int> ddm_class::get_N(double dt, double ta2, double eps)
{
    int N_large = static_cast<int>(std::ceil(1.0 / (M_PI * std::sqrt(dt))));
    if (M_PI * ta2 * eps < 1.0)
    {
        double term0 = -2 * std::log(M_PI * ta2 * eps) / (M_PI * M_PI * ta2);
        double term1 = std::ceil(std::sqrt(term0));
        N_large = std::max(N_large, static_cast<int>(term1));
    }

    int N_small;
    if (2.0 * std::sqrt(M_2PI * ta2) * eps < 1.0)
    {
        double term0 = std::log(2.0 * eps * std::sqrt(M_2PI * ta2));
        double term1 = 2.0 + std::sqrt(-2.0 * ta2 * term0);
        double term2 = std::sqrt(ta2) + 1;
        N_small = static_cast<int>(std::ceil(std::max(term2, term1)));
    }
    else
    {
        N_small = 2;
    }
    return {N_small, N_large};
}

double ddm_class::g_no_var(double dt, double a, double zr, double v)
{
    if (dt <= 0)
    {
        return 0.0;
    }

    double ta2 = dt / (a * a);
    const double factor = compute_g_factor(dt, a, zr, v, true);

    if (factor == 0)
    {
        return 0.0;
    }

    const double eps = EPSILON / factor;
    auto [N_small, N_large] = get_N(dt, ta2, eps);
    const bool use_small_time = (N_small < N_large);
    const int N = use_small_time ? N_small : N_large;
    return factor * compute_g_series(ta2, zr, use_small_time, N);
}

double ddm_class::integral_v(double dt, double a, double zr, double v)
{
    if (dt <= 0)
    {
        return 0;
    }
    if (m_sv == 0) // take m_sv from the member variable
    {
        return g_no_var(dt, a, zr, v);
    }

    const double ta2 = dt / a * a;
    // Boolean to false
    const double factor = compute_g_factor(dt, a, zr, v, false);
    if (factor == 0)
    {
        return 0.0;
    }

    const double eps = EPSILON / factor;
    auto [N_small, N_large] = get_N(dt, ta2, eps);
    const bool use_small_time = (N_small < N_large);
    const int N = use_small_time ? N_small : N_large;

    return factor * compute_g_series(ta2, zr, use_small_time, N);
}

double ddm_class::integrate_v_over_zr(double dt, double zr)
{
    double lower = zr - 0.5 * m_szr;
    double upper = zr + 0.5 * m_szr;
    double width = upper - lower;
    double tmp_N = width / TUNE_INT_Z;

    int N = std::max(4, static_cast<int>(tmp_N));
    double step = width / N;

    double result = 0;
    for (double x = lower + 0.5 * step; x < upper; x += step)
    {
        result += step * integral_v(dt, m_a, x, m_v);
    }

    return result / m_szr;
}

double ddm_class::integral_z(double dt)
{
    // integral over a uniform, fixed range of zr.
    return (m_szr < TUNE_SZ_EPSILON) ? integral_v(dt, m_a, m_zr, m_v)
                                     : integrate_v_over_zr(dt, m_zr);
}

double ddm_class::integrate_z_over_t(double dt)
{
    double lower = dt - 0.5 * m_st0;
    double upper = dt + 0.5 * m_st0;

    double width = upper - lower;
    double tmp_N = width / TUNE_INT_T0;
    int N = std::max(4, static_cast<int>(tmp_N));
    double step = width / N;

    double result = 0;

    for (double x = lower + 0.5 * step; x < upper; x += step)
    {
        result += step * integral_z(x);
    }

    return result / m_st0;
}

double ddm_class::integral_t0(double dt)
{
    // integral over a uniform, fixed range of t, depending on
    // integrate_z_over_t, and integral_z_g, which depends
    // on integrate_v_over_zr and integral_v_g.
    return (m_st0 < TUNE_ST0_EPSILON) ? integral_z(dt) : integrate_z_over_t(dt);
}

double ddm_class::g(double rt)
{
    return integral_t0(rt - m_t_offset);
}

// DDM CDF solved by the PDE method (PDE Solver)----------
void f_plain::solve(double left, double mid, double right)
{
    // Ensure tmp_buffer is large enough (resize only if necessary)
    if (m_tmp_buffer.size() < static_cast<size_t>(m_n_solver - 1))
    {
        m_tmp_buffer.resize(m_n_solver - 1);
    }

    // Forward substitution (Thomas algorithm)
    double inv_pivot = 1.0 / mid;
    m_tmp_buffer[0] = right * inv_pivot;
    m_data.f[1] = m_tmp_vector[1] * inv_pivot;

    for (int i = 2; i <= m_n_solver; ++i)
    {
        inv_pivot = 1.0 / (mid - left * m_tmp_buffer[i - 2]);
        m_data.f[i] = (m_tmp_vector[i] - left * m_data.f[i - 1]) * inv_pivot;
        if (i < m_n_solver)
        { // Don't compute for last element (superdiagonal ends at n-1)
            m_tmp_buffer[i - 1] = right * inv_pivot;
        }
    }

    // Backward substitution
    for (int i = m_n_solver - 1; i >= 1; --i)
    {
        m_data.f[i] -= m_tmp_buffer[i - 1] * m_data.f[i + 1];
    }
}

void f_plain::make_step(double dt)
{
    // std::vector<double> tmp_vector(m_N + 1);
    double dz2 = m_data.dz * m_data.dz;
    double dzv = m_data.dz * m_data.v;
    double left = (1 - dzv) / (2 * dz2);
    double mid = -1 / dz2;
    double right = (1 + dzv) / (2 * dz2);

    // Prepare right-hand side
    m_tmp_vector[1] = dt * left * m_data.f[0] +
                      (1 + 0.5 * dt * mid) * m_data.f[1] +
                      0.5 * dt * right * m_data.f[2];

    for (int i = 2; i < m_N - 1; ++i)
    {
        m_tmp_vector[i] = 0.5 * dt * left * m_data.f[i - 1] +
                          (1 + 0.5 * dt * mid) * m_data.f[i] +
                          0.5 * dt * right * m_data.f[i + 1];
    }

    m_tmp_vector[m_N - 1] = 0.5 * dt * left * m_data.f[m_N - 2] +
                            (1 + 0.5 * dt * mid) * m_data.f[m_N - 1] +
                            dt * right * m_data.f[m_N];

    solve(-0.5 * dt * left, 1 - 0.5 * dt * mid, -0.5 * dt * right);
}

void f_plain::advance_to(double t1)
{
    // t1 is the decision time
    // advance_to move dt step at every iteration until exceed the t1
    // (ie the decision time entered by the user)
    while (m_data.curr_t < t1)
    {
        double dt = std::min(m_data.dt_min + m_data.dt_scale * m_data.curr_t,
                             m_data.dt_max);
        dt = std::min(dt, t1 - m_data.curr_t);
        make_step(dt);
        m_data.curr_t += dt;
    }
}

/*  f_plain ********************************************** */
double f_plain::Data::f_limit(double z)
{
    if (std::fabs(v) < 1e-8)
    {
        return 1 - z / a;
    }
    else
    {
        return (std::exp(-2 * v * z) - std::exp(-2 * v * a)) /
               (1 - std::exp(-2 * v * a));
    }
}

f_plain::~f_plain()
{
    // Rcpp::Rcout << "f_plain destructor\n";
}

f_plain::f_plain(const ddm_class &p) : m_data{}
{
    // V&V's note "N must be even, otherwise the case szr == 1 fails"
    // TODO: test the case of szr = 1, when m_N is not even.
    m_N = std::max(4, 2 * static_cast<int>(p.m_a * 0.5 / p.TUNE_DZ + 0.5));

    m_data.f.resize(m_N + 1, 0.0);
    m_data.a = p.m_a;
    m_data.d = p.m_d;
    m_data.t0 = p.m_t0;

    m_data.dt_min = p.TUNE_PDE_DT_MIN;
    m_data.dt_max = p.TUNE_PDE_DT_MAX;
    m_data.dt_scale = p.TUNE_PDE_DT_SCALE;
    m_data.dz = p.m_a / m_N;

    // Reverse back the operation in ddm_class constructor
    // to be in line with the old version.
    m_data.v = (p.m_plus) ? -p.m_v : p.m_v;
    m_n_solver = m_N - 1;

    // Rcpp::Rcout << "f_plain: dz, nF, v = [" << m_data.dz << ", " << m_N + 1
    //             << ", " << m_data.v << "]\n";
}
void f_plain::start(bool plus)
{
    m_plus = plus;
    // TODO: How could we remove the plus operation?

    m_data.t_offset = m_data.t0 - m_data.d * (m_plus ? 0.5 : -0.5);
    m_data.curr_t = 0.0;

    // The CDF over the space z at the time curr_t + t_offset
    // Rcpp::Rcout << "a, v = " << m_data.a << ", " << m_data.v << std::endl;
    // Rcpp::Rcout << "z & m_data.f[i] = \n";

    m_data.f[0] = (m_plus) ? 1.0 : 0.0;
    for (int i = 1; i < m_N; ++i)
    {
        double z = get_z(i);
        m_data.f[i] = m_data.f_limit(z);
        // Rcpp::Rcout << z << ", ";
        // Rcpp::Rcout << m_data.f[i] << ", ";
    }
    m_data.f[m_N] = (m_plus) ? 1.0 : 0.0;

    // std::vector<double>
    m_tmp_vector.resize(m_N + 1); // 35
    m_tmp_buffer.resize(m_N - 2);
    // Rcpp::Rcout << "\nf_plain_start: m_plus t_offset, curr_t = [" << m_plus
    //             << ", " << m_data.t_offset << "," << m_data.curr_t << "]\n";
}

const std::vector<double> f_plain::get_f(double t)
{
    t -= m_data.t_offset; // convert RT TO DT
    // Rcpp::Rcout << "t = " << t << ", " << m_data.curr_t << "\n";

    if (t > m_data.curr_t)
    {
        advance_to(t);
        // Finally the curr_t will be set to the decision time
        m_data.curr_t = t;
    }

    // Rcpp::Rcout << "m_N = " << m_N << ", " << m_data.f.size() << "\n";
    // N = f_ptr->m_N;
    // for (size_t i = 0; i < m_data.f.size(); ++i)
    // {
    //     Rcpp::Rcout << m_data.f[i] << ", ";
    // }
    // Rcpp::Rcout << std::endl;

    return m_data.f;
}
double f_plain::get_z(int i) const
{
    // Rcpp::Rcout << "f_plain get_z\n";
    return i * m_data.dz;
}
/*  sz ********************************************** */
f_sz::~f_sz()
{
    // Rcpp::Rcout << "\nf_sz destructor\n";
}

f_sz::f_sz(const ddm_class &p)
{
    // Rcpp::Rcout << "f_sz ...\n";
    m_data.f_plain_ptr = std::make_unique<f_plain>(p);

    int N = m_data.f_plain_ptr->m_N; // this is a separate m_N store in f_plain
    double dz = m_data.f_plain_ptr->get_z(1) - m_data.f_plain_ptr->get_z(0);
    double tmp = p.m_sz / (2.0 * dz);

    int k = static_cast<int>(std::ceil(p.m_sz / (2.0 * dz)) + 0.5);
    if (2.0 * k > N)
    {
        throw std::runtime_error("2 * k > N");
    }

    // The member variable m_N for f_sz is inherited from f_calculator
    m_N = N - 2 * k;

    m_data.k = k; // the average involves 2*k+1 cells
    m_data.q = k - tmp;
    m_data.dzoversz = dz / p.m_sz;
    m_data.avg.resize(m_N + 1, 0.0); // This must be latter

    // Rcpp::Rcout << "f_sz: [N, dz, k, q, f & navg] = [" << m_N << ", " << dz
    //             << ", " << k << ", " << m_data.q << ", " << m_data.dzoversz
    //             << ", " << m_N + 1 << "]\n";
}

const std::vector<double> f_sz::get_f(double t)
{
    // Rcpp::Rcout << "f_sz get_f\n";

    std::vector<double> F = m_data.f_plain_ptr->get_f(t);

    // Rcpp::Rcout << "F length = " << F.size() << "\n";

    // for (size_t i = 0; i < F.size(); ++i)
    // {
    //     Rcpp::Rcout << std::fixed << std::setprecision(2) << F[i];
    //     if (i < F.size() - 1)
    //     {
    //         Rcpp::Rcout << ", ";
    //     }
    // }
    // Rcpp::Rcout << "\n\n";

    int m = 2.0 * m_data.k;
    // double q = m_data.q;
    double q2 = m_data.q * m_data.q;
    double half_one_minus_q_sqrt = 0.5 * (1 - m_data.q) * (1 - m_data.q);
    double dzoversz = m_data.dzoversz;

    double tmp;
    if (m >= 3)
    {
        for (int i = 0; i <= m_N; ++i) // 0 to 30 but F has 35 elements
        {
            tmp = F[i] * half_one_minus_q_sqrt;
            tmp += F[i + 1] * (1 - 0.5 * q2);

            for (int j = i + 2; j < i + m - 1; ++j)
            {
                tmp += F[j];
            }

            tmp += F[i + m - 1] * (1 - 0.5 * q2);
            tmp += F[i + m] * half_one_minus_q_sqrt;

            m_data.avg[i] = tmp * dzoversz;
        }
    }
    else
    {
        /* m == 2 */
        for (int i = 0; i <= m_N; ++i)
        {
            tmp = F[i] * half_one_minus_q_sqrt;
            tmp += F[i + 1] * (1 - q2);
            tmp += F[i + 2] * half_one_minus_q_sqrt;
            m_data.avg[i] = tmp * dzoversz;
        }
    }
    // /* m == 1 is impossible here */

    // Rcpp::Rcout << "m_data.avg length = " << m_data.avg.size() << "\n";

    // for (size_t i = 0; i < m_data.avg.size(); ++i)
    // {
    //     Rcpp::Rcout << std::fixed << std::setprecision(2) << m_data.avg[i];
    //     if (i < m_data.avg.size() - 1)
    //     {
    //         Rcpp::Rcout << ", ";
    //     }
    // }
    // Rcpp::Rcout << "\n\n";

    // Rcpp::Rcout << "f_sz get f: N & m = " << m_N << ", " << m << "\n";
    return m_data.avg;
}
double f_sz::get_z(int i) const
{
    // Rcpp::Rcout << "f_sz get_z\n";
    return m_data.f_plain_ptr->get_z(i + m_data.k);
}
void f_sz::start(bool plus)
{
    m_data.f_plain_ptr->start(plus);
}
/*** sv **************************************************/
f_sv::~f_sv()
{
    // Rcpp::Rcout << "f_sv destructor\n";
}

f_sv::f_sv(const ddm_class &p)
{
    // Rcpp::Rcout << "f_sv ...\n";
    int nv = std::max(3, static_cast<int>(p.m_sv / p.TUNE_DV + 0.5));
    m_data.f_sz_ptrs.resize(nv);
    ddm_class p_copy = p;
    p_copy.m_sv = 0.0;

    double x;
    for (int j = 0; j < nv; ++j)
    {
        x = Rf_qnorm5((0.5 + j) / nv, 0.0, 1.0, 1, 0);
        p_copy.m_v = p.m_sv * x + p.m_v;
        // Rcpp::Rcout << "\n(sv*x + v) =" << p_copy.m_v << "\n";

        // Create nv different f_plain unique pointer
        m_data.f_sz_ptrs[j] = std::make_unique<f_sz>(p_copy);
    }
    // Rcpp::Rcout << std::endl;

    m_N = m_data.f_sz_ptrs[0]->m_N;

    m_data.avg.resize(m_N + 1, 0.0);
    m_data.nv = nv;

    // Rcpp::Rcout << "f_sv: [nv, drift rate, m_N & avg size] = [" << nv << ", "
    //             << p.m_v << ", " << m_N << ", " << m_N + 1 << "]:\n\n";
}

const std::vector<double> f_sv::get_f(double t)
{
    // Rcpp::Rcout << "f_sv get_f...\n";
    std::vector<double> F = m_data.f_sz_ptrs[0]->get_f(t);
    m_data.avg = F;
    // for (size_t i = 0; i < m_data.avg.size(); ++i)
    // {
    //     Rcpp::Rcout << std::fixed << std::setprecision(2) << m_data.avg[i];
    //     if (i < m_data.avg.size() - 1)
    //     {
    //         Rcpp::Rcout << ", ";
    //     }
    // }
    // Rcpp::Rcout << "\n\n";

    for (int j = 1; j < m_data.nv; ++j)
    {
        F = m_data.f_sz_ptrs[j]->get_f(t);
        for (int i = 0; i < m_N + 1; ++i)
        {
            m_data.avg[i] += F[i];
            // Rcpp::Rcout << std::fixed << std::setprecision(2) <<
            // m_data.avg[i]; if (i < m_data.avg.size())
            // {
            //     Rcpp::Rcout << ", ";
            // }
        }
        // Rcpp::Rcout << "\n";
    }
    // Rcpp::Rcout << "\n\n";

    // m_data.avg.size() = m_N + 1
    for (int i = 0; i < m_N + 1; ++i)
    {
        m_data.avg[i] /= m_data.nv;
        // Rcpp::Rcout << std::fixed << std::setprecision(2) << m_data.avg[i];
        // if (i < m_data.avg.size() - 1)
        // {
        //     Rcpp::Rcout << ", ";
        // }
    }
    // Rcpp::Rcout << "\n\n";

    return m_data.avg;
}

double f_sv::get_z(int i) const
{
    // Rcpp::Rcout << "f_sv get_z\n";
    return m_data.f_sz_ptrs[0]->get_z(i);
}

void f_sv::start(bool plus)
{
    for (int i = 0; i < m_data.nv; ++i)
    {
        m_data.f_sz_ptrs[i]->start(plus);
    }
}

/*** st0 **************************************************/
f_st0::~f_st0()
{
    // Rcpp::Rcout << "f_st0 destructor\n";
}

f_st0::f_st0(const ddm_class &p)
{
    // Rcpp::Rcout << "f_st0 ...\n";
    m_data.f_sv_ptr = std::make_unique<f_sv>(p);
    m_data.st0 = p.m_st0;
    m_data.M = std::max(3, static_cast<int>(p.m_st0 / p.TUNE_DT0 + 1.5));
    m_data.dt = p.m_st0 / (m_data.M - 2); // time difference
    m_N = m_data.f_sv_ptr->m_N;
    m_ncol = m_N + 1;

    // m_data.values.resize(m_data.M * (m_N + 1), 0.0); // M row x n_N + 1
    // column?
    m_data.valid.resize(m_data.M, false);
    m_data.avg.resize(m_ncol, 0.0);
    m_data.base = 0;

    if (m_data.M == 0)
    {
        throw std::runtime_error("M = 0");
    }

    // Rcpp::Rcout << "[f_st0 ...]: m_N & M = [" << m_N << ", " << m_data.M
    //             << "]. size of values & avg : " << m_data.values.size() << ",
    //             "
    //             << m_data.avg.size() << ":\n\n";
}

void f_st0::get_row(int j, double a)
{
    if (j < 0 || j > m_data.M)
    {
        std::cout << "j and M = " << j << ", " << m_data.M << std::endl;
        throw std::runtime_error("j must be within 0 ~ M (inclusive)");
    }

    const int idx = (m_data.base + j) % m_data.M;
    // Rcpp::Rcout << m_data.base << " + " << j << " % " << m_data.M << " = "
    //             << idx << ", do compute? " <<
    //             std::to_string(!m_data.valid[idx])
    //             << std::endl;

    // std::vector<double> out;
    // Rcpp::Rcout << "t (in get row): " << std::endl;

    if (!m_data.valid[idx])
    {
        m_data.value_matrix[idx] =
            m_data.f_sv_ptr->get_f(m_data.start_time + j * m_data.dt);
        m_data.valid[idx] = true;
    }

    const auto &row = m_data.value_matrix[idx];
    if (a == 1)
    {
        for (int i = 0; i < m_ncol; ++i)
        {
            m_data.avg[i] += row[i];
        }
    }
    else
    {
        for (int i = 0; i < m_ncol; ++i)
        {
            m_data.avg[i] += a * row[i];
        }
    }

    // Rcpp::Rcout << "m_row size = " << m_data.value_matrix[j].size()
    //             << std::endl;

    // for (size_t i = 0; i < m_data.value_matrix[j].size(); ++i)
    // {
    //     Rcpp::Rcout << std::fixed << std::setprecision(2)
    //                 << m_data.value_matrix[j][i];
    //     if (i < m_data.value_matrix[j].size())
    //     {
    //         Rcpp::Rcout << ", ";
    //     }
    // }
    // Rcpp::Rcout << "\n\n";
}

const std::vector<double> f_st0::get_f(double t)
{
    // Calculate bounds and validate
    const double lower_bound = t - 0.5 * m_data.st0;
    const double upper_bound = t + 0.5 * m_data.st0;

    // Calculate shift with bounds checking
    const double time_diff = lower_bound - m_data.start_time;
    const int shift = (time_diff >= m_data.M * m_data.dt)
                          ? m_data.M
                          : static_cast<int>(time_diff / m_data.dt);

    if (shift < 0)
    {
        throw std::runtime_error("shift < 0; f_st0 get_f error");
    }

    // Invalidate shifted rows
    for (int j = 0; j < shift; ++j)
    {
        m_data.valid[(m_data.base + j) % m_data.M] = false;
    }

    // Update base and start time
    if (shift < m_data.M)
    {
        m_data.start_time += shift * m_data.dt;
        m_data.base = (m_data.base + shift) % m_data.M;
    }
    else
    {
        m_data.start_time = lower_bound;
    }

    // Calculate m, q, r in a more compact way
    const double tmp = (upper_bound - m_data.start_time) / m_data.dt;
    // const int m = std::min(m_data.M - 1, static_cast<int>(std::round(tmp)));
    const int m =
        std::min(m_data.M - 1, static_cast<int>(std::ceil(tmp) + 0.5));

    const double q = (lower_bound - m_data.start_time) / m_data.dt;
    const double r = m - tmp;

    // Reset and compute averages
    std::fill(m_data.avg.begin(), m_data.avg.end(), 0.0);

    if (m >= 3)
    {
        get_row(0, 0.5 * (1 - q) * (1 - q));
        get_row(1, 1 - 0.5 * q * q);
        for (int j = 2; j < m - 1; ++j)
        {
            get_row(j, 1);
        }
        get_row(m - 1, 1 - 0.5 * r * r);
        get_row(m, 0.5 * (1 - r) * (1 - r));
    }
    else if (m == 2)
    {
        get_row(0, 0.5 * (1 - q) * (1 - q));
        get_row(1, 1 - 0.5 * (q * q + r * r));
        get_row(2, 0.5 * (1 - r) * (1 - r));
    }
    else if (m == 1)
    {
        get_row(0, 0.5 * ((1 - q) * (1 - q) - r * r));
        get_row(1, 0.5 * ((1 - r) * (1 - r) - q * q));
    }
    else
    {
        throw std::runtime_error("Invalid m value");
    }

    // Normalize results
    const double norm_factor = m_data.dt / (upper_bound - lower_bound);
    for (auto &val : m_data.avg)
    {
        val *= norm_factor;
    }

    return m_data.avg;
}

double f_st0::get_z(int i) const
{
    // Rcpp::Rcout << "f_st0 get_z\n";
    return m_data.f_sv_ptr->get_z(i);
}
void f_st0::start(bool plus)
{
    m_data.f_sv_ptr->start(plus);

    m_data.start_time = std::numeric_limits<double>::lowest();
    std::fill(m_data.valid.begin(), m_data.valid.end(), false);
    m_data.value_matrix.resize(m_data.M, std::vector<double>(m_ncol));
}

// DDM simulation inverse method
std::unique_ptr<f_calculator> select_f_calculator(const ddm_class &p,
                                                  bool debug)
{
    std::unique_ptr<f_calculator> out;
    out.reset(); // Clear previous calculator

    bool is_small_st0 = p.m_st0 < p.TUNE_DT0 * 1e-6;
    bool is_small_sz = p.m_sz < p.TUNE_SZ_EPSILON;
    bool is_small_sv = p.m_sv < p.TUNE_SV_EPSILON;

    // Case 1: st0 ≈ 0 → Use f_sv (no non-decision time variability)
    if (is_small_st0 && is_small_sv && is_small_sz)
    {
        if (debug)
        {
            Rcpp::Rcout << "st0 = " << p.m_st0 << " < " << p.TUNE_DT0 * 1e-6
                        << ". ";
            Rcpp::Rcout << "sv = " << p.m_sv << " < " << p.TUNE_SV_EPSILON
                        << ". ";
            Rcpp::Rcout << "sz = " << p.m_sz << " < " << p.TUNE_SZ_EPSILON
                        << ". ";
            Rcpp::Rcout << "Selecting f_plain.\n\n";
        }
        out = std::make_unique<f_plain>(p);
    }
    // Case 2: sv ≈ 0 → Use f_sz (no drift rate variability)
    else if (is_small_st0 && is_small_sv)
    {
        if (debug)
        {
            Rcpp::Rcout << "st0 = " << p.m_st0 << " < " << p.TUNE_DT0 * 1e-6
                        << ". ";
            Rcpp::Rcout << "sv = " << p.m_sv << " < " << p.TUNE_SV_EPSILON
                        << ". ";
            Rcpp::Rcout << "sz = " << p.m_sz << ". ";
            Rcpp::Rcout << "Selecting f_sz.\n\n";
        }
        out = std::make_unique<f_sz>(p);
    }
    // Case 3: sz ≈ 0 → Use f_plain (no start-point variability)
    else if (is_small_st0)
    {
        if (debug)
        {
            Rcpp::Rcout << "st0 = " << p.m_st0 << " < " << p.TUNE_DT0 * 1e-6
                        << ". ";
            Rcpp::Rcout << "sv = " << p.m_sv << ". ";
            Rcpp::Rcout << "sz = " << p.m_sz << ". ";
            Rcpp::Rcout << "Selecting f_sv.\n\n";
        }
        out = std::make_unique<f_sv>(p);
    }
    // Default: Full model (f_st0)
    else
    {
        if (debug)
        {
            Rcpp::Rcout << "Full model (st0=" << p.m_st0 << ", sv=" << p.m_sv
                        << ", sz=" << p.m_sz << "). Selecting f_st0.\n\n";
        }
        out = std::make_unique<f_st0>(p);
    }
    return out;
}

double get_density(std::unique_ptr<f_calculator> &f_ptr, double t, double z)
{
    // get cdf
    const std::vector<double> F = f_ptr->get_f(t);

    double out;
    if (f_ptr->m_N == 0)
    {
        out = F[0];
    }
    else
    {
        double z0 = f_ptr->get_z(0);
        double z1 = f_ptr->get_z(f_ptr->m_N);
        int i = static_cast<int>(f_ptr->m_N * (z - z0) / (z1 - z0));

        if (i < f_ptr->m_N)
        {
            double z0 = f_ptr->get_z(i);
            double z1 = f_ptr->get_z(i + 1);
            double p = (z1 - z) / (z1 - z0);
            out = p * F[i] + (1 - p) * F[i + 1];
        }
        else
        {
            out = F[f_ptr->m_N];
        }
    }
    return out;
}

// Binary search to find slot - converted to iterative version for safety
int find_slot(double target, const std::vector<double> &values, int l, int r)
{
    while (true)
    {
        int m = l + (r - l) / 2; // Prevents potential overflow
        if (m == l)
        {
            return l;
        }
        else if (values[m] > target)
        {
            r = m;
        }
        else
        {
            l = m;
        }
    }
}

// Helper function to generate uniform random values
std::vector<double> generate_uniform_samples(unsigned int n, double min,
                                             double max)
{
    std::vector<double> samples(n);
    std::generate(samples.begin(), samples.end(),
                  [=]() { return Rf_runif(min, max); });
    return samples;
}
// Helper function to find time range boundaries
void find_time_range(std::unique_ptr<f_calculator> &f_ptr, double z,
                     double &t_min, double &t_max, double F_min, double F_max)
{
    // Upper boundary search
    f_ptr->start(true);
    while (get_density(f_ptr, t_max, z) < F_max)
    {
        t_max += 0.1;
    }

    // Lower boundary search
    f_ptr->start(false);
    while (get_density(f_ptr, -t_min, z) > F_min)
    {
        t_min -= 0.1;
    }
}

// Helper function to build F table
std::vector<double> build_F_table(std::unique_ptr<f_calculator> &f_ptr,
                                  double z, double t_min, double t_max, int N)
{
    std::vector<double> F(N + 1);
    const double dt = (t_max - t_min) / N;

    // Process upper boundary
    f_ptr->start(true);
    for (int i = 0; i <= N; ++i)
    {
        const double t = t_min + i * dt;
        if (t >= 0)
            F[i] = get_density(f_ptr, t, z);
    }

    // Process lower boundary
    f_ptr->start(false);
    for (int i = N; i >= 0; --i)
    {
        const double t = -(t_min + i * dt);
        if (t >= 0)
            F[i] = get_density(f_ptr, t, z);
    }

    // Ensure values are in [0,1] range
    std::transform(F.begin(), F.end(), F.begin(),
                   [](double val) { return std::clamp(val, 0.0, 1.0); });

    std::sort(F.begin(), F.end());
    return F;
}

// Helper function to calculate RTs and boundaries
void calculate_results(const std::vector<double> &Fs,
                       const std::vector<double> &F, double t_min, double dt,
                       std::vector<std::pair<unsigned int, double>> &out)
{
    // Fs are random uniform number 0 ~ 1;
    // F is the CDF
    unsigned int n = out.size();
    for (unsigned int i = 0; i < n; ++i)
    {
        const double y = Fs[i];
        const int N = F.size() - 1;
        const int k = find_slot(y, F, 0, N - 1); // N is F.size()-1

        if (F[k] > y || y > F[k + 1])
        {
            Rcpp::Rcout << "F[k], F[k + 1], & y " << F[k] << ", " << F[k + 1]
                        << ", " << y << std::endl;
            throw std::runtime_error("y not in the valid interpolation range");
        }

        const double t = t_min + (k + (y - F[k]) / (F[k + 1] - F[k])) * dt;
        unsigned int out_bound = t >= 0;
        double out_rt = std::fabs(t);

        out[i] = {out_bound, out_rt};
    }
}

std::vector<std::pair<unsigned int, double>>
simulate_trials(const ddm_class &ddm_obj, unsigned int n, bool debug)
{
    auto Fs = generate_uniform_samples(n);
    const auto [Fs_min, Fs_max] = std::minmax_element(Fs.begin(), Fs.end());

    double t_min = -0.5, t_max = 0.5;
    auto f_ptr = select_f_calculator(ddm_obj, debug);
    const double z = ddm_obj.m_zr * ddm_obj.m_a;

    find_time_range(f_ptr, z, t_min, t_max, *Fs_min, *Fs_max);

    const int N = static_cast<int>((t_max - t_min) / 0.001 + 0.5);
    auto F = build_F_table(f_ptr, z, t_min, t_max, N);

    F.front() = std::min(F.front(), *Fs_min);
    F.back() = std::max(F.back(), *Fs_max);

    const double dt = (t_max - t_min) / N;

    std::vector<std::pair<unsigned int, double>> results(n);
    calculate_results(Fs, F, t_min, dt, results);

    return results;
}

} // namespace ddm