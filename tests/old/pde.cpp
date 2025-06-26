#include "pde.h"

/*  PDE Solver Implementation  */
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

f_plain::f_plain(const ddm::ddm_class &p) : m_data{}
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

f_sz::f_sz(const ddm::ddm_class &p)
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

f_sv::f_sv(const ddm::ddm_class &p)
{
    // Rcpp::Rcout << "f_sv ...\n";
    int nv = std::max(3, static_cast<int>(p.m_sv / p.TUNE_DV + 0.5));
    m_data.f_sz_ptrs.resize(nv);
    ddm::ddm_class p_copy = p;
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

f_st0::f_st0(const ddm::ddm_class &p)
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