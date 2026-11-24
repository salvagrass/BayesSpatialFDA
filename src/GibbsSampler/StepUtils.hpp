#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>


inline Eigen::VectorXd fast_sampling_GMRF(
    Eigen::Ref<const Eigen::VectorXd> b,
    const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &Q_factor,
    std::mt19937_64 &rng){
    Eigen::VectorXd sample = Q_factor.permutationP() * b;
    sample = Q_factor.matrixL().solve(sample);
    std::normal_distribution<double> std_normal(0.0, 1.0);
    for (int i = 0; i < sample.size(); ++i)
    {
        sample(i) += std_normal(rng);
    }
    Eigen::VectorXd buffer = Q_factor.matrixU().solve(sample);
    return Q_factor.permutationPinv() * buffer;
}; //! Multiple allocations and return by value (Change it!)
//TODO See if inlining is the best practice 

// Proposal adaptation for the theta block (eta-space, d = 4).
// Supplies the proposal Cholesky L (read by the factorization recipe to form
// eta_prop = eta_curr + L * z) and adapts it from the realized chain.
class ProposalAdapter{
public:
    ProposalAdapter() = default;
    virtual ~ProposalAdapter() = default;

    // Lower-triangular Cholesky L of the proposal covariance.
    virtual const Eigen::Matrix4d &proposal_chol() const = 0;
    // Adapt after one MH iteration, given the chain's current position (eta)
    // and the acceptance probability just computed.
    virtual void update(const Eigen::Vector4d &eta, double accept_prob) = 0;
};

// Non-adaptive: constant proposal covariance.
class FixedProposal : public ProposalAdapter{
public:
    explicit FixedProposal(const Eigen::Matrix4d &L) : L_(L) {}
    const Eigen::Matrix4d &proposal_chol() const override { return L_; }
    void update(const Eigen::Vector4d&, double) override {}


private:
    Eigen::Matrix4d L_;
};



class ThetaPrior{
public:
    ThetaPrior() = default;
    virtual ~ThetaPrior() = default;

    virtual double log_pdf(const Eigen::Vector4d &theta) const = 0;
};

class PCPrior : public ThetaPrior{
public:
    PCPrior(double lk, double lv, double ls) : ThetaPrior(), lambda_kappa(lk), lambda_v(lv), lambda_sigma(ls) {}
    double log_pdf(const Eigen::Vector4d &theta) const override;

private:
    double lambda_kappa, lambda_v, lambda_sigma;
    static double f(double v_norm);
    static double f_prime(double v_norm);
    static double f_prime_over_v(double v_norm);
};




inline double PCPrior::f(double v_norm){
    return std::sqrt(EIGEN_PI * (3 * std::cosh(2 * v_norm) + 1) / 3);
};

inline double PCPrior::f_prime(double v_norm){
    return std::sqrt(EIGEN_PI) * std::sinh(2 * v_norm) + 1 / std::sqrt(std::cosh(2 * v_norm) + 1 / 3);
};

inline double PCPrior::f_prime_over_v(double v_norm){
    // Taylor expansion of third order of f'(v)/v
    // centered in 0 is f''(0)
    // f''(0) = lim_{v->0} f'(v)/v = sqrt(3*pi)
    static const double f2_0 = std::sqrt(3.0 * EIGEN_PI);
    if (v_norm < 1e-8)
        return f2_0; //! When returning in an if usually you leave out the else to not ident too much
    return f_prime(v_norm) / v_norm;
};

inline double PCPrior::log_pdf(const Eigen::Vector4d &theta) const{
    double v_norm = std::sqrt(theta(2) * theta(2) + theta(3) * theta(3));
    double res = log(lambda_kappa) + log(lambda_v) + log(f(v_norm)) + log(f_prime_over_v(v_norm)) - log(2 * EIGEN_PI) - lambda_v * (f(v_norm) - f(0)) - lambda_kappa * f(v_norm) * theta(0) + log(lambda_sigma) - theta(1) * lambda_sigma;
    return res;
};


