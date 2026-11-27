#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "SamplingUtils.hpp"

inline void fast_sampling_GMRF( 
    const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& Q_factor,
    Eigen::VectorXd& b,
    Eigen::VectorXd& buffer,
    std::mt19937_64 &rng
){
    buffer.noalias() = Q_factor.permutationP() * b;
    Q_factor.matrixL().solveInPlace(buffer);
    std::normal_distribution<double> std_normal(0.0, 1.0);
    for (int i = 0; i < buffer.size(); ++i)
    {
        buffer(i) += std_normal(rng);
    }
    Q_factor.matrixU().solveInPlace(buffer);
    b.noalias() = Q_factor.permutationPinv() * buffer;
}


class ProposalStrategy{
public:
    ProposalStrategy() = default;
    virtual ~ProposalStrategy() = default;

    // propose the new eta vector to explore
    virtual Eigen::Vector4d propose(
        const Eigen::Vector4d& eta_curr,
        std::mt19937_64& rng
    ) = 0;

    // adapt the proposal based on the last acceptance probability 
    virtual void adapt(
        const Eigen::Vector4d& eta_curr, 
        double accept_prob
    ){return;}
};

// Non-adaptive: constant proposal covariance.
class FixedProposal : public ProposalStrategy{
public:
    FixedProposal(Eigen::Matrix4d cov_matrix) : solver(cov_matrix){
        if(solver.info() != Eigen::Success)
            throw(std::runtime_error("Can't factorize covariance matrix"));
    }
    Eigen::Vector4d propose(
        const Eigen::Vector4d& eta_curr, 
        std::mt19937_64& rng
    ) override {
        Eigen::Vector4d noise;
        std::normal_distribution<double> normal(0.0,1.0);
        for(int i = 0; i<eta_curr.size();i++)
            noise[i] = normal(rng);
        return eta_curr + solver.matrixL()*noise;
    }

private:
    const Eigen::LLT<Eigen::Matrix4d> solver;
};


class ProposalHaario : public ProposalStrategy {
public:
    ProposalHaario(size_t num_saved_, size_t num_skipped_,Eigen::Matrix4d fixed_cov)
        : fix_prop(fixed_cov),
          end_burnin(num_saved_),
          num_skipped(num_skipped_),
          online_mean(Eigen::Vector4d::Zero()),
          online_cov(Eigen::Matrix4d::Zero()),
          log_scale(std::log(2.38*2.38/4)) {}

    Eigen::Vector4d propose(
        const Eigen::Vector4d& eta_curr,
        std::mt19937_64& rng
    ) override{
        if(iter <= num_skipped)
            return fix_prop.propose(eta_curr,rng);

        // Adding 5% adaptive fixed mixtures
        std::uniform_real_distribution unif(0.0,1.0);
        if(unif(rng)<0.05 || iter-num_skipped < 10 || end_burnin < num_skipped+10)
            return fix_prop.propose(eta_curr,rng);
    
        Eigen::Vector4d noise;
        std::normal_distribution<double> norm_dist(0.0,1.0);
        for(int i=0;i<noise.size();i++)
            noise[i] = norm_dist(rng);
        Eigen::Vector4d eta_prop = eta_curr + solver.matrixL()*noise;
        return eta_prop;
    }
    
    
    void adapt(const Eigen::Vector4d& eta_curr,double acc_prob) override{
        iter++;
        constexpr double epsilon = 1e-6;
        if(iter <= num_skipped)
            return;
        if(iter <= end_burnin){
            update_stats(eta_curr,iter-num_skipped);

            const double gamma_n = std::pow(static_cast<double>(iter-num_skipped),-0.6);
            log_scale += gamma_n*(acc_prob - 0.234);
            log_scale = std::clamp(log_scale,-12.0,4.0);

            const double s_d = std::exp(log_scale);
            prop_cov = s_d*online_cov;
            prop_cov.diagonal().array() += s_d*epsilon;
            solver.compute(prop_cov);
        }
    }

private:
    void update_stats(const Eigen::Vector4d& eta_curr, size_t n) {
        if (n == 1) {
            online_mean = eta_curr;
            return;
        }

        double denom = static_cast<double>(n - 1);
        Eigen::Vector4d diff_prev = eta_curr - online_mean;
        online_mean += diff_prev / static_cast<double>(n);
        Eigen::Vector4d diff_curr = eta_curr - online_mean;

        online_cov = ((n - 2.0) / denom) * online_cov +
                     (1.0 / denom) * (diff_prev * diff_curr.transpose());
    }
    
    FixedProposal fix_prop;
    
    size_t end_burnin;
    size_t num_skipped;
    size_t iter = 0;

    Eigen::Vector4d online_mean;
    Eigen::Matrix4d online_cov;
    Eigen::Matrix4d prop_cov;
    Eigen::LLT<Eigen::Matrix4d> solver;

    double log_scale; // Used for ASWAM adaptation

};


class ThetaPrior{
public:
    ThetaPrior() = default;
    virtual ~ThetaPrior() = default;

    virtual double log_pdf(const Eigen::Vector4d& theta) const = 0;
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
    return std::sqrt(EIGEN_PI * (3.0 * std::cosh(2.0 * v_norm) + 1.0) / 3.0);
};

inline double PCPrior::f_prime(double v_norm){
    return std::sqrt(EIGEN_PI) * std::sinh(2 * v_norm) / std::sqrt(std::cosh(2 * v_norm) + 1.0 / 3.0);
};

inline double PCPrior::f_prime_over_v(double v_norm){
    // Taylor expansion of third order of f'(v)/v
    // centered in 0 is f''(0)
    // f''(0) = lim_{v->0} f'(v)/v = sqrt(3*pi)
    static const double f2_0 = std::sqrt(3.0 * EIGEN_PI);
    if (v_norm < 1e-8)
        return f2_0;
    return f_prime(v_norm) / v_norm;
};

inline double PCPrior::log_pdf(const Eigen::Vector4d &theta) const{
    double v_norm = std::sqrt(theta(2) * theta(2) + theta(3) * theta(3));
    double res = log(lambda_kappa) + log(lambda_v) + log(f(v_norm)) + log(f_prime_over_v(v_norm)) - log(2 * EIGEN_PI) - lambda_v * (f(v_norm) - f(0)) - lambda_kappa * f(v_norm) * theta(0) + log(lambda_sigma) - theta(1) * lambda_sigma;
    return res;
};


class EGPrior : public ThetaPrior{
public:
    EGPrior(double lambda_kappa,double lambda_sigma,double prec_v) 
    : lambda_k(lambda_kappa),lambda_sigma(lambda_sigma),prec_v(prec_v){}

    double log_pdf(const Eigen::Vector4d& theta) const override{
        double v1 = theta[ThetaIdx::V1];
        double v2 = theta[ThetaIdx::V2];
        return -lambda_k*theta[ThetaIdx::KAPPA] 
        -  lambda_sigma*theta[ThetaIdx::SIGMA_BETA]
        - 0.5*prec_v*(v1*v1 + v2*v2);
    }
private:
    const double lambda_k, lambda_sigma,prec_v;
};

