#ifndef SAMPLER_CPP
#define SAMPLER_CPP

#define EIGEN_PATH "C:/Users/rizhi/Desktop/experiments/cpp/eigen-3.4.0/Eigen/Dense"
//#include <Eigen/Dense>
#include EIGEN_PATH
#include <random>
#include <vector>
#include <optional>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <fstream>

#include <iostream>
#include <chrono>
#include <algorithm>
#include <emscripten.h>

// Tolerances matching cobra defaults (model.tolerance)
const double FEASIBILITY_TOL = 1e-7;
const double BOUNDS_TOL = 1e-7;
const int MAX_TRIES = 100;

struct SamplerProblem {
    // Equalities: Aeq x = beq (stoichiometric matrix in variable space)
    Eigen::MatrixXd equalities;   // (m x 2n)
    Eigen::VectorXd beq;          // (m)

    // Variable bounds in variable space (all non-negative)
    Eigen::VectorXd variable_lb;  // (2n)
    Eigen::VectorXd variable_ub;  // (2n)
    std::vector<bool> variable_fixed; // (2n)

    bool homogeneous;
};

class OptGPSampler {
public:
    OptGPSampler(const SamplerProblem& prob, const Eigen::MatrixXd& warmup,
                 int thinning, unsigned int seed)
        : prob_(prob), thinning_(thinning), rng_(seed),
          uni01_(0.0, 1.0), n_vars_(prob.variable_lb.size()),
          retries_(0), n_samples_(0)
    {
        // Compute nproj like Python: min(n_vars^3, 1e6)
        long long nv3 = (long long)n_vars_ * n_vars_ * n_vars_;
        nproj_ = (int)std::min(nv3, (long long)1000000);

        // Remove redundant warmup points (matches Python hr_sampler._is_redundant)
        warmup_ = removeRedundant(warmup);
        n_warmup_ = warmup_.rows();

        if (n_warmup_ < 2)
            throw std::runtime_error("Warmup has fewer than 2 non-redundant points");

        // Precompute decomposition for reprojection
        decomp_ = (prob_.equalities * prob_.equalities.transpose()).fullPivLu();

        // Initialize center as mean of warmup (matches Python OptGPSampler.__init__)
        center_ = warmup_.colwise().mean();
    }

    // Matches cobra/sampling/core.py step() exactly
    Eigen::VectorXd step(const Eigen::VectorXd& x, const Eigen::VectorXd& delta,
                         double fraction = -1.0, int tries = 0)
    {
        // Python: valid = (abs(delta) > feasibility_tol) & ~variable_fixed
        // Python: valphas = ((1.0 - bounds_tol) * variable_bounds - x)[:, valid] / delta[valid]
        //
        // variable_bounds is (2, n_vars): row 0 = lower, row 1 = upper
        // So alphas come from both lower and upper bound constraints per variable
        std::vector<double> alphas;

        for (int i = 0; i < n_vars_; i++) {
            if (std::abs(delta(i)) <= FEASIBILITY_TOL || prob_.variable_fixed[i])
                continue;

            // Shrink bounds slightly inward: (1.0 - bounds_tol) * bounds
            double shrunk_lb = (1.0 - BOUNDS_TOL) * prob_.variable_lb(i);
            double shrunk_ub = (1.0 - BOUNDS_TOL) * prob_.variable_ub(i);

            alphas.push_back((shrunk_lb - x(i)) / delta(i));
            alphas.push_back((shrunk_ub - x(i)) / delta(i));
        }

        // No inequality constraints in metabolic models (prob.bounds.shape[0] == 0)

        // Python: pos_alphas = alphas[alphas > 0.0]
        //         neg_alphas = alphas[alphas <= 0.0]
        //         alpha_range = [neg_alphas.max() if len(neg_alphas) > 0 else 0,
        //                        pos_alphas.min() if len(pos_alphas) > 0 else 0]
        double max_alpha = 0.0;
        double min_alpha = 0.0;
        bool has_pos = false, has_neg = false;

        for (double a : alphas) {
            if (a > 0.0) {
                max_alpha = has_pos ? std::min(max_alpha, a) : a;
                has_pos = true;
            } else {
                min_alpha = has_neg ? std::max(min_alpha, a) : a;
                has_neg = true;
            }
        }

        if (!has_pos) max_alpha = 0.0;
        if (!has_neg) min_alpha = 0.0;

        // Python: if fraction: alpha = alpha_range[0] + fraction * (...)
        //         else: alpha = np.random.uniform(alpha_range[0], alpha_range[1])
        double alpha;
        if (fraction >= 0.0)
            alpha = min_alpha + fraction * (max_alpha - min_alpha);
        else
            alpha = min_alpha + uni01_(rng_) * (max_alpha - min_alpha);

        Eigen::VectorXd p = x + alpha * delta;

        // Check bounds violation (matches Python _bounds_dist)
        double lb_dist = std::numeric_limits<double>::infinity();
        double ub_dist = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n_vars_; i++) {
            lb_dist = std::min(lb_dist, p(i) - prob_.variable_lb(i));
            ub_dist = std::min(ub_dist, prob_.variable_ub(i) - p(i));
        }
        bool bounds_violated = (lb_dist < -BOUNDS_TOL) || (ub_dist < -BOUNDS_TOL);

        // Stuck detection: abs(abs(alpha_range).max() * delta).max() < bounds_tol
        double alpha_abs_max = std::max(std::abs(min_alpha), std::abs(max_alpha));
        double max_step_contribution = 0.0;
        for (int i = 0; i < n_vars_; i++)
            max_step_contribution = std::max(max_step_contribution, std::abs(alpha_abs_max * delta(i)));
        bool stuck = max_step_contribution < BOUNDS_TOL;

        if (bounds_violated || stuck) {
            if (tries > MAX_TRIES) {
                throw std::runtime_error(
                    "Cannot escape sampling region, model seems numerically unstable.");
            }
            retries_++;
            // Python: reset to CENTER, step toward random warmup point
            // return step(sampler, sampler.center, newdir - sampler.center, None, tries + 1)
            int pi = static_cast<int>(uni01_(rng_) * n_warmup_);
            Eigen::VectorXd newdir = warmup_.row(pi).transpose();
            return step(center_, newdir - center_, -1.0, tries + 1);
        }

        return p;
    }

    // Matches cobra/sampling/optgp.py _sample_chain()
    std::vector<Eigen::VectorXd> sample(int N) {
        if (N <= 0) return {};

        std::vector<Eigen::VectorXd> out(N);

        // Python: pi = np.random.randint(sampler.n_warmup)
        //         prev = sampler.warmup[pi, :]
        //         prev = step(sampler, center, prev - center, 0.95)
        int pi = static_cast<int>(uni01_(rng_) * n_warmup_);
        Eigen::VectorXd prev = warmup_.row(pi).transpose();
        prev = step(center_, prev - center_, 0.95);

        // Python: n_samples = max(sampler.n_samples, 1)
        n_samples_ = std::max(n_samples_, 1);

        // Python: for i in range(1, sampler.thinning * n + 1):
        for (int s = 0; s < N; s++) {
            for (int k = 0; k < thinning_; k++) {
                // Python: pi = np.random.randint(sampler.n_warmup)
                //         delta = sampler.warmup[pi, :] - center
                pi = static_cast<int>(uni01_(rng_) * n_warmup_);
                Eigen::VectorXd delta = warmup_.row(pi).transpose() - center_;

                // Python: prev = step(sampler, prev, delta)
                prev = step(prev, delta);

                // Reprojection: only when homogeneous (matches Python exactly)
                // Python: if sampler.problem.homogeneous and
                //             (n_samples * sampler.thinning % sampler.nproj == 0):
                if (prob_.homogeneous && ((n_samples_ * thinning_) % nproj_ == 0)) {
                    prev = reproject(prev);
                    center_ = reproject(center_);
                }

                // Python: center = (n_samples * center) / (n_samples + 1) + prev / (n_samples + 1)
                center_ = (n_samples_ * center_) / (n_samples_ + 1) + prev / (n_samples_ + 1);
                n_samples_++;
            }
            // Python: if i % sampler.thinning == 0: samples[...] = prev
            out[s] = prev;
        }

        return out;
    }

private:
    SamplerProblem prob_;
    Eigen::MatrixXd warmup_;
    int n_warmup_;
    int thinning_;
    int nproj_;
    int n_vars_;
    int retries_;
    int n_samples_;
    Eigen::VectorXd center_;
    Eigen::FullPivLU<Eigen::MatrixXd> decomp_;

    std::mt19937 rng_;
    std::uniform_real_distribution<double> uni01_;

    // Matches Python hr_sampler._reproject()
    Eigen::VectorXd reproject(const Eigen::VectorXd& v) const {
        if (prob_.equalities.rows() == 0) return v;
        Eigen::VectorXd residual = prob_.equalities * v - prob_.beq;
        // Python: if np.allclose(equalities.dot(p), self.problem.b, ...): new = p
        if (residual.cwiseAbs().maxCoeff() < FEASIBILITY_TOL)
            return v;
        // Python: new = nulls.dot(nulls.T.dot(p))
        // Equivalent: v' = v - Aeq^T * (Aeq * Aeq^T)^-1 * (Aeq * v - b)
        return v - prob_.equalities.transpose() * decomp_.solve(residual);
    }

    // Matches Python hr_sampler._is_redundant() — removes correlated warmup rows
    Eigen::MatrixXd removeRedundant(const Eigen::MatrixXd& warmup) const {
        if (warmup.rows() <= 2) return warmup;

        int n = warmup.rows();
        int d = warmup.cols();

        // Center rows for correlation
        Eigen::MatrixXd centered(n, d);
        Eigen::VectorXd norms(n);
        for (int i = 0; i < n; i++) {
            double mean = warmup.row(i).mean();
            centered.row(i) = warmup.row(i).array() - mean;
            norms(i) = centered.row(i).norm();
        }

        // Mark redundant rows (abs correlation > 0.9999)
        std::vector<bool> keep(n, true);
        for (int i = 0; i < n; i++) {
            if (!keep[i]) continue;
            if (norms(i) < 1e-12) { keep[i] = false; continue; }
            for (int j = i + 1; j < n; j++) {
                if (!keep[j]) continue;
                if (norms(j) < 1e-12) { keep[j] = false; continue; }
                double corr = centered.row(i).dot(centered.row(j)) / (norms(i) * norms(j));
                if (std::abs(corr) > 0.9999)
                    keep[j] = false;
            }
        }

        int kept = 0;
        for (bool k : keep) if (k) kept++;

        Eigen::MatrixXd result(kept, d);
        int idx = 0;
        for (int i = 0; i < n; i++) {
            if (keep[i]) result.row(idx++) = warmup.row(i);
        }

        std::cout << "Warmup: removed " << (n - kept) << " redundant points, keeping " << kept << std::endl;
        return result;
    }
};


// ##### Entry point — same extern "C" interface, no TS changes needed #####
const double VALIDATION_TOL = 1e-2;

extern "C"{
    double sample(const int samplesCount, const int thinning, const int reactionsCount, const int metabolitesCount, float* lbs, float* ubs, float* sData, int initOptRows, float* initOptData,
        float* resultSamples);
}

EMSCRIPTEN_KEEPALIVE
double sample(const int samplesCount, const int thinning, const int reactionsCount, const int metabolitesCount, float* lbs, float* ubs, float* sData, int initOptRows, float* initOptData,
float* resultSamples) {
    const int n_vars = 2 * reactionsCount;

    // === Transform from flux space (n) to variable space (2n) ===
    // Layout: [forward_0..forward_{n-1}, reverse_0..reverse_{n-1}]
    // flux_i = forward_i - reverse_i

    // Variable bounds (all non-negative in variable space)
    Eigen::VectorXd lb_var(n_vars), ub_var(n_vars);
    for (int i = 0; i < reactionsCount; i++) {
        double rxn_lb = lbs[i];
        double rxn_ub = ubs[i];
        // Forward variable: carries positive flux
        lb_var(i) = std::max(0.0, rxn_lb);
        ub_var(i) = std::max(0.0, rxn_ub);
        // Reverse variable: carries negative flux
        lb_var(reactionsCount + i) = std::max(0.0, -rxn_ub);
        ub_var(reactionsCount + i) = std::max(0.0, -rxn_lb);
    }

    // Stoichiometric matrix in variable space: [S, -S]
    // S * flux = 0  =>  S * (fwd - rev) = 0  =>  [S, -S] * [fwd; rev] = 0
    Eigen::MatrixXd S_var(metabolitesCount, n_vars);
    for (int i = 0; i < metabolitesCount; i++) {
        for (int j = 0; j < reactionsCount; j++) {
            double val = sData[i * reactionsCount + j];
            S_var(i, j) = val;                      // forward columns
            S_var(i, reactionsCount + j) = -val;     // reverse columns
        }
    }

    // Identify fixed variables (lb == ub)
    std::vector<bool> variable_fixed(n_vars, false);
    for (int i = 0; i < n_vars; i++) {
        if (std::abs(ub_var(i) - lb_var(i)) < FEASIBILITY_TOL)
            variable_fixed[i] = true;
    }

    // Add fixed non-zero variables as equality constraints
    // (matches Python hr_sampler.__build_problem: fixed_non_zero handling)
    std::vector<int> fixed_nonzero;
    for (int i = 0; i < n_vars; i++) {
        if (variable_fixed[i] && std::abs(ub_var(i)) > FEASIBILITY_TOL)
            fixed_nonzero.push_back(i);
    }

    Eigen::MatrixXd equalities;
    Eigen::VectorXd beq;
    bool homogeneous;

    if (fixed_nonzero.empty()) {
        equalities = S_var;
        beq = Eigen::VectorXd::Zero(metabolitesCount);
        homogeneous = true;
    } else {
        int n_extra = fixed_nonzero.size();
        equalities = Eigen::MatrixXd::Zero(metabolitesCount + n_extra, n_vars);
        equalities.topRows(metabolitesCount) = S_var;
        for (int k = 0; k < n_extra; k++)
            equalities(metabolitesCount + k, fixed_nonzero[k]) = 1.0;

        beq = Eigen::VectorXd::Zero(metabolitesCount + n_extra);
        for (int k = 0; k < n_extra; k++)
            beq(metabolitesCount + k) = ub_var(fixed_nonzero[k]);

        homogeneous = false;
    }

    // Transform warmup (InitOpt) from flux space to variable space
    Eigen::MatrixXd warmup(initOptRows, n_vars);
    for (int i = 0; i < initOptRows; i++) {
        for (int j = 0; j < reactionsCount; j++) {
            double flux = initOptData[i * reactionsCount + j];
            warmup(i, j) = std::max(0.0, flux);                   // forward
            warmup(i, reactionsCount + j) = std::max(0.0, -flux);  // reverse
        }
    }

    // Build sampler problem
    SamplerProblem prob;
    prob.equalities = equalities;
    prob.beq = beq;
    prob.variable_lb = lb_var;
    prob.variable_ub = ub_var;
    prob.variable_fixed = variable_fixed;
    prob.homogeneous = homogeneous;

    // Use fixed seed matching Python's OptGPSampler(model, thinning, seed=42)
    // In Python _sample_chain: np.random.seed((sampler._seed + idx) % np.iinfo(np.int32).max)
    // With idx=0 (single process), effective seed = 42
    unsigned int seed = 42;

    OptGPSampler sampler(prob, warmup, thinning, seed);

    auto start = std::chrono::high_resolution_clock::now();
    auto samples = sampler.sample(samplesCount);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    // === Transform results back to flux space and validate ===
    Eigen::MatrixXd S(metabolitesCount, reactionsCount);
    for (int i = 0; i < metabolitesCount; i++) {
        for (int j = 0; j < reactionsCount; j++)
            S(i, j) = sData[i * reactionsCount + j];
    }

    for (int i = 0; i < samplesCount; i++) {
        // Convert: flux_j = forward_j - reverse_j
        Eigen::VectorXd flux(reactionsCount);
        for (int j = 0; j < reactionsCount; j++)
            flux(j) = samples[i](j) - samples[i](reactionsCount + j);

        // Validate S*flux ≈ 0
        double mad = (S * flux).norm();
        if (mad > VALIDATION_TOL) {
            std::cout << "Sample " << i << " has S*v norm: " << mad << "\n";
            throw std::runtime_error("Invalid sample: " + std::to_string(i) + ", too big deviation: " + std::to_string(mad));
        }

        // Validate flux within reaction bounds
        for (int k = 0; k < reactionsCount; k++) {
            if ((flux(k) - ubs[k] > VALIDATION_TOL) || (flux(k) - lbs[k] < -VALIDATION_TOL)) {
                std::cout << i << ", " << flux(k) - ubs[k] << ", " << flux(k) - lbs[k] << "\n";
                throw std::runtime_error("Invalid sample: value is out of range");
            }
        }

        // Store flux result
        for (int k = 0; k < reactionsCount; k++)
            resultSamples[i * reactionsCount + k] = flux(k);
    }

    return elapsed.count();
}


#endif // !SAMPLER_CPP
