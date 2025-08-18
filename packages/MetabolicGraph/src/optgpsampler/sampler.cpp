// OptGPSampler (C++/Eigen)
// ------------------------------------------------------------
// A practical, single-file implementation inspired by COBRApy's
// OptGPSampler, using Eigen for linear algebra. Includes an
// automatic feasible-point finder so x0 is optional.
//
// Dependencies: Eigen 3
//   #include <Eigen/Dense>
//
// Compile (example):
//   g++ -O3 -std=c++17 -I/path/to/eigen3 optgp_sampler.hpp test.cpp -o test
//
// Notes:
//  • Feasible region:
//        A x <= b
//        Aeq x = beq
//        lb <= x <= ub
//  • Equality feasibility is preserved by projecting directions
//    into Null(Aeq).
//  • Uniform line sampling inside feasible segment per step.
//  • "thinning" and short chains (nproj) reduce autocorrelation.
//  • RNG is seedable.
//  • The feasible-point finder is a Phase‑I style heuristic: least-squares
//    equality solve + alternating projections (eq/bounds) + inequality
//    violation reduction with projected gradient/backtracking.
//    It is not a full LP solver, but works well in practice.

#ifndef SAMPLER_CPP
#define SAMPLER_CPP

//#include <Eigen/Dense>
#include "../../Eigen/Eigen/Dense"
#include <random>
#include <vector>
#include <optional>
#include <stdexcept>
#include <limits>
#include <cmath>

//#include <iostream>
//#include <chrono>

struct LinearPolytope {
    // Inequalities: A x <= b
    Eigen::MatrixXd A;      // (mi x n) or empty
    Eigen::VectorXd b;      // (mi)

    // Equalities: Aeq x = beq
    Eigen::MatrixXd Aeq;    // (me x n) or empty
    Eigen::VectorXd beq;    // (me)

    // Bounds: lb <= x <= ub
    Eigen::VectorXd lb;     // (n) (use -inf for no lower bound)
    Eigen::VectorXd ub;     // (n) (use +inf for no upper bound)

    // Optional feasible starting point
    std::optional<Eigen::VectorXd> x0;
};

struct OptGPOptions {
    int thinning = 100;     // steps between effective moves per internal step
    int nproj = 8;          // short chain length per returned sample
    int warmupSteps = 5000; // steps used to estimate a center
    unsigned int seed = 1234;

    // Feasible-point finder
    double feasTol = 1e-7;
    int feasMaxIter = 5000;
    double feasIneqWeight = 1.0;
    double feasStep = 1.0;  // initial step size
};

class OptGPSampler {
public:
    OptGPSampler(const LinearPolytope& model, const OptGPOptions& opt = {})
        : A_(model.A), b_(model.b), Aeq_(model.Aeq), beq_(model.beq),
        lb_(model.lb), ub_(model.ub), opt_(opt), rng_(opt.seed), norm01_(0.0, 1.0)
    {
        const int n = inferN_(model);
        n_ = n;
        // Dimension checks
        if (lb_.size() == 0) lb_ = Eigen::VectorXd::Constant(n_, -std::numeric_limits<double>::infinity());
        if (ub_.size() == 0) ub_ = Eigen::VectorXd::Constant(n_, std::numeric_limits<double>::infinity());
        if (A_.rows() > 0 && A_.cols() != n_) throw std::invalid_argument("A cols != n");
        if (A_.rows() > 0 && b_.size() != A_.rows()) throw std::invalid_argument("b size != A rows");
        if (Aeq_.rows() > 0 && Aeq_.cols() != n_) throw std::invalid_argument("Aeq cols != n");
        if (Aeq_.rows() > 0 && beq_.size() != Aeq_.rows()) throw std::invalid_argument("beq size != Aeq rows");

        decomp = (Aeq_ * Aeq_.transpose()).ldlt();

        // Initialize x: use x0 if feasible, else find one
        if (model.x0 && isFeasible(*model.x0)) {
            x_ = *model.x0;
        }
        else {
            x_ = findFeasiblePoint_(model.x0);
        }
    }

    // Generate N samples (each is an Eigen::VectorXd of size n)
    std::vector<Eigen::VectorXd> sample(int N) {
        if (N <= 0) return {};

        //auto start = std::chrono::high_resolution_clock::now();

        warmup_();

        //auto end = std::chrono::high_resolution_clock::now();

        //std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds

        //std::cout << "Warmup time: " << elapsed.count() << " ms\n";
        

        //std::vector<Eigen::VectorXd> out; out.reserve(N);

        std::vector<Eigen::VectorXd> out(N);

        for (int s = 0; s < N; ++s) {
            //auto start = std::chrono::high_resolution_clock::now();

            for (int j = 0; j < opt_.nproj; ++j) {
                for (int k = 0; k < opt_.thinning; ++k) x_ = singleStep_(x_);
            }

            //auto end = std::chrono::high_resolution_clock::now();

            //std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds

            //std::cout << s << "-th point: " << elapsed.count() << " ms\n";
            //out.push_back(x_);

            out[s] = x_;
        }
        return out;
    }

    int dim() const { return n_; }
    const Eigen::VectorXd& state() const { return x_; }

private:
    // Problem data
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    Eigen::MatrixXd Aeq_;
    Eigen::VectorXd beq_;
    Eigen::VectorXd lb_;
    Eigen::VectorXd ub_;

    Eigen::LDLT<Eigen::MatrixXd> decomp;

    // Options and state
    OptGPOptions opt_;
    int n_ = 0;
    Eigen::VectorXd x_;            // current point
    std::optional<Eigen::VectorXd> center_; // estimated center

    // RNG
    std::mt19937 rng_;
    std::normal_distribution<double> norm01_;
    std::uniform_real_distribution<double> uni01_{ 0.0,1.0 };

private:
    static int inferN_(const LinearPolytope& M) {
        if (M.x0 && M.x0->size() > 0) return int(M.x0->size());
        if (M.lb.size() > 0) return int(M.lb.size());
        if (M.ub.size() > 0) return int(M.ub.size());
        if (M.Aeq.cols() > 0) return int(M.Aeq.cols());
        if (M.A.cols() > 0) return int(M.A.cols());
        throw std::invalid_argument("Cannot infer variable dimension n");
    }

    bool isFeasible(const Eigen::VectorXd& x) const {
        // bounds
        for (int i = 0; i < n_; ++i) {
            if (x[i] < lb_[i] - 1e-8 || x[i] > ub_[i] + 1e-8) return false;
        }
        // inequalities
        if (A_.rows() > 0) {
            Eigen::VectorXd Ax = A_ * x;
            for (int i = 0; i < Ax.size(); ++i) if (Ax[i] - b_[i] > 1e-8) return false;
        }
        // equalities
        if (Aeq_.rows() > 0) {
            Eigen::VectorXd r = Aeq_ * x - beq_;
            if (r.cwiseAbs().maxCoeff() > 1e-8) return false;
        }
        return true;
    }

    // Nullspace projection: v' = v - Aeq^T (Aeq Aeq^T)^{-1} (Aeq v)
    Eigen::VectorXd projectToNull_(const Eigen::VectorXd& v) const {
        if (Aeq_.rows() == 0) return v;
        Eigen::VectorXd Av = Aeq_ * v;                      // m
        // Solve (Aeq Aeq^T) y = Av; use robust solver
        //Eigen::MatrixXd AAT = Aeq_ * Aeq_.transpose();
        //Eigen::VectorXd y = AAT.ldlt().solve(Av);
        return v - Aeq_.transpose() * decomp.solve(Av);
    }

    // Particular solution of Aeq x = beq (minimum-norm)
    Eigen::VectorXd eqParticular_() const {
        if (Aeq_.rows() == 0) return Eigen::VectorXd::Zero(n_);
        // Solve min ||x|| s.t. Aeq x = beq → x = Aeq^T (Aeq Aeq^T)^{-1} beq
        Eigen::MatrixXd AAT = Aeq_ * Aeq_.transpose();
        Eigen::VectorXd y = AAT.ldlt().solve(beq_);
        return Aeq_.transpose() * y;
    }

    // Alternating projection onto equalities and bounds
    Eigen::VectorXd projectEqAndBounds_(const Eigen::VectorXd& v) const {
        Eigen::VectorXd x = v;
        auto clamp = [&](Eigen::VectorXd& z) {
            for (int i = 0; i < n_; ++i) {
                if (std::isfinite(lb_[i]) && z[i] < lb_[i]) z[i] = lb_[i];
                if (std::isfinite(ub_[i]) && z[i] > ub_[i]) z[i] = ub_[i];
            }
            };

        for (int it = 0; it < std::min(200, opt_.feasMaxIter); ++it) {
            clamp(x);
            if (Aeq_.rows() > 0) {
                // project to equality manifold: x' = argmin ||x'-x|| s.t. Aeq x' = beq
                // This is: x' = x - Aeq^T (Aeq Aeq^T)^{-1} (Aeq x - beq)
                Eigen::VectorXd r = Aeq_ * x - beq_;
                Eigen::VectorXd y = (Aeq_ * Aeq_.transpose()).ldlt().solve(r);
                Eigen::VectorXd corr = Aeq_.transpose() * y;
                Eigen::VectorXd xnew = x - corr;
                if ((xnew - x).lpNorm<Eigen::Infinity>() < opt_.feasTol) { x = xnew; break; }
                x = xnew;
            }
            else {
                // Only bounds
                break;
            }
        }
        return x;
    }

    // Quantify total violation (for backtracking comparisons)
    double violation_(const Eigen::VectorXd& x) const {
        double v = 0.0;
        // bounds
        for (int i = 0; i < n_; ++i) {
            if (x[i] < lb_[i]) v += (lb_[i] - x[i]);
            if (x[i] > ub_[i]) v += (x[i] - ub_[i]);
        }
        // inequalities
        if (A_.rows() > 0) {
            Eigen::VectorXd Ax = A_ * x - b_;
            for (int i = 0; i < Ax.size(); ++i) if (Ax[i] > 0) v += Ax[i];
        }
        // equalities
        if (Aeq_.rows() > 0) {
            Eigen::VectorXd r = Aeq_ * x - beq_;
            v += r.cwiseAbs().sum();
        }
        return v;
    }

    Eigen::VectorXd findFeasiblePoint_(const std::optional<Eigen::VectorXd>& seed) const {
        // Start at seed, or equality particular solution, or mid-bounds
        Eigen::VectorXd x(n_);
        if (seed) {
            x = *seed;
        }
        else if (Aeq_.rows() > 0) {
            x = eqParticular_();
        }
        else {
            x.setZero();
            for (int i = 0; i < n_; ++i) {
                if (std::isfinite(lb_[i]) && std::isfinite(ub_[i])) x[i] = 0.5 * (lb_[i] + ub_[i]);
                else x[i] = 0.0;
            }
        }

        // Alternate projections eq/bounds a few times
        x = projectEqAndBounds_(x);

        // Reduce inequality violations by projected gradient in Null(Aeq)
        double step = opt_.feasStep;
        for (int it = 0; it < opt_.feasMaxIter; ++it) {
            bool violated = false;
            Eigen::VectorXd grad = Eigen::VectorXd::Zero(n_);
            if (A_.rows() > 0) {
                Eigen::VectorXd Ax = A_ * x - b_;
                for (int i = 0; i < Ax.size(); ++i) {
                    double val = Ax[i];
                    if (val > 0) {
                        violated = true;
                        grad.noalias() += val * A_.row(i).transpose();
                    }
                }
            }
            if (!violated) {
                x = projectEqAndBounds_(x);
                if (isFeasible(x)) return x;
            }
            if (grad.norm() < 1e-16) break;
            Eigen::VectorXd d = -opt_.feasIneqWeight * grad;
            d = projectToNull_(d);

            // backtracking with clipping
            bool accepted = false;
            for (int bt = 0; bt < 20; ++bt) {
                Eigen::VectorXd cand = x + step * d;
                // clip to box
                for (int i = 0; i < n_; ++i) {
                    if (std::isfinite(lb_[i]) && cand[i] < lb_[i]) cand[i] = lb_[i];
                    if (std::isfinite(ub_[i]) && cand[i] > ub_[i]) cand[i] = ub_[i];
                }
                double vCand = violation_(cand);
                if (vCand < violation_(x) - 1e-12) { x = projectEqAndBounds_(cand); accepted = true; break; }
                step *= 0.5;
            }
            if (!accepted) break;
            if (violation_(x) < opt_.feasTol) { x = projectEqAndBounds_(x); if (isFeasible(x)) return x; }
        }

        if (!isFeasible(x)) throw std::runtime_error("Failed to find a feasible starting point. Provide x0 or loosen constraints.");
        return x;
    }

    void warmup_() {
        if (center_) return;
        Eigen::VectorXd x = x_;
        Eigen::VectorXd mean = Eigen::VectorXd::Zero(n_);
        int cnt = 0;
        for (int t = 0; t < opt_.warmupSteps; ++t) {
            x = singleStep_(x);
            if ((t + 1) % opt_.thinning == 0) {
                ++cnt;
                mean += (x - mean) / double(cnt);
            }
        }
        center_ = (cnt > 0) ? mean : x;
        x_ = x;
    }

    // Compute feasible step interval along direction d from x
    std::pair<double, double> stepInterval_(const Eigen::VectorXd& x, const Eigen::VectorXd& d) const {
        double amin = -std::numeric_limits<double>::infinity();
        double amax = std::numeric_limits<double>::infinity();

        // bounds
        for (int i = 0; i < n_; ++i) {
            if (std::abs(d[i]) < 1e-14) continue;
            double a1 = (lb_[i] - x[i]) / d[i];
            double a2 = (ub_[i] - x[i]) / d[i];
            double lo = std::min(a1, a2), hi = std::max(a1, a2);
            if (lo > amin) amin = lo;
            if (hi < amax) amax = hi;
        }

        // inequalities
        if (A_.rows() > 0) {
            Eigen::VectorXd Ad = A_ * d;
            Eigen::VectorXd Ax = A_ * x;
            for (int i = 0; i < Ad.size(); ++i) {
                double c = Ad[i];
                double rhs = b_[i] - Ax[i];
                if (std::abs(c) < 1e-14) { if (rhs < -1e-12) return { 1.0, 0.0 }; continue; }
                double a = rhs / c;
                if (c > 0) { if (a < amax) amax = a; }
                else { if (a > amin) amin = a; }
            }
        }
        return { amin, amax };
    }

    Eigen::VectorXd singleStep_(const Eigen::VectorXd& x) {
        Eigen::VectorXd c = center_.value_or(x);

        // direction: (c - x) + normal noise, then project to Null(Aeq)

        Eigen::VectorXd g(n_);

        //auto start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < n_; ++i)
            g[i] = norm01_(rng_);

        Eigen::VectorXd d = (c - x) + g;

        d = projectToNull_(d);

        //auto end = std::chrono::high_resolution_clock::now();

        //std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds

        //std::cout << "projectToNull_: " << elapsed.count() << " ms\n";

        //start = std::chrono::high_resolution_clock::now();
        
        auto seg = stepInterval_(x, d);

        //end = std::chrono::high_resolution_clock::now();

        //elapsed = end - start; // In milliseconds

        //std::cout << "stepInterval_: " << elapsed.count() << " ms\n\n";


        

        if (!(seg.first < seg.second) || !std::isfinite(seg.first) || !std::isfinite(seg.second)) {
            // fallback: pure random nullspace direction
            for (int i = 0; i < n_; ++i)
                g[i] = norm01_(rng_);

            d = projectToNull_(g);

            seg = stepInterval_(x, d);

            if (!(seg.first < seg.second)) 
                return x; // degenerate
        }

        double a = seg.first + uni01_(rng_) * (seg.second - seg.first);

        return x + a * d;
    }
};

/* -------------------------- Example usage -----------------------------
#include "optgp_sampler.hpp"
#include <iostream>
int main(){
  // Simple chain Sv=0 with bounds 0<=v<=10
  Eigen::MatrixXd S(2,3); S << 1,-1,0, 0,1,-1; // Aeq
  Eigen::VectorXd beq(2); beq << 0,0;
  Eigen::VectorXd lb(3); lb << 0,0,0;
  Eigen::VectorXd ub(3); ub << 10,10,10;

  LinearPolytope P;
  P.A = Eigen::MatrixXd();
  P.b = Eigen::VectorXd();
  P.Aeq = S; P.beq = beq;
  P.lb = lb; P.ub = ub;
  P.x0 = std::nullopt; // let the sampler find one

  OptGPOptions opt; opt.seed=42; opt.thinning=50; opt.nproj=8; opt.warmupSteps=2000;

  OptGPSampler sampler(P, opt);
  auto samples = sampler.sample(1000);
  std::cout << "N=" << samples.size() << " first sample: " << samples[0].transpose() << "\n";
}
*/

#endif // !SAMPLER_CPP