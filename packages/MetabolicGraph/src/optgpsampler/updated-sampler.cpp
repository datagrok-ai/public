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

#include <iostream>
#include <chrono>

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

    Eigen::MatrixXd InitOpt;
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

    double spreadRatio = 0.7;
};

class OptGPSampler {
public:
    OptGPSampler(const LinearPolytope& model, const OptGPOptions& opt = {})
        : A_(model.A), b_(model.b), Aeq_(model.Aeq), beq_(model.beq),
        lb_(model.lb), ub_(model.ub), opt_(opt), rng_(opt.seed), norm01_(0.0, 1.0), InitOpt(model.InitOpt)
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

        //decomp = (Aeq_ * Aeq_.transpose()).ldlt();

        decomp = (Aeq_ * Aeq_.transpose()).fullPivLu();

        std::cout << "rank(A * A_tr): " << decomp.rank() << std::endl;

        // Initialize x: use x0 if feasible, else find one
        if (model.x0 && isFeasible(*model.x0)) {
            x_ = *model.x0;
        }
        else {
            //Eigen::VectorXd x(n_);

            // mean
            //x << 2.88696837413413, 2.82100848365597, 2.79721931959462, 2.51076098626129, 5.30798030585591, 1.65593085470754, 8.29781909375042, 0.070493213472858, 0.399393690745588, 0.0868055555555555, 9.14106825347296, 8.60800609468012, 1.36517010719548, 1.45142386336462, -7.4870821177433, 8.30331575017093, 0.0759898698933729, 2.45193276260344, 2.74764865732274, 0.0694444444444445, -8.86251611121454, 0.416839941682045, 0.225093354979181, -8.41808157122955, -9.10049099589946, -3.09363120287671, 0.0184739731457025, 0.054780627725467, -1.23808106787509, 0.12230214225842, -9.65616263910463, -9.56820511753623, 1.32025712825345, -0.379311803799166, 0.114039081163163, 0.337093444163699, 0.119235592018175, 0.20125358471021, 0.0618872846693813, 0.0346926195800632, 0.173611111111111, 0.0347222222222222, 0.138888888888889, 0.069394276472952, 0.607777230499511, 0.377978384242748, 0.229798846256763, 0.528222214291086, 0.460453857955481, 0.0820179926920349, -9.54894282114226, -9.39626262299088, 0.043402777777774, 0.297095696236722, 3.6504811263257, -8.1396295149746, 0.185917832481748, 0.260589941682053, 0.347395497237601, -0.0845145940387369, 0.00984248483843235, 0.00164041413973873, 0.00765863224965369, 0.00765863224965369, 0.00601936048491425, 0.0247148014191648, 0.0129038196130459, 0.0118109818061188, 0.0300658263063023, 0.0300658263063023, 0.0144394682350605, 0.00313163823774526, 0.002202841844792, 0.00534306319800613, 0.00314022135321413, 0.194406666502026, 0.194374778477652, 0.000641291711640513, 0.104265018948946, 0.0184103131338416, 0.00621346483144762, 0.00467160513281025, -9.55327383997151, 0.00287933879139529, 0.00287933879139529, 0.000759427211332416, 0, 0.00200211655244578, 0.00814808795018311, 0.234344877105532, 2.72025666641734, 0.286458333333334, 9.81692366695694, 0.542780487800658, 0.00120461017076561, 0.0108750982876006, 0.0129140657737514, 0.001568324822668, 0.00226487791947547, 0.00102020276635806, 6.42040759194236E-06, 0.00117197683351213, 9.04810852842805, 0.720486111111111, 2.88696837413413, 9.72067882985945, -9.10049099589946, 0.342096119890277, -4.65124633602191, -3.28275713748243, 0.0312967019755127, -9.4261147726933, 0.0212385248716719, -9.8323985042735, 8.65848990724926, 9.37899953131425, 0.356071728881824, 0.167436616133692, 0.83984278812672, 0.291582147183778, 0.127361240012397, 0.630418271200858, -9.13345183792803, -3.25107807341801, 9.78369465834856, 0.0732629346501976, 0.0437787017516424, 0.407483946975114, 0.312588135190033, 0.759214220056152, 8.08126822968287, 0.444434539984996, 0.138888888888889, 0.229798846256763, 2.72025666641734, 0.656819568963821, -8.91043968150091, 9.81692366695694, -9.55327383997151, 0.0858547058151044, 0.185917832481748, 0.234344877105532, 0.00117197683351213, 1.38310185185185;

            //x << 10, 10, 10, 0, 10, 10, 10, 0, 2E-15, 0, 10, 10, 0, 1.2E-15, -10, 10, 0, -2.8E-14, -5.68E-14, 0, -10, 0, 0, -10, -10, 1.4E-15, 4E-16, -1.2E-15, -4E-16, 0, -10, -10, -1.3E-14, 2.84E-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2E-15, 1.2E-15, 0, 0, 0, 0, -10, -10, 0, 0, 0, -10, 0, 0, 0, 2E-16, -2E-16, 0, 0, 0, 0, -4E-16, -4E-16, -2E-16, -6E-16, -6E-16, -4E-16, 0, 0, 0, 0, -3.2E-15, -3.2E-15, 2E-15, 0, -4E-16, 0, 0, -8E-16, 0, 0, 0, 0, 0, -2E-16, -4.6E-15, 0, 10, 10, -2.92E-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, -10, 2.4E-15, -10, 6E-15, -1.4E-15, -10, -2.8E-15, -10, 10, 10, -6E-16, -1.2E-15, -1.2E-15, -2.4E-15, -6E-16, -1.2E-15, -10, 2.52E-14, 10, 1.8E-15, -8E-16, -2E-15, 2E-15, -2E-15, 10, 0, 0, 0, 0, -2.92E-14, -3E-14, 10, -8E-16, 4E-16, 0, -4.6E-15, 0, 0;

            //x_ = x;

            //x_ = findFeasiblePoint_(model.x0);

            x_ = InitOpt.transpose().rowwise().mean();
        }

        std::cout << "Initial point: " << x_.transpose() << std::endl;
    }

    // Generate N samples (each is an Eigen::VectorXd of size n)
    std::vector<Eigen::VectorXd> sample(int N) {
        if (N <= 0) return {};

        warmup_();

        std::vector<Eigen::VectorXd> out(N);

        for (int s = 0; s < N; ++s) {
            for (int j = 0; j < opt_.nproj; ++j) {
                for (int k = 0; k < opt_.thinning; ++k) x_ = singleStep_(x_);
            }

            out[s] = x_;

            /*if (uni01_(rng_) <= opt_.spreadRatio)
                x_ = InitOpt.row(static_cast<int>(uni01_(rng_) * InitOpt.rows())).transpose();*/

            double coef = -0.5 + static_cast<int>(uni01_(rng_));

            x_ = InitOpt.transpose().rowwise().mean();
            x_ = x_ + (InitOpt.row(static_cast<int>(uni01_(rng_) * InitOpt.rows())).transpose() - x_) * coef;
                
                ; //coef * x_ + (1 - coef) * InitOpt.row(static_cast<int>(uni01_(rng_) * InitOpt.rows())).transpose();
        }

        /*for (int s = 0; s < N - InitOpt.rows(); ++s) {
            for (int j = 0; j < opt_.nproj; ++j) {
                for (int k = 0; k < opt_.thinning; ++k) x_ = singleStep_(x_);
            }

            out[s] = x_;
        }

        for (int s = 0; s < InitOpt.rows(); ++s) {
            x_ = InitOpt.row(s).transpose();

            for (int j = 0; j < opt_.nproj; ++j) {
                for (int k = 0; k < opt_.thinning; ++k) x_ = singleStep_(x_);
            }

            out[N - InitOpt.rows() + s] = x_;
        }*/

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

    Eigen::MatrixXd InitOpt;

    //Eigen::LDLT<Eigen::MatrixXd> decomp;
    Eigen::FullPivLU<Eigen::MatrixXd> decomp;

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
        //Eigen::VectorXd y = AAT.ldlt().solve(beq_);
        Eigen::VectorXd y = AAT.fullPivLu().solve(beq_);
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
                Eigen::VectorXd y = (Aeq_ * Aeq_.transpose()).fullPivLu().solve(r);
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

        std::cout << "center: " << center_.value() << "\n";

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