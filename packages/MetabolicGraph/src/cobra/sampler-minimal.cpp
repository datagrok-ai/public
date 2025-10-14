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
#include <emscripten.h>

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
    int thinning = 20;     // steps between effective moves per internal step
    int nproj = 1000000;          // short chain length per returned sample
    int warmupSteps = 5000; // steps used to estimate a center
    unsigned int seed = 1234; // random seed based on time

    // Feasible-point finder
    double feasTol = 1e-7;
    int feasMaxIter = 5000;
    double feasIneqWeight = 1.0;
    double feasStep = 0.95;  // initial step size
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


        // Initialize x: use x0 if feasible, else find one
        if (model.x0 && isFeasible(*model.x0)) {
            x_ = *model.x0;
        }
        else {
            x_ = InitOpt.transpose().rowwise().mean();
        }
    }

    void remove_redundent_from_initOpt(const Eigen::MatrixXd &M) {
        // this will remove correlated points from initOpt
        Eigen::MatrixXd corr = rowCorrelation(M);
        std::vector<int> to_remove;
        int n = corr.rows();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(corr(i, j)) > 0.995) {
                    to_remove.push_back(j);
                }
            }
        }
        std::sort(to_remove.begin(), to_remove.end());
        to_remove.erase(std::unique(to_remove.begin(), to_remove.end()), to_remove.end());
        Eigen::MatrixXd newM(n - to_remove.size(), M.cols());
        int idx = 0;
        for (int i = 0; i < n; ++i) {
            if (std::find(to_remove.begin(), to_remove.end(), i) == to_remove.end()) {
                newM.row(idx++) = M.row(i);
            }
        }
        std::cout << "Removed " << to_remove.size() << " correlated points from InitOpt." << std::endl;
        InitOpt = newM;

    }

    Eigen::MatrixXd rowCorrelation(const Eigen::MatrixXd &M) {
        int n = M.rows();
        Eigen::MatrixXd corr(n, n);
        corr.setZero();

        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd xi = M.row(i).transpose();
            double mean_i = xi.mean();
            double std_i = std::sqrt((xi.array() - mean_i).square().sum());

            for (int j = 0; j <= i; ++j) { // compute only lower-triangular
                Eigen::VectorXd xj = M.row(j).transpose();
                double mean_j = xj.mean();
                double std_j = std::sqrt((xj.array() - mean_j).square().sum());

                if (std_i == 0.0 || std_j == 0.0) {
                    corr(i, j) = 0.0; // handle constant rows
                } else {
                    double dot = ((xi.array() - mean_i) * (xj.array() - mean_j)).sum();
                    corr(i, j) = dot / (std_i * std_j);
                }
                corr(j, i) = corr(i, j); // symmetry
            }
        }
        return corr;
    }

    Eigen::VectorXd step(Eigen::VectorXd x, Eigen::VectorXd delta) {
        double maxAlpha = std::numeric_limits<double>::infinity();
        double minAlpha = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n_; ++i) {
            double di = delta(i);
            double ci = x(i);
            double li = lb_(i);
            double ui = ub_(i);
            if (std::abs(di) < 1e-12) {
                // direction component is zero, skip
                continue;
            }
            double maxUpperBoundCoef = (ui - ci) / di; // here we know for sure that di != 0
            double minLowerBoundCoef = (li - ci) / di;
            if (maxUpperBoundCoef > 0)
                maxAlpha = std::min(maxAlpha, maxUpperBoundCoef);
            else
                minAlpha = std::max(minAlpha, maxUpperBoundCoef); // although this is not possible since min <= max

            if (minLowerBoundCoef > 0)// although this is not possible since min <= max
                maxAlpha = std::min(maxAlpha, minLowerBoundCoef);
            else
                minAlpha = std::max(minAlpha, minLowerBoundCoef);
        }
                // make sure min and max and not infinite or NaN
        if (std::isinf(minAlpha)) minAlpha = 0;
        if (std::isinf(maxAlpha)) maxAlpha = 0;
        if (std::isnan(minAlpha)) minAlpha = 0;
        if (std::isnan(maxAlpha)) maxAlpha = 0;
        maxAlpha = std::max(0.0, maxAlpha);
        minAlpha = std::min(0.0, minAlpha);
        if (maxAlpha == 0 && minAlpha == 0) {
            // direction is zero, jump to random point
            // although again, this is not possible :D
            std::cout << "Warning: direction is zero!" << std::endl;
            maxAlpha = 1.0; 
        }
        
        // now sample coef in [minAlpha, maxAlpha]
        double rnd = uni01_(rng_);
        double coef = minAlpha + (maxAlpha - minAlpha) * rnd;
        return x + delta * coef;

    }

    std::vector<Eigen::VectorXd> sampleCB(int N) {
        if (N <= 0) return {};

        //warmup_();

        std::vector<Eigen::VectorXd> out(N);
        // remove redundant points from InitOpt
        remove_redundent_from_initOpt(InitOpt);

        Eigen::VectorXd center = InitOpt.transpose().rowwise().mean();
        Eigen::VectorXd prev = InitOpt.row(static_cast<int>(uni01_(rng_) * InitOpt.rows())).transpose();
        // reseed random
        rng_.seed(N);
        prev = step(center, prev - center);
        int n_samples = 1;
        for (int s = 0; s < N; ++s) {
            for (int k = 0; k < opt_.thinning; ++k) {
                Eigen::VectorXd targetPoint = InitOpt.row(static_cast<int>(uni01_(rng_) * InitOpt.rows())).transpose();
                auto delta = targetPoint - center;
                prev = step(prev, delta);
                int iter = s * opt_.thinning + k;
                if ((iter * opt_.thinning) % opt_.nproj == 0 && iter > 0) {
                    std::cout << "Updating center, sample: " << s << ", step: " << k << "  " << opt_.nproj <<  std::endl;
                    // update center estimate
                    prev = projectToNull_(prev);
                    center = projectToNull_(center);
                }
 
                center = (n_samples * center) / (n_samples + 1) + prev / (n_samples + 1);
                n_samples += 1; 
            }
            out[s] = prev;
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

};


// ##### sampler caller function #####
const double DEFAULT_UPPER_BOUND = 10;
const double TOL = 1e-2;

extern "C"{
    double sample(const int samplesCount, const int thinning, const int reactionsCount, const int metabolitesCount, float* lbs, float* ubs, float* sData, int initOptRows, float* initOptData,
        float* resultSamples);
}

EMSCRIPTEN_KEEPALIVE
double sample(const int samplesCount, const int thinning, const int reactionsCount, const int metabolitesCount, float* lbs, float* ubs, float* sData, int initOptRows, float* initOptData,
float* resultSamples) {
    Eigen::VectorXd beq(metabolitesCount);
	beq.setConstant(0.0);
    Eigen::VectorXd lb(reactionsCount);

    double val = 0;

	for (int i = 0; i < reactionsCount; ++i) {		
		lb(i) = lbs[i];
	}

	Eigen::VectorXd ub(reactionsCount);
	for (int i = 0; i < reactionsCount; ++i) {		
		ub(i) = ubs[i];
	}
	//ub.setConstant(2.0);
    // S matrix or the stoichiometric matrix. each row represents a metabolite, each column a reaction
	Eigen::MatrixXd S(metabolitesCount, reactionsCount);
	
	for (int i = 0; i < metabolitesCount; ++i) {
		for (int j = 0; j < reactionsCount; ++j) {
			S(i, j) = sData[i * reactionsCount + j];	
        }	
    }

	int initOptRowCount = 288;
	int initOptColCount = 144;

	Eigen::MatrixXd InitOpt(initOptRows, reactionsCount);
	//lb(24) = 0.0; // lactate exchange reaction
	for (int i = 0; i < initOptRows; ++i) {
		for (int j = 0; j < reactionsCount; ++j) {
			InitOpt(i, j) = initOptData[i * reactionsCount + j];	
		}
	}

    LinearPolytope P;
	P.A = Eigen::MatrixXd();
	P.b = Eigen::VectorXd();
	P.Aeq = S; P.beq = beq;
	P.lb = lb; P.ub = ub;
	P.InitOpt = InitOpt;

	OptGPOptions opt;
	//opt.seed = 1983;
	opt.seed = 1234 * time(0) % std::numeric_limits<int>::max();
	opt.thinning = thinning;
	opt.nproj = 1000000;
	// opt.warmupSteps = WARM_STEPS;
	// opt.feasTol = FEAS_TOL;
	// opt.feasMaxIter = FEAS_MAX_ITER;
	// opt.spreadRatio = 0.5;
    OptGPSampler sampler(P, opt);

    auto start = std::chrono::high_resolution_clock::now();

	auto samples = sampler.sampleCB(samplesCount);

	auto end = std::chrono::high_resolution_clock::now();


    std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds	

	double mad = 0;

	for (int i = 0; i < samplesCount; ++i) {
		mad = (S * samples[i]).norm();
		if (mad > TOL) {
            std::cout << "Sample " << i << " has S*v norm: " << mad << "\n";
			throw std::runtime_error("Invalid sample: " + std::to_string(i) + ", too big deviation: " + std::to_string(mad));
		}

		for (int k = 0; k < reactionsCount; ++k) {
			if ((samples[i](k) - ub(k) > 1e-5) || (samples[i](k) - lb(k) < -1e-5)) {
				std::cout << samples[i](k) - ub(k) << ", " << samples[i](k) - lb(k) << "\n";
				throw std::runtime_error("Invalid sample: value is out of range");
			}
		}
        // store the sample
        for (int k = 0; k < reactionsCount; ++k) {
            resultSamples[i * reactionsCount + k] = samples[i](k);
        }
	}

	return elapsed.count();
}


// int main() {
//     std::cout << samplerChompact(100000)/* samplerRandomTest(10, 20, 1) */ << " ms" << " --- " << "SUCCESS!\n";
//     return 0;
// }

#endif // !SAMPLER_CPP