// Flux sampler for the MetabolicGraph package, compiled to WebAssembly (see emscc.bat).
//
// It reproduces cobrapy's OptGPSampler (cobra 0.29 - 0.32) with processes=1 exactly: for the same
// model, seed, sample count and thinning it returns the same samples as
//     OptGPSampler(model, thinning=thinning, processes=1, seed=seed).sample(samplesCount)
// That takes three parts, each a port of the matching cobra/optlang/numpy code path:
//   1. warmup points: cobra's FVA warmup, solved with GLPK 5.0 exactly the way cobra + optlang call
//      it (cobraWarmup below);
//   2. the OptGP chain: hit-and-run steps along (warmup point - running center) directions in the
//      split space of forward/reverse variables, with cobra's step, retry and reprojection rules;
//   3. random numbers: numpy's legacy RandomState draws over the same MT19937 stream.
// Like cobra, it samples the space of forward and reverse variables (fwd_i, rev_i >= 0,
// S (fwd - rev) = 0) and reports net fluxes fwd - rev.
#include <glpk.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#ifndef SAMPLER_OUT_T
#define SAMPLER_OUT_T float
#endif
typedef SAMPLER_OUT_T OutT;

namespace {

// cobra: feasibility_tol = bounds_tol = model.tolerance (1e-7 with GLPK)
const double TOL = 1e-7;
const int MAX_TRIES = 100;

// Error codes returned by the entry points (non-negative results are elapsed milliseconds)
const double ERR_SINGLE_POINT = -1;   // the warmup collapses to one point: nothing to sample
const double ERR_TWO_DIRECTIONS = -2; // inhomogeneous problem with only 2 search directions
const double ERR_UNSTABLE = -3;       // cannot escape the sampling region (numerically unstable)
const double ERR_INFEASIBLE = -4;     // no warmup LP could be solved: infeasible bounds
const double ERR_GLPK = -5;           // GLPK failed to copy the problem

// numpy's legacy RandomState draws (np.random.seed / randint / uniform) over the same MT19937
// stream as std::mt19937, which is what makes the chain identical to cobra's
class NumpyRandom {
public:
    explicit NumpyRandom(uint32_t seed) : mt_(seed) {}
    // mt19937_next_double: a 53-bit double from two 32-bit draws
    double uniform01() {
        const uint32_t a = mt_() >> 5, b = mt_() >> 6;
        return (a * 67108864.0 + b) / 9007199254740992.0;
    }
    // RandomState.randint(n): masked rejection sampling on 32-bit draws
    int randint(int n) {
        const uint32_t rng = (uint32_t)(n - 1);
        if (rng == 0) return 0;
        uint32_t mask = rng;
        mask |= mask >> 1; mask |= mask >> 2; mask |= mask >> 4; mask |= mask >> 8; mask |= mask >> 16;
        uint32_t v;
        while ((v = (uint32_t)mt_() & mask) > rng) {}
        return (int)v;
    }
private:
    std::mt19937 mt_;
};

// cobra Reaction.update_variable_bounds: bounds of the forward and reverse variables
void splitBounds(double lb, double ub, double& fLo, double& fHi, double& rLo, double& rHi) {
    if (lb <= 0 && 0 <= ub) { fLo = 0; fHi = ub; rLo = 0; rHi = -lb; }
    else if (lb > 0) { fLo = lb; fHi = ub; rLo = 0; rHi = 0; }
    else { fLo = 0; fHi = 0; rLo = -ub; rHi = -lb; }
}

// ---------------------------------------------------------------------------------------------
// 1. Warmup: cobra HRSampler.generate_fva_warmup through optlang's GLPK interface
// ---------------------------------------------------------------------------------------------
namespace cobra_warmup {

enum Status { OPTIMAL, OTHER, UNDEFINED };

// optlang glpk_interface Model._run_glp_simplex (only optimal / undefined / other matter here)
Status runSimplex(glp_prob* P, glp_smcp* parm) {
    if (glp_simplex(P, parm) == GLP_ETMLIM) return OTHER;
    const int st = glp_get_status(P);
    return st == GLP_OPT ? OPTIMAL : st == GLP_UNDEF ? UNDEFINED : OTHER;
}

// optlang glpk_interface Model._optimize, with presolve "auto" (off) or True (on)
Status optimizeOnce(glp_prob* P, glp_smcp* parm, bool presolveOn) {
    glp_scale_prob(P, GLP_SF_AUTO);
    parm->presolve = presolveOn ? GLP_ON : GLP_OFF;
    Status st = runSimplex(P, parm);
    if (st == UNDEFINED && !presolveOn) {
        glp_adv_basis(P, 0);
        st = runSimplex(P, parm);
    }
    if (st == UNDEFINED && presolveOn) {
        parm->presolve = GLP_OFF;
        st = runSimplex(P, parm);
    }
    parm->presolve = GLP_OFF;
    return st;
}

// optlang interface Model.optimize: a non-optimal result is retried once with presolve on
Status optimize(glp_prob* P, glp_smcp* parm) {
    const Status st = optimizeOnce(P, parm, false);
    return st == OPTIMAL ? st : optimizeOnce(P, parm, true);
}

// optlang's column bounds: GLPK has no infinite bounds, so those become one-sided/free columns
void setColBounds(glp_prob* P, int j, double lb, double ub) {
    const bool noLb = std::isinf(lb), noUb = std::isinf(ub);
    if (noLb && noUb) glp_set_col_bnds(P, j, GLP_FR, 0., 0.);
    else if (noLb) glp_set_col_bnds(P, j, GLP_UP, 0., ub);
    else if (noUb) glp_set_col_bnds(P, j, GLP_LO, lb, 0.);
    else if (lb == ub) glp_set_col_bnds(P, j, GLP_FX, lb, lb);
    else glp_set_col_bnds(P, j, GLP_DB, lb, ub);
}

// Appends the raw warmup points (before the redundancy filter) to `out`, in the split layout
// [fwd_0..fwd_{n-1}, rev_0..rev_{n-1}]. Returns the number of points, or ERR_GLPK.
//  - The problem is built as optlang builds it from a cobra model: columns fwd_0, rev_0, fwd_1, ...
//    and one fixed row per metabolite (S rows in model.metabolites order), each set in ascending
//    column order. GLPK's pivoting depends on that layout, so it has to match.
//  - It is then copied through GLPK's file format, as cobra's model.copy() pickles optlang models;
//    the copy rounds all numbers to 15 significant digits and cobra solves on the copy.
//  - Every non-fixed reaction is minimized, then maximized, as fwd - rev, with warm-started
//    simplex runs (scaling before each run, presolve off, optlang's fallbacks).
int generate(int n, int m, const double* lbs, const double* ubs, const double* S, std::vector<double>& out) {
#ifdef __EMSCRIPTEN__
    const char* copyPath = "/tmp/cobra-model.glpk";
#else
    const char* copyPath = "cobra-model.glpk";
#endif
    glp_term_out(GLP_OFF);
    glp_prob* orig = glp_create_prob();
    glp_add_rows(orig, m);
    glp_add_cols(orig, 2 * n);
    for (int k = 0; k < n; k++) {
        double fLo, fHi, rLo, rHi;
        splitBounds(lbs[k], ubs[k], fLo, fHi, rLo, rHi);
        setColBounds(orig, 2 * k + 1, fLo, fHi);
        setColBounds(orig, 2 * k + 2, rLo, rHi);
    }
    std::vector<int> ind(2 * (size_t)n + 1);
    std::vector<double> val(2 * (size_t)n + 1);
    for (int i = 0; i < m; i++) {
        glp_set_row_bnds(orig, i + 1, GLP_FX, 0., 0.);
        int len = 0;
        for (int k = 0; k < n; k++) {
            const double s = S[(size_t)i * n + k];
            if (s == 0) continue;
            ind[++len] = 2 * k + 1; val[len] = s;
            ind[++len] = 2 * k + 2; val[len] = -s;
        }
        glp_set_mat_row(orig, i + 1, len, ind.data(), val.data());
    }
    glp_set_obj_dir(orig, GLP_MAX);

    glp_prob* P = glp_create_prob();
    const bool copied = glp_write_prob(orig, 0, copyPath) == 0 && glp_read_prob(P, 0, copyPath) == 0;
    glp_delete_prob(orig);
    std::remove(copyPath);
    if (!copied) { glp_delete_prob(P); return (int)ERR_GLPK; }

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    parm.presolve = GLP_OFF;
    parm.tm_lim = INT_MAX;

    int rows = 0;
    glp_set_obj_dir(P, GLP_MAX); // model.objective = Zero
    for (int sense : {GLP_MIN, GLP_MAX}) {
        glp_set_obj_dir(P, sense);
        for (int k = 0; k < n; k++) {
            if (ubs[k] - lbs[k] < TOL) continue; // cobra skips fixed reactions
            glp_set_obj_coef(P, 2 * k + 1, 1.);
            glp_set_obj_coef(P, 2 * k + 2, -1.);
            if (optimize(P, &parm) != OPTIMAL) continue; // cobra leaves these coefficients set
            const size_t base = out.size();
            out.resize(base + 2 * (size_t)n);
            for (int j = 0; j < n; j++) {
                out[base + j] = glp_get_col_prim(P, 2 * j + 1);
                out[base + n + j] = glp_get_col_prim(P, 2 * j + 2);
            }
            rows++;
            glp_set_obj_coef(P, 2 * k + 1, 0.);
            glp_set_obj_coef(P, 2 * k + 2, 0.);
        }
    }
    glp_delete_prob(P);
    return rows;
}

} // namespace cobra_warmup

// ---------------------------------------------------------------------------------------------
// 2. The OptGP chain
// ---------------------------------------------------------------------------------------------

// cobra HRSampler.__build_problem for a model without extra constraints
struct Problem {
    int nVars = 0;                 // 2 * reactions, layout [fwd..., rev...]
    int rows = 0;
    std::vector<double> A;         // equalities, row-major (rows x nVars)
    std::vector<double> b;
    std::vector<double> lb, ub;    // variable bounds
    std::vector<double> slb, sub;  // (1 - TOL) * bounds: the targets of a step
    std::vector<char> fixed;       // ub - lb < TOL
    bool homogeneous = true;
};

Problem buildProblem(int reactions, int metabolites, const double* lbs, const double* ubs, const double* S) {
    Problem p;
    p.nVars = 2 * reactions;
    p.lb.resize(p.nVars); p.ub.resize(p.nVars);
    for (int j = 0; j < reactions; j++)
        splitBounds(lbs[j], ubs[j], p.lb[j], p.ub[j], p.lb[reactions + j], p.ub[reactions + j]);
    p.fixed.resize(p.nVars);
    p.slb.resize(p.nVars); p.sub.resize(p.nVars);
    std::vector<int> fixedNonZero;
    for (int i = 0; i < p.nVars; i++) {
        p.fixed[i] = p.ub[i] - p.lb[i] < TOL;
        p.slb[i] = (1 - TOL) * p.lb[i];
        p.sub[i] = (1 - TOL) * p.ub[i];
        if (p.fixed[i] && std::abs(p.ub[i]) > TOL) fixedNonZero.push_back(i);
    }
    // fixed non-zero variables become equalities (and make the problem inhomogeneous)
    p.rows = metabolites + (int)fixedNonZero.size();
    p.A.assign((size_t)p.rows * p.nVars, 0.0);
    p.b.assign(p.rows, 0.0);
    for (int r = 0; r < metabolites; r++) {
        for (int j = 0; j < reactions; j++) {
            const double s = S[(size_t)r * reactions + j];
            p.A[(size_t)r * p.nVars + j] = s;
            p.A[(size_t)r * p.nVars + reactions + j] = -s;
        }
    }
    for (size_t k = 0; k < fixedNonZero.size(); k++) {
        const int r = metabolites + (int)k, i = fixedNonZero[k];
        p.A[(size_t)r * p.nVars + i] = 1;
        p.b[r] = p.ub[i];
    }
    p.homogeneous = fixedNonZero.empty();
    return p;
}

class OptGP {
public:
    OptGP(const Problem& prob, std::vector<double> warmup, int nWarmup, int thinning, unsigned seed)
        : p_(prob), n_(prob.nVars), thinning_(thinning), rng_(seed),
          delta_(n_), next_(n_), prev_(n_), center_(n_), center0_(n_) {
        nproj_ = (int64_t)std::min<double>(std::pow((double)n_, 3), 1e6); // cobra: min(n_vars^3, 1e6)
        warmup_ = std::move(warmup);
        nWarmup_ = nWarmup;
        dropRedundant();
        if (nWarmup_ <= 1) { status_ = ERR_SINGLE_POINT; return; }
        if (nWarmup_ == 2) {
            if (!p_.homogeneous) { status_ = ERR_TWO_DIRECTIONS; return; }
            // cobra: "all search directions on a line, adding another one"
            warmup_.resize((size_t)3 * n_);
            for (int i = 0; i < n_; i++) warmup_[2 * (size_t)n_ + i] = 0.25 * (warmup_[i] + warmup_[n_ + i]);
            nWarmup_ = 3;
        }
        // numpy's mean over rows: sequential sums, then one division
        for (int k = 0; k < nWarmup_; k++)
            for (int i = 0; i < n_; i++) center0_[i] += warmup_[(size_t)k * n_ + i];
        for (int i = 0; i < n_; i++) center0_[i] /= nWarmup_;
    }

    double status() const { return status_; }
    int retries() const { return retries_; }
    int warmupPoints() const { return nWarmup_; }

    // cobra optgp._sample_chain (chain 0); writes the net flux of every recorded sample
    bool sample(int count, int reactions, OutT* out) {
        center_ = center0_;
        std::copy_n(&warmup_[(size_t)randint(nWarmup_) * n_], n_, prev_.begin());
        for (int i = 0; i < n_; i++) delta_[i] = prev_[i] - center_[i];
        if (!step(&center_, 0.95)) return false;
        prev_.swap(next_);

        int64_t nSamples = 1;
        const int64_t total = (int64_t)thinning_ * count;
        for (int64_t it = 1; it <= total; it++) {
            const double* w = &warmup_[(size_t)randint(nWarmup_) * n_];
            for (int i = 0; i < n_; i++) delta_[i] = w[i] - center_[i];
            if (!step(&prev_, -1)) return false;
            prev_.swap(next_);
            if (p_.homogeneous && (nSamples * thinning_) % nproj_ == 0) {
                reproject(prev_);
                reproject(center_);
            }
            if (it % thinning_ == 0) {
                OutT* row = out + (size_t)(it / thinning_ - 1) * reactions;
                for (int j = 0; j < reactions; j++) row[j] = (OutT)(prev_[j] - prev_[reactions + j]);
            }
            // cobra's operation order, so the chain stays bit-identical
            const double ns = (double)nSamples, ns1 = (double)(nSamples + 1);
            for (int i = 0; i < n_; i++) center_[i] = (ns * center_[i]) / ns1 + prev_[i] / ns1;
            nSamples++;
        }
        return true;
    }

private:
    const Problem& p_;
    int n_;
    int thinning_;
    int64_t nproj_ = 0;
    std::vector<double> warmup_;
    int nWarmup_ = 0;
    int retries_ = 0;
    double status_ = 0;
    NumpyRandom rng_;
    std::vector<double> delta_, next_, prev_, center_, center0_;

    int randint(int n) { return rng_.randint(n); }

    // cobra core.step: hit-and-run move from x along delta_ into next_. fraction < 0 means a uniform
    // position on the feasible segment. An infeasible or stuck move restarts from the initial
    // center towards a random warmup point, as cobra's recursive retry does.
    bool step(const std::vector<double>* x, double fraction) {
        for (int tries = 0;; tries++) {
            double lo = 0, hi = 0;
            bool hasLo = false, hasHi = false;
            for (int i = 0; i < n_; i++) {
                const double d = delta_[i];
                if (std::abs(d) <= TOL || p_.fixed[i]) continue;
                const double xi = (*x)[i];
                const double a1 = (p_.slb[i] - xi) / d, a2 = (p_.sub[i] - xi) / d;
                if (a1 > 0) { if (!hasHi || a1 < hi) hi = a1; hasHi = true; }
                else { if (!hasLo || a1 > lo) lo = a1; hasLo = true; }
                if (a2 > 0) { if (!hasHi || a2 < hi) hi = a2; hasHi = true; }
                else { if (!hasLo || a2 > lo) lo = a2; hasLo = true; }
            }
            const double alpha = fraction >= 0 ? lo + fraction * (hi - lo) : lo + (hi - lo) * rng_.uniform01();
            const double reach = std::max(std::abs(lo), std::abs(hi));
            double lbDist = std::numeric_limits<double>::infinity(), ubDist = lbDist, move = 0;
            for (int i = 0; i < n_; i++) {
                const double v = (*x)[i] + alpha * delta_[i];
                next_[i] = v;
                lbDist = std::min(lbDist, v - p_.lb[i]);
                ubDist = std::min(ubDist, p_.ub[i] - v);
                move = std::max(move, std::abs(reach * delta_[i]));
            }
            if (lbDist >= -TOL && ubDist >= -TOL && move >= TOL) return true;
            if (tries > MAX_TRIES) { status_ = ERR_UNSTABLE; return false; }
            retries_++;
            const double* w = &warmup_[(size_t)randint(nWarmup_) * n_];
            for (int i = 0; i < n_; i++) delta_[i] = w[i] - center0_[i];
            x = &center0_;
            fraction = -1;
        }
    }

    // cobra _reproject: a point that violates the equalities is replaced by a random point
    // (the mean of two random warmup points); the nullspace projection it computes is never kept
    void reproject(std::vector<double>& v) {
        for (int r = 0; r < p_.rows; r++) {
            double s = -p_.b[r];
            for (int i = 0; i < n_; i++) s += p_.A[(size_t)r * n_ + i] * v[i];
            if (std::abs(s) > TOL) {
                const double* w1 = &warmup_[(size_t)randint(nWarmup_) * n_];
                const double* w2 = &warmup_[(size_t)randint(nWarmup_) * n_];
                for (int i = 0; i < n_; i++) v[i] = (w1[i] + w2[i]) / 2;
                return;
            }
        }
    }

    // cobra _is_redundant: drop a warmup point whose Pearson correlation with any earlier point
    // exceeds 1 - TOL (an extra column keeps constant and all-zero points comparable)
    void dropRedundant() {
        const int m = n_ + 1;
        std::vector<double> X((size_t)nWarmup_ * m), norm(nWarmup_);
        for (int k = 0; k < nWarmup_; k++) {
            const double* w = &warmup_[(size_t)k * n_];
            double* x = &X[(size_t)k * m];
            double sum = 0;
            for (int i = 0; i < n_; i++) { x[i] = w[i]; sum += w[i]; }
            x[n_] = sum == 0 ? 2 : w[0] + 1;
            double mean = 0;
            for (int i = 0; i < m; i++) mean += x[i];
            mean /= m;
            double ss = 0;
            for (int i = 0; i < m; i++) { x[i] -= mean; ss += x[i] * x[i]; }
            norm[k] = std::sqrt(ss);
        }
        std::vector<char> redundant(nWarmup_, 0);
        for (int k = 1; k < nWarmup_; k++) {
            for (int j = 0; j < k && !redundant[k]; j++) {
                if (norm[k] == 0 || norm[j] == 0) continue; // numpy yields NaN, which never compares true
                double dot = 0;
                for (int i = 0; i < m; i++) dot += X[(size_t)k * m + i] * X[(size_t)j * m + i];
                if (std::abs(dot / (norm[k] * norm[j])) > 1 - TOL) redundant[k] = 1;
            }
        }
        int kept = 0;
        for (int k = 0; k < nWarmup_; k++) {
            if (redundant[k]) continue;
            if (kept != k) std::copy_n(&warmup_[(size_t)k * n_], n_, &warmup_[(size_t)kept * n_]);
            kept++;
        }
        nWarmup_ = kept;
        warmup_.resize((size_t)kept * n_);
    }
};

double run(int samplesCount, int thinning, int seed, int reactions, int metabolites, const double* lbs,
           const double* ubs, const double* S, std::vector<double> warmup, int warmupRows, OutT* result, int* stats) {
    const auto start = std::chrono::steady_clock::now();
    const Problem prob = buildProblem(reactions, metabolites, lbs, ubs, S);
    OptGP sampler(prob, std::move(warmup), warmupRows, std::max(1, thinning), (unsigned)seed);
    if (stats) { stats[1] = warmupRows; stats[2] = sampler.warmupPoints(); }
    if (sampler.status() < 0) return sampler.status();
    const bool ok = sampler.sample(samplesCount, reactions, result);
    if (stats) stats[0] = sampler.retries();
    if (!ok) return sampler.status();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

} // namespace

extern "C" {

// Samples `samplesCount` flux vectors into `result` (samplesCount x reactionsCount, row-major),
// exactly as cobra's OptGPSampler(model, thinning, processes=1, seed).sample(samplesCount).
//   lbs, ubs: reaction bounds (reactionsCount)
//   S:        stoichiometric matrix, row-major (metabolitesCount x reactionsCount), rows in the
//             order of the model's metabolites (the warmup LPs depend on that order)
//   stats:    optional int[3] receiving [retries, raw warmup points, warmup points used]
// Returns the elapsed time in ms, or a negative error code.
EMSCRIPTEN_KEEPALIVE
double sample(int samplesCount, int thinning, int seed, int reactionsCount, int metabolitesCount,
              const double* lbs, const double* ubs, const double* S, OutT* result, int* stats) {
    std::vector<double> warmup;
    const int rows = cobra_warmup::generate(reactionsCount, metabolitesCount, lbs, ubs, S, warmup);
    if (rows < 0) return rows;
    if (rows == 0) return ERR_INFEASIBLE;
    return run(samplesCount, thinning, seed, reactionsCount, metabolitesCount, lbs, ubs, S,
               std::move(warmup), rows, result, stats);
}

// Same, with precomputed warmup points in the split layout [fwd..., rev...] (warmupRows x
// 2*reactionsCount), e.g. cobra's own warmup from the ComputeExtremePoints Python script.
EMSCRIPTEN_KEEPALIVE
double sample_with_warmup(int samplesCount, int thinning, int seed, int reactionsCount, int metabolitesCount,
                          const double* lbs, const double* ubs, const double* S, int warmupRows,
                          const double* warmup, OutT* result, int* stats) {
    if (warmupRows <= 0) return ERR_INFEASIBLE;
    std::vector<double> w(warmup, warmup + (size_t)warmupRows * 2 * reactionsCount);
    return run(samplesCount, thinning, seed, reactionsCount, metabolitesCount, lbs, ubs, S,
               std::move(w), warmupRows, result, stats);
}

}
