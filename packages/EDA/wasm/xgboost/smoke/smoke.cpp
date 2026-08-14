// Smoke test for the XGBOOST_WASM_MINIMAL build of XGBoost
// (see ../patches/xgboost-wasm-minimal-v3.3.0.patch):
// the three kept objectives must train and predict, the trimmed ones must fail.
//
// Build against a native minimal build of XGBoost (see CMakeLists.txt here),
// or compile with emcc (add -DSMOKE_SKIP_NEGATIVE) and run under node to test
// the wasm build itself.
#include <xgboost/c_api.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#define CHECK_OK(call)                                                    \
  do {                                                                    \
    if ((call) != 0) {                                                    \
      std::printf("FAIL %s:%d: %s\n  %s\n", __FILE__, __LINE__, #call,    \
                  XGBGetLastError());                                     \
      return 1;                                                           \
    }                                                                     \
  } while (0)

namespace {

constexpr int kRows = 120;
constexpr int kCols = 4;

// Deterministic pseudo-random features; label depends on feature sums.
void MakeData(std::vector<float>* x, std::vector<float>* y, int n_class) {
  x->resize(kRows * kCols);
  y->resize(kRows);
  unsigned s = 42;
  for (int i = 0; i < kRows; ++i) {
    float sum = 0.f;
    for (int j = 0; j < kCols; ++j) {
      s = s * 1103515245u + 12345u;
      float v = static_cast<float>((s >> 16) & 0x7fff) / 32768.0f;
      (*x)[i * kCols + j] = v;
      sum += v;
    }
    if (n_class <= 0)
      (*y)[i] = sum;  // regression target
    else
      (*y)[i] = static_cast<float>(static_cast<int>(sum * n_class / kCols) % n_class);
  }
}

int RunObjective(char const* objective, int n_class, float* out_first_pred) {
  std::vector<float> x, y;
  MakeData(&x, &y, n_class);

  DMatrixHandle dmat;
  CHECK_OK(XGDMatrixCreateFromMat(x.data(), kRows, kCols, -1.0f, &dmat));
  CHECK_OK(XGDMatrixSetFloatInfo(dmat, "label", y.data(), kRows));

  BoosterHandle booster;
  DMatrixHandle mats[1] = {dmat};
  CHECK_OK(XGBoosterCreate(mats, 1, &booster));
  CHECK_OK(XGBoosterSetParam(booster, "objective", objective));
  if (n_class > 2)
    CHECK_OK(XGBoosterSetParam(booster, "num_class", std::to_string(n_class).c_str()));
  CHECK_OK(XGBoosterSetParam(booster, "max_depth", "3"));

  for (int it = 0; it < 10; ++it)
    CHECK_OK(XGBoosterUpdateOneIter(booster, it, dmat));

  // Save to buffer and reload: round trip must work.
  char const config[] = "{\"format\": \"ubj\"}";
  bst_ulong len = 0;
  char const* buf = nullptr;
  CHECK_OK(XGBoosterSaveModelToBuffer(booster, config, &len, &buf));
  std::vector<char> model(buf, buf + len);
  CHECK_OK(XGBoosterFree(booster));

  BoosterHandle loaded;
  CHECK_OK(XGBoosterCreate(nullptr, 0, &loaded));
  CHECK_OK(XGBoosterLoadModelFromBuffer(loaded, model.data(), len));

  bst_ulong out_len = 0;
  float const* preds = nullptr;
  CHECK_OK(XGBoosterPredict(loaded, dmat, 0, 0, 0, &out_len, &preds));

  int bad = 0;
  for (bst_ulong i = 0; i < out_len; ++i) {
    float p = preds[i];
    if (p != p) ++bad;                                     // NaN
    if (n_class == 2 && (p < 0.f || p > 1.f)) ++bad;       // probability range
    if (n_class > 2 && (p < 0.f || p > n_class - 1)) ++bad;  // class index range
  }
  *out_first_pred = out_len > 0 ? preds[0] : -1.f;
  std::printf("OK   %-18s rows=%d preds=%llu model=%llu bytes%s\n", objective, kRows,
              static_cast<unsigned long long>(out_len),
              static_cast<unsigned long long>(len), bad ? " [BAD VALUES]" : "");

  CHECK_OK(XGBoosterFree(loaded));
  CHECK_OK(XGDMatrixFree(dmat));
  return bad ? 1 : 0;
}

// A trimmed component must fail: error either on SetParam or on first update.
int ExpectFailure(char const* param, char const* value) {
  std::vector<float> x, y;
  MakeData(&x, &y, 0);
  DMatrixHandle dmat;
  CHECK_OK(XGDMatrixCreateFromMat(x.data(), kRows, kCols, -1.0f, &dmat));
  CHECK_OK(XGDMatrixSetFloatInfo(dmat, "label", y.data(), kRows));
  BoosterHandle booster;
  DMatrixHandle mats[1] = {dmat};
  CHECK_OK(XGBoosterCreate(mats, 1, &booster));

  int rc = XGBoosterSetParam(booster, param, value);
  if (rc == 0)
    rc = XGBoosterUpdateOneIter(booster, 0, dmat);

  XGBoosterFree(booster);
  XGDMatrixFree(dmat);
  if (rc == 0) {
    std::printf("FAIL %s=%s unexpectedly available (should be trimmed)\n", param, value);
    return 1;
  }
  std::printf("OK   %s=%s correctly unavailable\n", param, value);
  return 0;
}

}  // namespace

int main() {
  int failures = 0;
  float p;

  failures += RunObjective("reg:squarederror", 0, &p);
  failures += RunObjective("binary:logistic", 2, &p);
  failures += RunObjective("multi:softmax", 3, &p);

// The negative checks rely on C++ exceptions reaching the C API boundary.
// The production wasm build has exceptions disabled (any throw aborts the
// instance), so run them only where exceptions work (native build).
#ifndef SMOKE_SKIP_NEGATIVE
  failures += ExpectFailure("objective", "rank:ndcg");
  failures += ExpectFailure("objective", "reg:quantileerror");
  failures += ExpectFailure("objective", "binary:hinge");
  failures += ExpectFailure("booster", "gblinear");
  failures += ExpectFailure("tree_method", "exact");
  failures += ExpectFailure("tree_method", "approx");
#endif

  std::printf(failures ? "\nSMOKE TEST FAILED (%d)\n" : "\nSMOKE TEST PASSED\n", failures);
  return failures;
}
