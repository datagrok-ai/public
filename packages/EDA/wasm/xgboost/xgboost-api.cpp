// Thin wasm wrapper over the XGBoost C API for the EDA package.
//
// Design constraints (see xgboost-update-plan.md, phase 2):
//  - Input is columnar float32, exactly as Datagrok stores columns: one
//    contiguous buffer where column j occupies [j*nRows, (j+1)*nRows).
//    Columns are passed to XGBoost via the array-interface JSON, so no
//    transposition and no extra copy happens here.
//  - The wasm build has C++ exceptions disabled: any internal XGBoost error
//    aborts the instance. All argument/label validation MUST happen on the
//    TypeScript side before calling into this module.
//  - Models live in a handle table and survive between calls; prediction
//    never re-parses model bytes.
#include <xgboost/c_api.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

namespace {

// Objective codes shared with the TypeScript side.
enum Objective : int {
  kRegression = 0,  // reg:squarederror
  kBinary = 1,      // binary:logistic (predictions are probabilities)
  kMulticlass = 2,  // multi:softmax   (predictions are class indices)
};

struct ModelSlot {
  BoosterHandle booster{nullptr};
  std::string bytes;  // serialized model cache (filled by xgbModelSize)
};

// Slot 0 is never used so that a zero handle always means "error".
std::vector<ModelSlot>& Slots() {
  static std::vector<ModelSlot> slots{ModelSlot{}};
  return slots;
}

ModelSlot* GetSlot(int handle) {
  auto& slots = Slots();
  if (handle <= 0 || static_cast<size_t>(handle) >= slots.size()) return nullptr;
  if (slots[handle].booster == nullptr) return nullptr;
  return &slots[handle];
}

int StoreBooster(BoosterHandle booster) {
  auto& slots = Slots();
  for (size_t i = 1; i < slots.size(); ++i) {
    if (slots[i].booster == nullptr) {
      slots[i] = ModelSlot{booster, {}};
      return static_cast<int>(i);
    }
  }
  slots.push_back(ModelSlot{booster, {}});
  return static_cast<int>(slots.size() - 1);
}

// JSON array-interface for one float32 column at `data` with `nRows` values.
void AppendColumnInterface(std::string* out, const float* data, int nRows) {
  char buf[160];
  std::snprintf(buf, sizeof(buf),
                "{\"data\":[%lu,true],\"shape\":[%d],\"typestr\":\"<f4\",\"version\":3}",
                static_cast<unsigned long>(reinterpret_cast<uintptr_t>(data)), nRows);
  *out += buf;
}

// JSON list of array-interfaces: one per column of the col-major buffer.
std::string ColumnarInterface(const float* colMajorData, int nRows, int nCols) {
  std::string json;
  json.reserve(static_cast<size_t>(nCols) * 96 + 2);
  json += '[';
  for (int j = 0; j < nCols; ++j) {
    if (j > 0) json += ',';
    AppendColumnInterface(&json, colMajorData + static_cast<size_t>(j) * nRows, nRows);
  }
  json += ']';
  return json;
}

std::string FloatToJson(float v) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.9g", v);
  return buf;
}

int SetParam(BoosterHandle booster, const char* name, const std::string& value) {
  return XGBoosterSetParam(booster, name, value.c_str());
}

}  // namespace

extern "C" {

// Train a model. Returns a model handle (> 0) on success, 0 on failure.
// colMajorData: nRows*nCols float32, column j at [j*nRows, (j+1)*nRows).
// labels: nRows float32. For kBinary labels are 0/1; for kMulticlass labels
// are integral values in [0, numClass) — validated on the TypeScript side.
EMSCRIPTEN_KEEPALIVE
int xgbTrain(const float* colMajorData, int nRows, int nCols, float missing,
             const float* labels,
             int objective, int numClass,
             int iterations, float eta, int maxDepth, float lambda, float alpha) {
  const std::string data = ColumnarInterface(colMajorData, nRows, nCols);
  const std::string config = "{\"missing\":" + FloatToJson(missing) + ",\"nthread\":1}";

  DMatrixHandle dmat = nullptr;
  if (XGDMatrixCreateFromColumnar(data.c_str(), config.c_str(), &dmat) != 0) return 0;
  if (XGDMatrixSetFloatInfo(dmat, "label", labels, nRows) != 0) {
    XGDMatrixFree(dmat);
    return 0;
  }

  BoosterHandle booster = nullptr;
  DMatrixHandle mats[1] = {dmat};
  int rc = XGBoosterCreate(mats, 1, &booster);

  if (rc == 0) {
    switch (objective) {
      case kBinary:
        rc = XGBoosterSetParam(booster, "objective", "binary:logistic");
        break;
      case kMulticlass:
        rc = XGBoosterSetParam(booster, "objective", "multi:softmax");
        if (rc == 0) rc = SetParam(booster, "num_class", std::to_string(numClass));
        break;
      default:
        rc = XGBoosterSetParam(booster, "objective", "reg:squarederror");
    }
  }
  if (rc == 0) rc = SetParam(booster, "eta", FloatToJson(eta));
  if (rc == 0) rc = SetParam(booster, "max_depth", std::to_string(maxDepth));
  if (rc == 0) rc = SetParam(booster, "lambda", FloatToJson(lambda));
  if (rc == 0) rc = SetParam(booster, "alpha", FloatToJson(alpha));

  for (int it = 0; rc == 0 && it < iterations; ++it)
    rc = XGBoosterUpdateOneIter(booster, it, dmat);

  XGDMatrixFree(dmat);
  if (rc != 0) {
    if (booster != nullptr) XGBoosterFree(booster);
    return 0;
  }
  return StoreBooster(booster);
}

// Serialize the model (UBJSON) into the slot cache; returns byte count or -1.
EMSCRIPTEN_KEEPALIVE
int xgbModelSize(int modelHandle) {
  ModelSlot* slot = GetSlot(modelHandle);
  if (slot == nullptr) return -1;

  bst_ulong len = 0;
  const char* buf = nullptr;
  if (XGBoosterSaveModelToBuffer(slot->booster, "{\"format\":\"ubj\"}", &len, &buf) != 0)
    return -1;
  slot->bytes.assign(buf, len);
  return static_cast<int>(len);
}

// Copy the serialized model (cached by xgbModelSize) into dst; returns byte
// count or -1 (no cache / dst too small).
EMSCRIPTEN_KEEPALIVE
int xgbModelCopy(int modelHandle, uint8_t* dst, int dstSize) {
  ModelSlot* slot = GetSlot(modelHandle);
  if (slot == nullptr || slot->bytes.empty()) return -1;
  if (dstSize < static_cast<int>(slot->bytes.size())) return -1;
  std::memcpy(dst, slot->bytes.data(), slot->bytes.size());
  return static_cast<int>(slot->bytes.size());
}

// Load a model from serialized bytes. Returns a model handle (> 0) or 0.
EMSCRIPTEN_KEEPALIVE
int xgbLoadModel(const uint8_t* bytes, int size) {
  BoosterHandle booster = nullptr;
  if (XGBoosterCreate(nullptr, 0, &booster) != 0) return 0;
  if (XGBoosterLoadModelFromBuffer(booster, bytes, size) != 0) {
    XGBoosterFree(booster);
    return 0;
  }
  return StoreBooster(booster);
}

// Inplace prediction over columnar data; no DMatrix is built.
// out must hold nRows floats (regression: value; binary: probability of
// class 1; multiclass: class index). Returns 0 on success, -1 on failure.
EMSCRIPTEN_KEEPALIVE
int xgbPredict(int modelHandle, const float* colMajorData, int nRows, int nCols,
               float missing, float* out, int outLen) {
  ModelSlot* slot = GetSlot(modelHandle);
  if (slot == nullptr || outLen < nRows) return -1;

  const std::string data = ColumnarInterface(colMajorData, nRows, nCols);
  const std::string config =
      "{\"type\":0,\"training\":false,\"iteration_begin\":0,\"iteration_end\":0,"
      "\"strict_shape\":false,\"missing\":" + FloatToJson(missing) + "}";

  bst_ulong const* outShape = nullptr;
  bst_ulong outDim = 0;
  float const* result = nullptr;
  if (XGBoosterPredictFromColumnar(slot->booster, data.c_str(), config.c_str(), nullptr,
                                   &outShape, &outDim, &result) != 0)
    return -1;

  bst_ulong n = 1;
  for (bst_ulong d = 0; d < outDim; ++d) n *= outShape[d];
  if (n != static_cast<bst_ulong>(nRows)) return -1;

  std::memcpy(out, result, sizeof(float) * n);
  return 0;
}

EMSCRIPTEN_KEEPALIVE
void xgbFreeModel(int modelHandle) {
  ModelSlot* slot = GetSlot(modelHandle);
  if (slot == nullptr) return;
  XGBoosterFree(slot->booster);
  slot->booster = nullptr;
  slot->bytes.clear();
  slot->bytes.shrink_to_fit();
}

}  // extern "C"
