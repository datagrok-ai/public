// Thin wasm wrapper over the LIBSVM C API. Two LIBSVM-specific points:
//  1. No in-memory (de)serialization: svm_save_model / svm_load_model are
//     file-based, so we round-trip through Emscripten MEMFS (kModelPath).
//  2. svm_train does NOT copy the support vectors (model->SV point into the
//     input buffers), so right after training we save+load through MEMFS to
//     obtain a self-contained model (free_sv=1) before freeing the inputs.
#include "svm.h"

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

// Sequential MEMFS scratch path.
const char* const kModelPath = "/svm.model";

// Silence LIBSVM's SMO progress output (stdout is meaningless in wasm).
void svmQuiet(const char*) {}

struct ModelSlot {
  svm_model* model{nullptr};
  std::string bytes;  // serialized model cache (LIBSVM text format)
};

// Slot 0 is never used so that a zero handle always means "error".
std::vector<ModelSlot>& Slots() {
  static std::vector<ModelSlot> slots(1);
  return slots;
}

ModelSlot* GetSlot(int handle) {
  auto& slots = Slots();
  if (handle <= 0 || static_cast<size_t>(handle) >= slots.size()) return nullptr;
  if (slots[handle].model == nullptr) return nullptr;
  return &slots[handle];
}

int StoreModel(svm_model* model, std::string bytes) {
  auto& slots = Slots();
  for (size_t i = 1; i < slots.size(); ++i) {
    if (slots[i].model == nullptr) {
      slots[i].model = model;
      slots[i].bytes = std::move(bytes);
      return static_cast<int>(i);
    }
  }
  slots.push_back(ModelSlot{model, std::move(bytes)});
  return static_cast<int>(slots.size() - 1);
}

std::string readFile(const char* path) {
  FILE* f = fopen(path, "rb");
  if (f == nullptr) return {};
  std::string s;
  if (fseek(f, 0, SEEK_END) == 0) {
    long n = ftell(f);
    if (n > 0) {
      s.resize(static_cast<size_t>(n));
      fseek(f, 0, SEEK_SET);
      size_t rd = fread(&s[0], 1, static_cast<size_t>(n), f);
      s.resize(rd);
    }
  }
  fclose(f);
  return s;
}

bool writeFile(const char* path, const uint8_t* data, int n) {
  FILE* f = fopen(path, "wb");
  if (f == nullptr) return false;
  size_t wr = fwrite(data, 1, static_cast<size_t>(n), f);
  return (fclose(f) == 0) && (wr == static_cast<size_t>(n));
}

}  // namespace

extern "C" {

// Train a model → handle (>0) or 0. Columnar float32 in; unused params ignored.
EMSCRIPTEN_KEEPALIVE
int svmTrain(const float* colMajor, int nRows, int nCols, const float* labels,
             int svmType, int kernelType, double C, double gamma, int degree,
             double coef0, double epsilon, double nu, double eps,
             double cacheSize, int shrinking, int probability) {
  if (nRows <= 0 || nCols <= 0) return 0;

  svm_set_print_string_function(&svmQuiet);

  const int stride = nCols + 1;  // features + (-1) terminator
  std::vector<svm_node> xSpace(static_cast<size_t>(nRows) * stride);
  std::vector<svm_node*> x(nRows);
  std::vector<double> y(nRows);
  for (int i = 0; i < nRows; ++i) {
    svm_node* row = &xSpace[static_cast<size_t>(i) * stride];
    for (int j = 0; j < nCols; ++j) {
      row[j].index = j + 1;  // LIBSVM indices are 1-based
      row[j].value = colMajor[static_cast<size_t>(j) * nRows + i];
    }
    row[nCols].index = -1;
    x[i] = row;
    y[i] = labels[i];
  }

  svm_problem prob;
  prob.l = nRows;
  prob.y = y.data();
  prob.x = x.data();

  svm_parameter param;
  param.svm_type = svmType;
  param.kernel_type = kernelType;
  param.degree = degree;
  param.gamma = (gamma == 0.0) ? 1.0 / nCols : gamma;  // 0 => 1/num_features
  param.coef0 = coef0;
  param.cache_size = cacheSize;  // MB
  param.eps = eps;               // stopping tolerance
  param.C = C;
  param.nr_weight = 0;
  param.weight_label = nullptr;
  param.weight = nullptr;
  param.nu = nu;
  param.p = epsilon;  // epsilon in the SVR loss
  param.shrinking = shrinking;
  param.probability = probability;

  if (svm_check_parameter(&prob, &param) != nullptr) return 0;

  svm_model* trained = svm_train(&prob, &param);
  if (trained == nullptr) return 0;

  // Round-trip so the stored model owns its SVs (svm_train references xSpace,
  // which is freed when this function returns). save() must run while xSpace
  // is still alive; free_sv=0 means destroy does NOT free the SV nodes.
  std::string bytes;
  const bool saved = (svm_save_model(kModelPath, trained) == 0);
  if (saved) bytes = readFile(kModelPath);
  svm_free_and_destroy_model(&trained);
  if (!saved || bytes.empty()) return 0;

  svm_model* model = svm_load_model(kModelPath);
  if (model == nullptr) return 0;
  return StoreModel(model, std::move(bytes));
}

// Serialize the model into the slot cache; returns byte count or -1.
EMSCRIPTEN_KEEPALIVE
int svmModelSize(int handle) {
  ModelSlot* slot = GetSlot(handle);
  if (slot == nullptr) return -1;
  if (slot->bytes.empty()) {
    if (svm_save_model(kModelPath, slot->model) != 0) return -1;
    slot->bytes = readFile(kModelPath);
    if (slot->bytes.empty()) return -1;
  }
  return static_cast<int>(slot->bytes.size());
}

// Copy serialized model into dst; returns bytes or -1.
EMSCRIPTEN_KEEPALIVE
int svmModelCopy(int handle, uint8_t* dst, int dstSize) {
  ModelSlot* slot = GetSlot(handle);
  if (slot == nullptr || slot->bytes.empty()) return -1;
  if (dstSize < static_cast<int>(slot->bytes.size())) return -1;
  std::memcpy(dst, slot->bytes.data(), slot->bytes.size());
  return static_cast<int>(slot->bytes.size());
}

// Load model from LIBSVM text bytes → handle (>0) or 0; bytes cached.
EMSCRIPTEN_KEEPALIVE
int svmLoadModel(const uint8_t* bytes, int size) {
  if (bytes == nullptr || size <= 0) return 0;
  if (!writeFile(kModelPath, bytes, size)) return 0;
  svm_model* model = svm_load_model(kModelPath);
  if (model == nullptr) return 0;
  return StoreModel(model,
                    std::string(reinterpret_cast<const char*>(bytes),
                                static_cast<size_t>(size)));
}

// Predict over columnar data; 0 on success, -1 on failure.
EMSCRIPTEN_KEEPALIVE
int svmPredict(int handle, const float* colMajor, int nRows, int nCols,
               float* out, int outLen) {
  ModelSlot* slot = GetSlot(handle);
  if (slot == nullptr || outLen < nRows || nRows <= 0 || nCols <= 0) return -1;

  std::vector<svm_node> row(nCols + 1);
  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      row[j].index = j + 1;
      row[j].value = colMajor[static_cast<size_t>(j) * nRows + i];
    }
    row[nCols].index = -1;
    out[i] = static_cast<float>(svm_predict(slot->model, row.data()));
  }
  return 0;
}

EMSCRIPTEN_KEEPALIVE
void svmFreeModel(int handle) {
  ModelSlot* slot = GetSlot(handle);
  if (slot == nullptr) return;
  svm_free_and_destroy_model(&slot->model);  // frees and nulls slot->model
  slot->bytes.clear();
  slot->bytes.shrink_to_fit();
}

}  // extern "C"
