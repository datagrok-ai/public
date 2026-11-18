const SVM_TYPES = ['c_svc', 'nu_svc', 'epsilon_svr', 'nu_svr'];
const KERNELS = ['linear', 'polynomial', 'rbf', 'sigmoid'];

export interface LibsvmModel {
  svmType: string,
  kernel: string,
  gamma: number,
  degree: number,
  coef0: number,
  nrClass: number,
  label: number[],
  nrSv: number[],
  /** Offset of each class's first SV */
  starts: number[],
  rho: number[],
  svCoef: Float64Array[],
  /** SVs in CSR form: SV `s` has the 1-based feature indices `svIndex[svOffsets[s]..svOffsets[s + 1])` */
  svOffsets: Int32Array,
  svIndex: Int32Array,
  svValue: Float64Array,
}

function isClassifier(model: LibsvmModel): boolean {
  return model.svmType === 'c_svc' || model.svmType === 'nu_svc';
}

/** A model in the text format of libsvm's `svm_save_model` */
export function parseLibsvmModel(text: string): LibsvmModel {
  const model: LibsvmModel = {svmType: '', kernel: '', gamma: 0, degree: 3, coef0: 0, nrClass: 2, label: [], nrSv: [],
    starts: [], rho: [], svCoef: [], svOffsets: new Int32Array(0), svIndex: new Int32Array(0),
    svValue: new Float64Array(0)};
  const lines = text.split('\n');
  let totalSv = 0;
  let i = 0;
  for (; i < lines.length; ++i) {
    const [key, ...rest] = lines[i].trim().split(/\s+/);
    if (key === 'SV') {
      ++i;
      break;
    }
    const numbers = rest.map(Number);
    switch (key) {
    case 'svm_type': model.svmType = rest[0]; break;
    case 'kernel_type': model.kernel = rest[0]; break;
    case 'gamma': model.gamma = numbers[0]; break;
    case 'degree': model.degree = numbers[0]; break;
    case 'coef0': model.coef0 = numbers[0]; break;
    case 'nr_class': model.nrClass = numbers[0]; break;
    case 'total_sv': totalSv = numbers[0]; break;
    case 'rho': model.rho = numbers; break;
    case 'label': model.label = numbers; break;
    case 'nr_sv': model.nrSv = numbers; break;
    }
  }
  if (!SVM_TYPES.includes(model.svmType))
    throw new Error(`LIBSVM: unsupported svm_type ${model.svmType}`);
  if (!KERNELS.includes(model.kernel))
    throw new Error(`LIBSVM: unsupported kernel_type ${model.kernel}`);

  const k = isClassifier(model) ? model.nrClass - 1 : 1;
  model.svCoef = Array.from({length: k}, () => new Float64Array(totalSv));
  model.svOffsets = new Int32Array(totalSv + 1);
  const index: number[] = [];
  const value: number[] = [];
  for (let s = 0; s < totalSv; ++s, ++i) {
    const parts = lines[i].trim().split(/\s+/);
    for (let c = 0; c < k; ++c)
      model.svCoef[c][s] = Number(parts[c]);
    for (let p = k; p < parts.length; ++p) {
      const [idx, v] = parts[p].split(':');
      index.push(Number(idx));
      value.push(Number(v));
    }
    model.svOffsets[s + 1] = index.length;
  }
  model.svIndex = Int32Array.from(index);
  model.svValue = Float64Array.from(value);

  model.starts = [0];
  for (let c = 1; c < model.nrClass; ++c)
    model.starts[c] = model.starts[c - 1] + (model.nrSv[c - 1] ?? 0);
  return model;
}

function dot(model: LibsvmModel, s: number, x: ArrayLike<number>): number {
  let sum = 0;
  for (let p = model.svOffsets[s]; p < model.svOffsets[s + 1]; ++p)
    sum += model.svValue[p] * x[model.svIndex[p] - 1];
  return sum;
}

/** The kernel values of every SV against the dense row `x`, as libsvm's `Kernel::k_function` */
function kernelValues(model: LibsvmModel, x: ArrayLike<number>): Float64Array {
  const total = model.svOffsets.length - 1;
  const kv = new Float64Array(total);
  switch (model.kernel) {
  case 'linear':
    for (let s = 0; s < total; ++s)
      kv[s] = dot(model, s, x);
    return kv;
  case 'polynomial':
    for (let s = 0; s < total; ++s)
      kv[s] = Math.pow(model.gamma * dot(model, s, x) + model.coef0, model.degree);
    return kv;
  case 'sigmoid':
    for (let s = 0; s < total; ++s)
      kv[s] = Math.tanh(model.gamma * dot(model, s, x) + model.coef0);
    return kv;
  default: {
    const seen = new Uint8Array(x.length);
    for (let s = 0; s < total; ++s) {
      // the SV's own features first, then x² of the ones it lacks (none: the old wasm saved dense SVs)
      seen.fill(0);
      let sum = 0;
      for (let p = model.svOffsets[s]; p < model.svOffsets[s + 1]; ++p) {
        const j = model.svIndex[p] - 1;
        const d = model.svValue[p] - x[j];
        sum += d * d;
        seen[j] = 1;
      }
      for (let j = 0; j < x.length; ++j) {
        if (!seen[j])
          sum += x[j] * x[j];
      }
      kv[s] = Math.exp(-model.gamma * sum);
    }
    return kv;
  }
  }
}

/** libsvm's `svm_predict` of the dense row `x`: the class label (one-vs-one voting) or the regression value */
export function svmPredict(model: LibsvmModel, x: ArrayLike<number>): number {
  const kv = kernelValues(model, x);
  if (!isClassifier(model)) {
    let sum = 0;
    for (let i = 0; i < kv.length; ++i)
      sum += model.svCoef[0][i] * kv[i];
    return sum - model.rho[0];
  }

  const nc = model.nrClass;
  const starts = model.starts;
  const votes = new Int32Array(nc);
  let p = 0;
  for (let i = 0; i < nc; ++i) {
    for (let j = i + 1; j < nc; ++j) {
      let sum = 0;
      const c1 = model.svCoef[j - 1];
      const c2 = model.svCoef[i];
      for (let k = 0; k < model.nrSv[i]; ++k)
        sum += c1[starts[i] + k] * kv[starts[i] + k];
      for (let k = 0; k < model.nrSv[j]; ++k)
        sum += c2[starts[j] + k] * kv[starts[j] + k];
      sum -= model.rho[p++];
      ++votes[sum > 0 ? i : j];
    }
  }
  let best = 0;
  for (let i = 1; i < nc; ++i) {
    if (votes[i] > votes[best])
      best = i;
  }
  return model.label[best];
}
