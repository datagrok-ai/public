const m = 3;
const n = 2;
const c = 4;

const X = new Array<Float32Array>(m);
X[0] = new Float32Array([1, 2]);
X[1] = new Float32Array([3, 4]);
X[2] = new Float32Array([5, 6]);
let xBuf: Float32Array;

const W = new Array<Float32Array>(c);
W[0] = new Float32Array([4, -5, 6]);
W[1] = new Float32Array([-4, 5, -6]);
W[2] = new Float32Array([1, 1, -2]);
W[3] = new Float32Array([2, 1, -3]);
let wBuf: Float32Array;

const Z = new Array<Float32Array>(m);
for (let i = 0; i < m; ++i)
  Z[i] = new Float32Array(c);
let zBuf: Float32Array;
let sum: number;
let sumExp: number;
let maxProb: number;
let maxProbIdx: number;

for (let j = 0; j < m; ++j) {
  xBuf = X[j];
  zBuf = Z[j];
  sum = 0;
  sumExp = 0;

  for (let i = 0; i < c; ++i) {
    wBuf = W[i];
    sum = wBuf[n];

    for (let k = 0; k < n; ++k)
      sum += wBuf[k] * xBuf[k];

    zBuf[i] = Math.exp(sum);
    sumExp += zBuf[i];
  }

  for (let i = 0; i < c; ++i)
    zBuf[i] /= sumExp;
}

const predClass = new Uint32Array(m);

for (let i = 0; i < m; ++i) {
  zBuf = Z[i];
  maxProb = zBuf[0];
  maxProbIdx = 0;

  for (let j = 1; j < c; ++j) {
    if (maxProb < zBuf[j]) {
      maxProb = zBuf[j];
      maxProbIdx = j;
    }
  }

  predClass[i] = maxProbIdx;
}


console.log(Z);
console.log('=====================');
console.log(predClass);
