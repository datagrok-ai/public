/*const m = 3;
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

// Predict
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

//console.log(Z);
//console.log('=====================');
//console.log(predClass);

const dZ = new Array<Float32Array>(c);
for (let i = 0; i < c; ++i)
  dZ[i] = new Float32Array(m);

const Y = new Array<Uint8Array>(m);
Y[0] = new Uint8Array([1, 0, 0, 0]);
Y[1] = new Uint8Array([0, 1, 0, 0]);
Y[2] = new Uint8Array([0, 0, 1, 0]);
let yTrue: Uint8Array;
let yPred: Float32Array;

//console.log(Y);
//console.log('=====================');

// dZ
for (let j = 0; j < m; ++j) {
  yPred = Z[j];
  yTrue = Y[j];

  for (let i = 0; i < c; ++i)
    dZ[i][j] = yPred[i] - yTrue[i];
}

//console.log(dZ);
//console.log('=====================');

const dW = new Array<Float32Array>(c);
for (let i = 0; i < c; ++i)
  dW[i] = new Float32Array(n + 1);
let dWbuf: Float32Array;

// dB
for (let i = 0; i < c; ++i) {
  sum = 0;
  zBuf = dZ[i];

  for (let j = 0; j < m; ++j)
    sum += zBuf[j];

  dW[i][n] = sum / m;
}

// x transposed
const transposedX = new Array<Float32Array>(n);
transposedX[0] = new Float32Array([1, 3, 5]);
transposedX[1] = new Float32Array([2, 4, 6]);

// dW
for (let i = 0; i < c; ++i) {
  zBuf = dZ[i];
  wBuf = dW[i];

  for (let j = 0; j < n; ++j) {
    xBuf = transposedX[j];

    sum = 0;
    for (let k = 0; k < m; ++k)
      sum += zBuf[k] * xBuf[k];

    wBuf[j] = sum / m;
  }
}

const alpha = 0.01;
const lambda = 0.1;

console.log(W);
console.log('--------------------------');
console.log(dW);
console.log('==========================');

// update W
for (let i = 0; i < c; ++i) {
  wBuf = W[i];
  dWbuf = dW[i];

  for (let j = 0; j < n; ++j)
    wBuf[j] = (1 - lambda * alpha / m) * wBuf[j] - alpha * dWbuf[j];

  wBuf[n] -= alpha * dWbuf[n];
}

console.log(W);
*/
