const m = 3;
const n = 2;
const c = 4;

const X = new Array<Float32Array>(m);
X[0] = new Float32Array([1, 2]);
X[1] = new Float32Array([3, 4]);
X[2] = new Float32Array([5, 6]);

const W = new Array<Float32Array>(c);
W[0] = new Float32Array([4, 5, 6]);
W[1] = new Float32Array([-4, -5, -6]);
W[2] = new Float32Array([1, -1, 2]);
W[3] = new Float32Array([2, -1, 3]);

const Z = new Array<Float32Array>(m);
for (let i = 0; i < m; ++i)
  X[i] = new Float32Array(c);

const S = new Float32Array(m);
