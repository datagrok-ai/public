import {dmLinearIndex} from '../distance-matrix';


export class TSNE {
  private perplexity: number;
  private dim: number;
  private epsilon: number;
  private iter = 0;
  private N!: number;
  private P!: Float32Array;
  private Y!: any[];
  private gains!: any[];
  private ystep!: any[];
  private random: () => number = Math.random;
  private assert(condition: any, message: any) {
    if (!condition) throw message || 'Assertion failed';
  }

  // syntax sugar
  private getopt(opt: {[_: string]: any}, field: string, defaultval: any) {
    if (opt.hasOwnProperty(field))
      return opt[field];
    else
      return defaultval;
  }

    // return 0 mean unit standard deviation random number
    private returnV = false;
    private vValue = 0.0;
    gaussRandom(): number {
      if (this.returnV) {
        this.returnV = false;
        return this.vValue;
      }
      const u = 2 * this.random() - 1;
      const v = 2 * this.random() - 1;
      const r = u * u + v * v;
      if (r === 0 || r > 1) return this.gaussRandom();
      const c = Math.sqrt(-2 * Math.log(r) / r);
      this.vValue = v * c; // cache this for next function call for efficiency
      this.returnV = true;
      return u * c;
    }

    // return random normal number
    private randn(mu: number, std: number) {
      return mu + this.gaussRandom() * std;
    }

    // utilitity that creates contiguous vector of zeros of size n
    private zeros(n: number) {
      if (typeof (n) === 'undefined' || isNaN(n)) return new Float32Array();
      if (typeof ArrayBuffer === 'undefined') {
        // lacking browser support
        const arr = new Float32Array(n);
        for (let i = 0; i < n; i++) arr[i] = 0;
        return arr;
      } else {
        return new Float32Array(n); // typed arrays are faster
      }
    }

    // utility that returns 2d array filled with random numbers
    // or with value s, if provided
    private randn2d(n: number, d: number, s?: number) {
      const uses = typeof s !== 'undefined';
      const x = (new Array(n)).fill(null).map(() => new Float32Array(d));
      if (uses) {
        for (let i = 0; i < n; i++)
          x[i] = new Float32Array(d).fill(s);
      } else {
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < d; j++)
            x[i][j] = this.randn(0.0, 1e-4);
        }
      }
      return x;
    }

    // compute L2 distance between two vectors
    private l2(x1: number[], x2: number[]) {
      const D = x1.length;
      let d = 0;
      for (let i = 0; i < D; i++) {
        const x1i = x1[i];
        const x2i = x2[i];
        d += (x1i - x2i) * (x1i - x2i);
      }
      return d;
    }

    // compute pairwise distance in all vectors in X
    private xtod(X: Array< Array<number> >) {
      const N = X.length;
      const dist = this.zeros(N * N); // allocate contiguous array
      for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
          const d = this.l2(X[i], X[j]);
          dist[i * N + j] = d;
          dist[j * N + i] = d;
        }
      }
      return dist;
    }

    public getIterSize(rowSize: number) {
      if (rowSize <= 2000)
        return 100;
      else if (rowSize <= 3000)
        return 90;
      else if (rowSize <= 5000)
        return 80;
      else
        return 70;
    }

    // compute (p_{i|j} + p_{j|i})/(2n)
    private d2p(D: Float32Array, perplexity: number, tol: number, rowSize: number) {
      // const nf = Math.sqrt(D.length); // this better be an integer
      // const n = Math.floor(nf);
      // this.assert(n === nf, 'D should have square number of elements.');
      const indexer = dmLinearIndex(rowSize);
      const distances = new Float32Array(rowSize * rowSize);
      for (let i = 0; i < rowSize; i++) {
        for (let j = i + 1; j < rowSize; j++)
          distances[i * rowSize + j] = distances[j * rowSize + i] = D[indexer(i, j)];
      }
      const n = rowSize;
      const hTarget = Math.log(perplexity); // target entropy of distribution
      const P = this.zeros(n * n); // temporary probability matrix
      const prow = this.zeros(n); // a temporary storage compartment
      for (let i = 0; i < n; i++) {
        let betamin = -Infinity;
        let betamax = Infinity;
        let beta = 1; // initial value of precision
        let done = false;
        const maxtries = Math.floor(this.getIterSize(rowSize) / 5);

        // perform binary search to find a suitable precision beta
        // so that the entropy of the distribution is appropriate
        let num = 0;
        while (!done) {
          //debugger;

          // compute entropy and kernel row with beta precision
          let psum = 0.0;
          for (let j = 0; j < n; j++) {
            const pj = i === j ? 0 : Math.exp(- distances[i * n + j] * beta);
            prow[j] = pj;
            psum += pj;
          }
          // normalize p and compute entropy
          let nHere = 0.0;
          for (let j = 0; j < n; j++) {
            let pj;
            if (psum === 0)
              pj = 0;
            else
              pj = prow[j] / psum;

            prow[j] = pj;
            if (pj > 1e-7) nHere -= pj * Math.log(pj);
          }

          // adjust beta based on result
          if (nHere > hTarget) {
            // entropy was too high (distribution too diffuse)
            // so we need to increase the precision for more peaky distribution
            betamin = beta; // move up the bounds
            if (betamax === Infinity) beta = beta * 2; else beta = (beta + betamax) / 2;
          } else {
            // converse case. make distrubtion less peaky
            betamax = beta;
            if (betamin === -Infinity) beta = beta / 2; else beta = (beta + betamin) / 2;
          }

          // stopping conditions: too many tries or got a good precision
          num++;
          if (Math.abs(nHere - hTarget) < tol) done = true;
          if (num >= maxtries) done = true;
        }

        // console.log('data point ' + i + ' gets precision ' + beta + ' after ' + num + ' binary search steps.');
        // copy over the final prow to P at row i
        for (let j = 0; j < n; j++) P[i * n + j] = prow[j];
      } // end loop over examples i

      // symmetrize P and normalize it to sum to 1 over all ij
      const pOut = this.zeros(n * n);
      const N2 = n * 2;
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++)
          pOut[i * n + j] = Math.max((P[i * n + j] + P[j * n + i]) / N2, 1e-100);
      }

      return pOut;
    }

    // helper function
    private sign(x: number) { return x > 0 ? 1 : x < 0 ? -1 : 0; }

    constructor(opt: {[_: string]: any}) {
      opt = opt || {};
      this.perplexity = this.getopt(opt, 'perplexity', 30); // effective number of nearest neighbors
      this.dim = this.getopt(opt, 'dim', 2); // by default 2-D tSNE
      this.epsilon = this.getopt(opt, 'epsilon', 10); // learning rate
      this.random = this.getopt(opt, 'random', Math.random); // random number generator
    }

    // this function takes a set of high-dimensional points
    // and creates matrix P from them using gaussian kernel
    public initDataRaw(X: Array<Array<any>>) {
      const N = X.length;
      const D = X[0].length;
      this.assert(N > 0, ' X is empty? You must have some data!');
      this.assert(D > 0, ' X[0] is empty? Where is the data?');
      const dists = this.xtod(X); // convert X to distances using gaussian kernel
      this.P = this.d2p(dists, this.perplexity, 1e-4, D); // attach to object
      this.N = N; // back up the size of the dataset
      this.initSolution(); // refresh this
    }

    // this function takes a given distance matrix and creates
    // matrix P from them.
    // D is assumed to be provided as a list of lists, and should be symmetric
    public initDataDist(D: Float32Array, rowSize: number) {
      const N = D.length;
      this.assert(N > 0, ' X is empty? You must have some data!');
      // convert D to a (fast) typed array version
      // dists.forEach((d) => console.log(d));
      console.time('distances to matrix');
      this.P = this.d2p(D, this.perplexity, 1e-4, rowSize);
      console.timeEnd('distances to matrix');
      this.N = rowSize;
      this.initSolution(); // refresh this
    }

    // (re)initializes the solution to random
    public initSolution() {
      // generate random solution to t-SNE
      this.Y = this.randn2d(this.N, this.dim); // the solution
      this.gains = this.randn2d(this.N, this.dim, 1.0); // step gains to accelerate progress in unchanging directions
      this.ystep = this.randn2d(this.N, this.dim, 0.0); // momentum accumulator
      this.iter = 0;
    }

    // return pointer to current solution
    public getSolution() {
      return this.Y;
    }

    // perform a single step of optimization to improve the embedding
    public step() {
      this.iter += 1;
      const N = this.N;
      const cg = this.costGrad(this.Y); // evaluate gradient
      const cost = cg.cost;
      const grad = cg.grad;

      // perform gradient step
      const ymean = this.zeros(this.dim);
      for (let i = 0; i < N; i++) {
        for (let d = 0; d < this.dim; d++) {
          const gid = grad[i][d];
          const sid = this.ystep[i][d];
          const gainid = this.gains[i][d];

          // compute gain update
          let newgain = this.sign(gid) === this.sign(sid) ? gainid * 0.8 : gainid + 0.2;
          if (newgain < 0.01) newgain = 0.01; // clamp
          this.gains[i][d] = newgain; // store for next turn

          // compute momentum step direction
          const momval = this.iter < 250 ? 0.5 : 0.8;
          const newsid = momval * sid - this.epsilon * newgain * grad[i][d];
          this.ystep[i][d] = newsid; // remember the step we took

          // step!
          this.Y[i][d] += newsid;

          ymean[d] += this.Y[i][d]; // accumulate mean so that we can center later
        }
      }

      // reproject Y to be zero mean
      for (let i = 0; i < N; i++) {
        for (let d = 0; d < this.dim; d++)
          this.Y[i][d] -= ymean[d] / N;
      }

      //if(this.iter%100===0) console.log('iter ' + this.iter + ', cost: ' + cost);
      return cost; // return current cost
    }

    // for debugging: gradient check
    public debugGrad() {
      const N = this.N;

      const cg = this.costGrad(this.Y); // evaluate gradient
      //const cost = cg.cost;
      const grad = cg.grad;

      const e = 1e-5;
      for (let i = 0; i < N; i++) {
        for (let d = 0; d < this.dim; d++) {
          const yold = this.Y[i][d];

          this.Y[i][d] = yold + e;
          const cg0 = this.costGrad(this.Y);

          this.Y[i][d] = yold - e;
          const cg1 = this.costGrad(this.Y);

          const analytic = grad[i][d];
          const numerical = (cg0.cost - cg1.cost) / (2 * e);
          console.log(i + ',' + d + ': gradcheck analytic: ' + analytic + ' vs. numerical: ' + numerical);

          this.Y[i][d] = yold;
        }
      }
    }

    private quArr?: Float32Array;
    // return cost and gradient, given an arrangement
    public costGrad(Y: number[][]) {
      const N = this.N;
      const dim = this.dim; // dim of output space
      const P = this.P;

      const pmul = this.iter < 100 ? 4 : 1; // trick that helps with local optima

      // compute current Q distribution, unnormalized first
      this.quArr ??= this.zeros(N * N);
      let qsum = 0.0;
      for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
          const dsum = new Array(dim).reduce((prev, _, d) => prev + Math.pow(Y[i][d] - Y[j][d], 2), 0);
          const qu = 1.0 / (1.0 + dsum); // Student t-distribution
          this.quArr[i * N + j] = qu;
          this.quArr[j * N + i] = qu;
          qsum += 2 * qu;
        }
      }
      // normalize Q distribution to sum to 1
      //const NN = N * N;
      //const Q = this.zeros(NN);
      let cost = 0.0;
      const grad = new Array(N).fill(null).map(() => (new Float32Array(dim)).fill(0.0));
      for (let i = 0; i < N; i++) {
        // const gsum = new Float32Array(dim).fill(0.0); // init grad for point i
        for (let j = i + 1; j < N; j++) {
          const Q = Math.max(this.quArr[i * N + j] / qsum, 1e-100);
          cost += (- P[i * N + j] * Math.log(Q)) * 2; // accumulate cost (the non-constant portion at least...)
          const premult = 4 * (pmul * P[i * N + j] - Q) * this.quArr[i * N + j];
          //console.log(premult, Q, i, j);
          for (let d = 0; d < dim; d++) {
            grad[i][d] += premult * (Y[i][d] - Y[j][d]);
            grad[j][d] += premult * (Y[j][d] - Y[i][d]);
          }
        }
        //grad[i] = gsum;
      }

      return {cost, grad};
    }
}
