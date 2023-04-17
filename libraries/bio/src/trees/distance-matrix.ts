import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Distance matrix class compatible with data structure of scipy.spatial.distance.pdist */
export class DistanceMatrix {
  _data: Float32Array;
  _size: number;

  get data(): Float32Array { return this._data; }

  get size(): number { return this._size; }

  /**
   * @param {Float64Array} data  Distance data
   * @param {number} m           Number of original observations
   */
  constructor(data?: Float32Array, size?: number) {
    if (size == undefined) {
      if (data == undefined) throw new Error('Arguments error: data or size is required.');

      this._data = data!;
      this._size = (1 + Math.sqrt(1 + 4 * 2 * this._data.length)) / 2;
      if (this._size != Math.floor(this._size))
        throw new Error(`Invalid data length ${this._data.length} leads to non integer size ${this._size}`);
    } else {
      this._size = size;
      const dataLength: number = size * (size - 1) / 2;
      if (data) {
        if (data.length != dataLength)
          throw new Error(`Invalid data length. Observations size ${size} requires data length ${dataLength}.`);
        this._data = data;
      } else {
        this._data = new Float32Array(dataLength);
      }
    }
  }

  private _linearizeIJ(i: number, j: number): number {
    if (!(i < j)) throw new Error('i must be less than j');
    return this._size * i + j - Math.floor(((i + 2) * (i + 1)) / 2);
  }

  get(i: number, j: number) {
    if (i == j)
      return 0;
    else if (i < j)
      return this._data[this._linearizeIJ(i, j)];
    else
      return this._data[this._linearizeIJ(j, i)];
  }

  set(i: number, j: number, value: number) {
    this._data[this._linearizeIJ(i, j)] = value;
  }

  static calc<TObj>(list: Indexable<TObj>, method: (a: TObj, b: TObj) => number): DistanceMatrix {
    const size: number = list.length;
    const res = new DistanceMatrix(undefined, size);
    for (let i = 0; i < size; i++) {
      for (let j = i + 1; j < size; j++)
        res.set(i, j, method(list[i], list[j]));
    }
    return res;
  }

  // Combines different matrices into one using simple vector math.
  // Used to combine distance matrices for different features like ones of numbers and ones of sequences.
  static combine(matrices: DistanceMatrix[], method: 'euclidean' | 'manhattan' = 'euclidean') {
    matrices.forEach((mat) => mat.normileze());
    const size = matrices[0].size;
    const res = new DistanceMatrix(undefined, size);
    for (let i = 0; i < size; i++) {
      for (let j = i + 1; j < size; j++) {
        let sum = 0;
        switch (method) {
        case 'manhattan':
          for (let k = 0; k < matrices.length; k++)
            sum += matrices[k].get(i, j);

          res.set(i, j, sum);
          break;
        default:
          for (let k = 0; k < matrices.length; k++)
            sum += matrices[k].get(i, j) ** 2;

          res.set(i, j, Math.sqrt(sum));
          break;
        }
      }
    }
    return res;
  }

  //normilze distance matrix in place
  normileze() {
    const max = Math.max(...this.data);
    const min = Math.min(...this.data);
    const range = max - min;
    for (let i = 0; i < this.data.length; i++)
      this.data[i] = range === 0 ? this.data[i] - min : (this.data[i] - min) / (max - min);
  }
}

export interface Indexable<TObj> {
  [index: number]: TObj;
  get length(): number;
}
