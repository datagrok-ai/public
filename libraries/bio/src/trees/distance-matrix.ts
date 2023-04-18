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
   * @param {Float32Array} data Distance data
   * @param {number} size Number of original observations
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

  // squares each value in matrix in place
  public square() {
    for (let i = 0; i < this._data.length; i++)
      this._data[i] = this._data[i] ** 2;
  }

  // adds another matrix to this one in place
  public add(other: DistanceMatrix) {
    if (this._size != other._size)
      throw new Error(`Matrices must have the same size. This size: ${this._size}, other size: ${other._size}`);
    for (let i = 0; i < this._data.length; i++)
      this._data[i] += other._data[i];
  }

  // square root each value in matrix in place
  public sqrt() {
    for (let i = 0; i < this._data.length; i++)
      this._data[i] = Math.sqrt(this._data[i]);
  }

  //normilze distance matrix in place
  public normalize() {
    const max = Math.max(...this._data);
    const min = Math.min(...this._data);
    const range = max - min;
    for (let i = 0; i < this._data.length; i++)
      this._data[i] = range === 0 ? this._data[i] - min : (this._data[i] - min) / (max - min);
  }
}

export interface Indexable<TObj> {
  [index: number]: TObj;
  get length(): number;
}
