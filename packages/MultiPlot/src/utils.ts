import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export class MPUtils {
  constructor() {}

  toStringColor(num : number) : string {
    num >>>= 0;
    const b = num & 0xFF;
    const g = (num & 0xFF00) >>> 8;
    const r = (num & 0xFF0000) >>> 16;
    // a = ( (num & 0xFF000000) >>> 24 ) / 255 ;
    const a = 1;
    return 'rgba(' + [r, g, b, a].join(',') + ')';
  }

  normalize100(column: DG.Column): DG.Row[] {
    const r = [];
    const delta = column.max - column.min;
    const rawData = column.getRawData();
    for (let i = 0; i < column.length; i++) {
      // avoid case when max == min
      const t = delta ? (rawData[i] - column.min) * 100 / delta : rawData[i];
      r.push(t);
    }
    return r;
  }

  // build 2d array for series.data of echart
  getUniversalData(table: DG.DataFrame, fieldsNames: string[], indexes: Int32Array, condition : any): any[] {
    const r = [];

    function getRowFields(row: DG.Row): any[] {
      const fields = [];
      for (let i = 0; i < fieldsNames.length; i++) {
        let cell = row[fieldsNames[i]];
        if (typeof cell === 'number' && cell === DG.INT_NULL) {
          cell = 0;
        }
        fields.push(cell);
      }
      return fields;
    };

    if (indexes) {
      for (let ind = 0; ind < indexes.length; ind++) {
        const row = table.row(indexes[ind]);
        const fields = getRowFields(row);
        if (!condition || row[condition.field] === condition.value) {
          r.push(fields);
        }
      }
    } else {
      for (let i = 0; i < table.rowCount; i++) {
        const row = table.row(i);
        const fields = getRowFields(row);
        r.push(fields);
      }
    }
    return r;
  }

  getBitByIndex32(b: any, index: number) {
    const b2 = b.getBuffer();
    const rez = !!(b2[~~ (index / 32)] & (1 << (index & 31)));
    return rez;
  }

  trimCategoryString(s: string, length : number) : string {
    return s.length > length ? s.substring(0, length) + '...' : s;
  }

  // trim to keep only entire words not used right now
  trimCategoryWord( str: string, n: number, useWordBoundary : boolean) : string {
    if (str.length <= n) {
      return str;
    }
    const subString = str.substr(0, n-1); // the original check
    return (useWordBoundary ?
      subString.substr(0, subString.lastIndexOf(' ')) + '...' :
      subString) + '...';
  }
}
