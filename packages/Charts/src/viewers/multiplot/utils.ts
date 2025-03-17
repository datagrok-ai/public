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

  // input: 2d array with table, column number to search categories
  // result: object with unique values as keys and array of rows as values
  getCategories(ar: any, columnNumber: number) : any {
    const r: {[key: string]: any} = {};
    for (let i=0; i<ar.length; i++) {
      const key = ar[i][columnNumber];
      if (!r[key]) r[key] = [];
      r[key].push(ar[i]);
    }
    return r;
  }

  // input: object with categories
  // output: array of categories sorted by number of rows in each category (descending)
  // ['categ1', 'categ2', ...]
  getCategoriesArray(categObject: any) : any {
    const keys = Object.keys(categObject);
    keys.sort((a, b) => categObject[a].length < categObject[b].length ? 1 : -1);
    return keys;
  }

  // convert one multiline describtion to the several plots
  // only objects with field "splitByColumnName"
  generatePlotsFromDescriptions(descr: any, categs: string[], allCats: string[]): any {
    const r = [];
    const categsMax = descr.maxLimit;
    const fieldsNumber = categs.length < categsMax ? categs.length : categsMax;
    for (let i = 0; i < fieldsNumber; i++) {
      const t: any = {};
      t.series = { type: descr.type, data: [], descr: 222 };
      t.x = descr.x;
      t.y = descr.y;
      t.condition = {
        field: descr.splitByColumnName,
        value: categs[i],
      };
      t.title = descr.title ?? categs[i];
      t.height = descr.height;
      t.show = descr.show;
      t.currentCat = categs[i];
      t.allCats = allCats;
      t.tableName = descr.tableName;
      t.yType = 'value';
      t.comboEdit = descr.comboEdit;
      t.multiEdit = descr.multiEdit;
      t.extraFields = descr.extraFields;
      t.multiLineFieldIndex = descr.multiLineFieldIndex;
      t.showLegend = descr.showLegend;
      r.push(t);
    }
    return r;
  }

  concatArrayUnique(ar1: string[], ar2: string[], limit: number) : string[] {
    const r: any[] = [];
    ar1.map((e, i) => (i<limit) ? r.push(e) : '');
    for (let i=0; i<ar2.length; i++) {
      if (!r.includes(ar2[i])) {
        r.push(ar2[i]);
        if (r.length >= limit) break;
      }
    }
    return r;
  }

  getPlotsFromParams(tables: any, plts : any) : any {
    let r: any[] = [];
    for (let i=0; i< plts.length; i++) {
      const p = plts[i];
      if (p.splitByColumnName) {
        const table = tables[p.tableName];
        const data = this.getUniversalData(table, [p.x, p.y, p.splitByColumnName]);
        const cats = this.getCategories(data, 2);
        const sortedCats = this.getCategoriesArray(cats);
        p.allCats = sortedCats;
        const activeCats = this.concatArrayUnique(p.categories ?? [], sortedCats, p.maxLimit);
        p.activeCats = activeCats;
        const plotsArray = this.generatePlotsFromDescriptions(p, activeCats, sortedCats);
        r = r.concat(plotsArray);
      } else {
        if (p.statusChart) {
          if (p.statusChart.splitByColumnName) {
            // generate array for the filter in one plot
            const table = tables[p.tableName];
            const fields = [p.x, p.y].concat(p.extraFields ?? []);
            const data = this.getUniversalData(table, fields);
            const cats = this.getCategories(data, 1);
            const sortedCats = this.getCategoriesArray(cats);
            const allCats = this.concatArrayUnique(p.statusChart.categories ?? [], sortedCats, p.statusChart.maxLimit);
            p.condition = {
              field: fields[1],
              value: allCats,
            };
            p.statusChart.value = allCats;
          }
        }
        r.push(p);
      }
    }
    return r;
  }

  normalize100(column: DG.Column): DG.Row[] {
    const r: any[] = [];
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
  getUniversalData(table: DG.DataFrame, fieldsNames: string[], indexes?: Int32Array, condition? : any): any[] {
    const r = [];
    function getRowFields(row: DG.Row): any[] {
      const fields = [];
      for (let i = 0; i < fieldsNames.length; i++) {
        let cell = row[fieldsNames[i]];
        if (typeof cell === 'number' && cell === DG.INT_NULL)
          cell = 0;

        fields.push(cell);
      }
      return fields;
    };

    if (indexes) {
      for (let ind = 0; ind < indexes.length; ind++) {
        const row = table.row(indexes[ind]);
        const fields = getRowFields(row);
        if (!condition || row[condition.field] === condition.value)
          r.push(fields);

        if (!condition || condition.value.includes(row[condition.field]))
          r.push(fields);
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
    if (str.length <= n)
      return str;

    const subString = str.substr(0, n-1); // the original check
    return (useWordBoundary ?
      subString.substr(0, subString.lastIndexOf(' ')) + '...' :
      subString) + '...';
  }

  getArrayOfColumnNames(colNames: any) {
    return Array.isArray(colNames) ? colNames : [colNames];
  }

  splitToMultipleSeries(colIndex: number) {
  }
}
