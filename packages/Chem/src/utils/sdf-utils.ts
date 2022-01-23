import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';


export function _saveAsSdf() {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: open dialog
  //todo: UI for choosing structure column if necessary
  //todo: UI for choosing columns with properties

  const table = grok.shell.t;
  const structureColumn = table.columns.bySemType('Molecule');
  if (structureColumn == null)
    return;

  let result = '';

  for (let i = 0; i < table.rowCount; i++) {
    try {
      const mol = OCL.Molecule.fromSmiles(structureColumn.get(i));
      result += `\n${mol.toMolfile()}\n`;

      // properties
      for (const col of table.columns) {
        if (col !== structureColumn)
          result += `>  <${col.name}>\n${col.get(i)}\n\n`;
      }

      result += '$$$$';
    } catch (error) {
      console.error(error);
    }
  }

  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}


class SDFReader {
  dataColls: { [_: string]: any };

  constructor() {
    this.dataColls = {'molecule': []};
  }

  get_colls(content: string) {
    this.read(content);
    return this.dataColls;
  }

  read(content: string) {
    let startIndex = content.indexOf('$$$$', 0);
    this.parse(content, 0, startIndex, (name: string, val: any) => { // TODO: type
      this.dataColls[name] = [];
      this.dataColls[name].push(val);
    });
    startIndex += 5;
    while (startIndex > -1 && startIndex < content.length)
      startIndex = this.readNext(content, startIndex);
  }

  readNext(content: string, startIndex: number) {
    const nextStartIndex = content.indexOf('$$$$', startIndex);
    if (nextStartIndex === -1)
      return -1;
    else {
      this.parse(content, startIndex, nextStartIndex,
        (name: string, val: number) => this.dataColls[name].push(val));
    }

    if (nextStartIndex > -1)
      return nextStartIndex + 5;

    return nextStartIndex;
  }

  parse(content: string, start: number, end: number, handler: any) {
    const molEnd = +content.indexOf('M  END\n', start) + 7;
    let localEnd = start;
    this.dataColls['molecule'].push(content.substr(start, molEnd - start));

    start = molEnd;
    while (localEnd < end) {
      start = content.indexOf('> <', localEnd);
      if (start === -1)
        return;

      start += 3;
      localEnd = content.indexOf('>\n', start);
      if (localEnd === -1)
        return;

      const propertyName = content.substr(start, localEnd - start);
      start = localEnd + 2;

      localEnd = content.indexOf('\n', start);
      if (localEnd === -1)
        localEnd = end;

      handler(propertyName, content.substr(start, localEnd - start));
      localEnd += 2;
    }
  }
}
