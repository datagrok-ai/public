/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class SubstructureFilter extends DG.Filter {

  molfile: string = '';
  column: DG.Column | null = null;
  _sketcher: HTMLElement;
  readonly WHITE_MOL = `
  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`;

  constructor() {
    super();
    this.root = ui.divV([]);
    this._sketcher = grok.chem.sketcher((_: any, molfile: string) => {
      this.molfile = molfile;
      this.dataFrame?.rows.requestFilter();
    });
    this.root.appendChild(this._sketcher);
  }

  attach(dFrame: DG.DataFrame) {
    this.dataFrame = dFrame;
    this.column = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.subs.push(this.dataFrame.onRowsFiltering.subscribe((_: any) => { this.applyFilter(); } ));
  }

  // TODO: this needs to be triggered
  reset() {
    // this.dataFrame.filter.setAll(true, false);
    if (this.column?.temp['chem-scaffold-filter'])
      delete this.column.temp['chem-scaffold-filter'];
    // this._sketcher?.setSmiles('');
    // this.dataFrame.filter.fireChanged();
  }

  applyFilter() {
    if (!this.molfile || this.molfile.endsWith(this.WHITE_MOL)) {
      this.reset();
      return;
    }
    grok.functions.call('Chem:searchSubstructure', {
      'molStringsColumn' : this.column,
      'molString': this.molfile,
      'substructLibrary': true,
      'molStringSmarts': ''
    })
    .then((bitset_col) => {
      this.dataFrame?.filter.and(bitset_col.get(0));
      this.column!.temp['chem-scaffold-filter'] = this.molfile; // not sure if !
      this.dataFrame?.filter.fireChanged();
    }).catch((e) => {
      console.warn(e);
      this.reset();
    })
  }
  
  detach() {
    this.subs.forEach((s) => s.unsubscribe());
  }
}