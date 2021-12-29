/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {chem} from "datagrok-api/grok";
import Sketcher = chem.Sketcher;

export class SubstructureFilter extends DG.Filter {

  column: DG.Column | null = null;
  sketcher: Sketcher = new Sketcher();
  readonly WHITE_MOL = '\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n';

  get filterSummary(): string {
    return this.sketcher.getSmiles();
  }

  get isFiltering(): boolean {
    return this.sketcher.getSmiles() !== '';
  }

  constructor() {
    super();
    this.root = ui.divV([]);
    this.sketcher.onChanged.subscribe((_) => this.dataFrame?.rows.requestFilter());
    this.root.appendChild(this.sketcher.root);
  }

  attach(dataFrame: DG.DataFrame) {
    this.dataFrame = dataFrame;
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
    if (!this.sketcher.getMolFile() || this.sketcher.getMolFile().endsWith(this.WHITE_MOL)) {
      this.reset();
      return;
    }
    grok.functions.call('Chem:searchSubstructure', {
      'molStringsColumn' : this.column,
      'molString': this.sketcher.getMolFile(),
      'substructLibrary': true,
      'molStringSmarts': ''
    })
    .then((bitset_col) => {
      this.dataFrame?.filter.and(bitset_col.get(0));
      this.column!.temp['chem-scaffold-filter'] = this.sketcher.getMolFile(); // not sure if !
      this.dataFrame?.filter.fireChanged();
    }).catch((e) => {
      console.warn(e);
      this.reset();
    })
  }
}