/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {chem} from 'datagrok-api/grok';
import {chemSubstructureSearchLibrary} from "./chem-searches";
import {initRdKitService} from "./chem-common-rdkit";
import Sketcher = chem.Sketcher;

export class SubstructureFilter extends DG.Filter {
  column: DG.Column | null = null;
  sketcher: Sketcher = new Sketcher();
  loader = ui.loader();
  readonly WHITE_MOL = '\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n';

  _indicateProgress(on = true) {
    this.loader.style.display = on ? 'block' : 'none';
  }

  get filterSummary(): string {
    return this.sketcher.getSmiles();
  }

  get isFiltering(): boolean {
    return this.sketcher.getSmiles() !== '';
  }

  constructor() {
    super();
    /* No await! */ initRdKitService();
    this.root = ui.divV([]);
    this._indicateProgress(false);
    this.loader.style.position = 'absolute';
    this.loader.style.right = '60px';
    this.loader.style.top = '4px';
    this.sketcher.onChanged.subscribe((_) => this.dataFrame?.rows.requestFilter());
    this.root.appendChild(this.sketcher.root);
    this.root!.firstChild!.firstChild!.firstChild!.appendChild(this.loader);
  }

  attach(dataFrame: DG.DataFrame) {
    this.dataFrame = dataFrame;
    this.column = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    /*
    const debounceMsec = 500;
    this.subs.push(DG.debounce(this.dataFrame.onRowsFiltering, debounceMsec)
      .subscribe((_: any) => this.applyFilter()));
    */
    this.subs.push(this.dataFrame.onRowsFiltering.subscribe((_: any) => this.applyFilter()));
  }

  // TODO: this needs to be triggered
  reset() {
    if (this.column?.temp['chem-scaffold-filter'])
      delete this.column.temp['chem-scaffold-filter'];
  }

  applyFilter() {
    if (!this.sketcher.getMolFile() || this.sketcher.getMolFile().endsWith(this.WHITE_MOL)) {
      this.reset();
      return;
    }

    this._indicateProgress();
    (async() => {
      await chemSubstructureSearchLibrary(this.column!, this.sketcher.getMolFile(), await this.sketcher.getSmarts())
        .then((bitset) => {
          // console.log('Substructure filtering produced results');
          this.dataFrame?.filter.and(bitset);
          this.column!.temp['chem-scaffold-filter'] = this.sketcher.getMolFile();
          this.dataFrame?.filter.fireChanged();
          this._indicateProgress(false);
        }).catch((e) => {
        console.warn(e);
        this.reset();
      });
    })();
  }
}
