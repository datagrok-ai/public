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
  sketcher: Sketcher = new Sketcher();
  bitset: DG.BitSet | null = null;
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
    /* No await */ initRdKitService();
    this.root = ui.divV([]);
    this._indicateProgress(false);
    this.loader.style.position = 'absolute';
    this.loader.style.right = '60px';
    this.loader.style.top = '4px';
    this.sketcher.onChanged.subscribe(async (_) => await this._onSketchChanged());
    this.root.appendChild(this.sketcher.root);
    this.root.firstChild!.firstChild!.firstChild!.appendChild(this.loader);
  }

  attach(dataFrame: DG.DataFrame) {
    this.column = dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    super.attach(dataFrame);
  }

  applyFilter(): void {
    if (this.bitset) {
      this.dataFrame?.filter.and(this.bitset);
      this.column!.temp['chem-scaffold-filter'] = this.sketcher.getMolFile();
    }
  }

  async _onSketchChanged(): Promise<void> {
    if (!this.sketcher.getMolFile() || this.sketcher.getMolFile().endsWith(this.WHITE_MOL)) {
      if (this.column?.temp['chem-scaffold-filter'])
        delete this.column.temp['chem-scaffold-filter'];
      this.dataFrame?.filter.setAll(true);
      this.bitset = null;
    } else {
      this._indicateProgress();
      this.bitset = await chemSubstructureSearchLibrary(
        this.column!, this.sketcher.getMolFile(), await this.sketcher.getSmarts());
      this._indicateProgress(false);
      this.dataFrame?.rows.requestFilter();
    }
  }
}
