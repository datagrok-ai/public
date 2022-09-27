import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import {PeptidesModel} from '../model';

export class LogoSummary extends DG.JsViewer {
  model!: PeptidesModel;
  viewerGrid!: DG.Grid;
  initialized: boolean = false;

  constructor() {
    super();
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();

    this.model = await PeptidesModel.getInstance(this.dataFrame);
    
    this.subs.push(this.model.onLogoSummaryGridChanged.subscribe((grid) => {
      this.viewerGrid = grid;
      this.render();
    }));
    this.model.updateDefault();
    this.viewerGrid = this.model.logoSummaryGrid;
    this.initialized = true;
    this.render();
  }

  detach(): void {this.subs.forEach(sub => sub.unsubscribe());}

  render(): void {
    if (this.initialized) {
      $(this.root).empty();
      this.root.appendChild(this.viewerGrid.root);
      this.viewerGrid.invalidate();
    }
  }
}
