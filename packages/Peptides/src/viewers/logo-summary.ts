import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {PeptidesModel} from '../model';

export class LogoSummary extends DG.JsViewer {
  _titleHost = ui.divText('Logo Summary Table', {id: 'pep-viewer-title'});
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

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  render(): void {
    if (this.initialized) {
      $(this.root).empty();
      this.viewerGrid.root.style.width = 'auto';
      this.root.appendChild(ui.divV([this._titleHost, this.viewerGrid.root]));
      this.viewerGrid.invalidate();
    }
  }
}
