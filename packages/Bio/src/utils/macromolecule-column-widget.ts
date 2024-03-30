import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as wlTAGS} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

import {WebLogoViewer} from '../viewers/web-logo-viewer';

import {_package} from '../package';

/** Used in Macromolecule column tooltip */
export class MacromoleculeColumnWidget extends DG.Widget {
  private viewed: boolean = false;

  private readonly seqCol: DG.Column<string>;

  private wlViewer: WebLogoViewer | null = null;

  constructor(seqCol: DG.Column<string>) {
    super(ui.divV([]));

    this.seqCol = seqCol;
  }

  async init(): Promise<void> {
    const sh = SeqHandler.forColumn(this.seqCol);
    const pkgTooltipWebLogo = _package.properties.TooltipWebLogo;
    const colTooltipWebLogo = this.seqCol.getTag(wlTAGS.tooltipWebLogo);

    if (pkgTooltipWebLogo !== false && !['false', 'off', 'disable', 'disabled'].includes(colTooltipWebLogo)) {
      this.wlViewer = await this.seqCol.dataFrame.plot.fromType('WebLogo', {
        sequenceColumnName: this.seqCol.name,
        backgroundColor: 0x00000000,
        positionHeight: 'Entropy',
        positionWidth: (sh.getAlphabetIsMultichar() ? 24 : 16),
        fixWidth: true,
        fitArea: false,
        // maxHeight: 100,
        // minHeight: 25,
      }) as unknown as WebLogoViewer;
      this.wlViewer.root.style.height = `50px`;

      this.root.appendChild(this.wlViewer.root);
      this.root.style.width = '100%';
    }
  }

  override detach() {
    if (this.wlViewer) {
      this.wlViewer.detach();
      this.wlViewer = null;
    }
    super.detach();
  }
}
