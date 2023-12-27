import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as wlTAGS} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {WebLogoViewer} from '../viewers/web-logo-viewer';

import {_package} from '../package';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

/** Used in Macromolecule column tooltip */
export class MacromoleculeColumnWidget extends DG.Widget {
  private viewed: boolean = false;

  private readonly seqCol: DG.Column<string>;

  private wlViewer: WebLogoViewer;

  constructor(seqCol: DG.Column<string>) {
    super(ui.divV([]));

    this.seqCol = seqCol;
  }

  async init(): Promise<void> {
    const uh = UnitsHandler.getOrCreate(this.seqCol);
    const pkgTooltipWebLogo = _package.properties.TooltipWebLogo;
    const colTooltipWebLogo = this.seqCol.getTag(wlTAGS.tooltipWebLogo);

    if (pkgTooltipWebLogo !== false && !['false', 'off', 'disable', 'disabled'].includes(colTooltipWebLogo)) {
      this.wlViewer = await this.seqCol.dataFrame.plot.fromType('WebLogo', {
        sequenceColumnName: this.seqCol.name,
        backgroundColor: 0x00000000,
        positionHeight: 'Entropy',
        positionWidth: (uh.getAlphabetIsMultichar() ? 24 : 16),
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
    this.wlViewer.detach();
    super.detach();
  }
}
