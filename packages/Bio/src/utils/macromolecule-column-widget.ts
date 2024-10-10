import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as wlTAGS} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {WebLogoViewer} from '../viewers/web-logo-viewer';

import {_package} from '../package';

/** Used in Macromolecule column tooltip */
export class MacromoleculeColumnWidget extends DG.Widget {
  private viewed: boolean = false;

  private wlViewer: WebLogoViewer | null = null;

  constructor(
    private readonly seqCol: DG.Column<string>,
    private readonly seqHelper: ISeqHelper,
  ) {
    super(ui.divV([]));
  }

  async init(): Promise<void> {
    const sh = this.seqHelper.getSeqHandler(this.seqCol);
    const pkgTooltipWebLogo = _package.properties.tooltipWebLogo;
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
        positionNames: '', // to ensure position names by default
        endPositionName: '50', // limit WebLogo for visible monomers
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
