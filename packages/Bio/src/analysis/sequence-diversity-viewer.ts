import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import $ from 'cash-dom';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import {SequenceSearchBaseViewer} from './sequence-search-base-viewer';
import {getMonomericMols} from '../calculations/monomerLevelMols';
import {updateDivInnerHTML} from '../utils/ui-utils';
import {Subject} from 'rxjs';

export class SequenceDiversityViewer extends SequenceSearchBaseViewer {
  diverseColumnLabel: string | null; // Use postfix Label to prevent activating table column selection editor

  renderMolIds: number[] | null = null;
  columnNames = [];
  computeCompleted = new Subject<boolean>();

  constructor() {
    super('diversity');
    this.diverseColumnLabel = this.string('diverseColumnLabel', null);
  }

  async render(computeData = true): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.dataFrame) {
      if (computeData && this.moleculeColumn) {
        const monomericMols = await getMonomericMols(this.moleculeColumn);
        //need to create df to calculate fingerprints
        const monomericMolsDf = DG.DataFrame.fromColumns([monomericMols]);
        this.renderMolIds = await grok.functions.call('Chem:callChemDiversitySearch', {
          col: monomericMols,
          metricName: this.distanceMetric,
          limit: this.limit,
          fingerprint: this.fingerprint
        });
        const diverseColumnName: string = this.diverseColumnLabel != null ? this.diverseColumnLabel :
          `diverse (${this.moleculeColumnName})`;
        const resCol = DG.Column.string(diverseColumnName, this.renderMolIds!.length)
          .init((i) => this.moleculeColumn?.get(this.renderMolIds![i]));
        resCol.semType = DG.SEMTYPE.MACROMOLECULE;
        this.tags.forEach((tag) => resCol.setTag(tag, this.moleculeColumn!.getTag(tag)));
        const resDf = DG.DataFrame.fromColumns([resCol]);
        updateDivInnerHTML(this.root, resDf.plot.grid().root);
        this.computeCompleted.next(true);
      }
    }
  }
}
