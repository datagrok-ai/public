import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getSplitterForColumn} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';

import {PeptidesModel} from '../model';
import * as C from '../utils/constants';

export function identityWidget(model: PeptidesModel): DG.Widget {
  let sequence = model.identityTemplate;
  const sequencesCol = model.df.getCol(model.settings.sequenceColumnName!);
  const templateInput = ui.stringInput('Template', sequence, async () => {
    if (isNaN(parseInt(templateInput.value))) {
      if (templateInput.value.length === 0) {
        calculateIdentityBtn.disabled = true;
        model.identityTemplate = '';
        return;
      }
      sequence = templateInput.value;
    } else {
      const rowIndex = parseInt(templateInput.value);
      if (rowIndex < 0 || rowIndex >= model.df.rowCount) {
        grok.shell.warning('Invalid row index');
        calculateIdentityBtn.disabled = true;
        return;
      }
      sequence = sequencesCol.get(rowIndex);
    }
    model.identityTemplate = templateInput.value;
    calculateIdentityBtn.disabled = false;
  }, {placeholder: 'Sequence or row index...'});
  templateInput.setTooltip('Template sequence. Can be row index, peptide ID or sequence.');

  const calculateIdentityBtn = ui.button('Calculate', async () => calculateIdentity(sequence, model));
  templateInput.fireChanged();

  return new DG.Widget(ui.div([templateInput, calculateIdentityBtn]));
}

export async function calculateIdentity(sourceSeq: string, model: PeptidesModel): Promise<DG.Column<number>> {
  const tempDf = DG.DataFrame.fromCsv(`sequence\n${sourceSeq}`);
  await grok.data.detectSemanticTypes(tempDf);
  const splitter = getSplitterForColumn(tempDf.getCol('sequence'));
  const template = splitter(sourceSeq);
  const splitSeqDf = model.splitSeqDf;
  const numPositions = splitSeqDf.columns.length;
  const positionCols: Uint32Array[] = new Array(numPositions);
  const positionEmptyCategories: number[] = new Array(numPositions);
  const categoryIndexesTemplate: number[] = new Array(numPositions);

  for (let posIdx = 0; posIdx < numPositions; ++posIdx) {
    const posCol = splitSeqDf.columns.byIndex(posIdx);
    positionCols[posIdx] = posCol.getRawData() as Uint32Array;
    positionEmptyCategories[posIdx] = posCol.categories.indexOf('');
    categoryIndexesTemplate[posIdx] = posCol.categories.indexOf(template[posIdx] ?? '');
  }

  const identityScoresCol: DG.Column<number> = model.df.columns.addNewInt(model.df.columns.getUnusedName('Identity'));
  const identityScoresData = identityScoresCol.getRawData();
  for (let rowIndex = 0; rowIndex < splitSeqDf.rowCount; ++rowIndex) {
    identityScoresData[rowIndex] = 0;
    for (let posIdx = 0; posIdx < numPositions; ++posIdx) {
      const categoryIndex = positionCols[posIdx][rowIndex];
      if (categoryIndex !== categoryIndexesTemplate[posIdx])
        ++identityScoresData[rowIndex];
    }
  }
  identityScoresCol.setTag(C.TAGS.IDENTITY_TEMPLATE, sourceSeq);
  return identityScoresCol;
}
