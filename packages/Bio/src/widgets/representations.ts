/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as mmcrTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';

import {MmcrTemps, rendererSettingsChangedState} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';

import {_package} from '../package';
import {NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getColumnAnnotations, clearAnnotations} from '../utils/annotations/annotation-manager';
import {getBioRendererDefaults} from '../utils/renderer-defaults';
import {ISeqGraphInfoHelm} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {AnnotationCategory} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';

/** Whether the sequence renderer can draw the column's connections: BILN, or a notation whose splitter gives some
 * rows connections, unless the column is drawn by the HELM renderer */
function canShowConnections(col: DG.Column<string>): boolean {
  if (col.meta.units === NOTATION.BILN) return true;
  if (col.getTag(DG.TAGS.CELL_RENDERER) === 'helm') return false;
  try {
    const sh = _package.seqHelper.getSeqHandler(col);
    for (let rowIdx = 0; rowIdx < Math.min(col.length, 20); ++rowIdx) {
      const gi = sh.getSplitted(rowIdx).graphInfo as ISeqGraphInfoHelm | undefined;
      if (gi?.connections?.length || gi?.hbonds?.length) return true;
    }
  } catch {
    // the column cannot be split: nothing to draw
  }
  return false;
}

/**
 * @export
 * @param {DG.Column} col macromolecule cell.
 * @return {Promise<DG.Widget>} Widget.
 */
export function getMacromoleculeColumnPropertyPanel(col: DG.Column): DG.Widget {
  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns as any).filter(
    (c: any) => c.semType === DG.SEMTYPE.MOLECULE).map((c: any) => c.name);
  const columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  let maxMonomerLength: number | null = (_package.properties ? _package.properties.maxMonomerLength : 4);
  if (mmcrTAGS.maxMonomerLength in col.tags) {
    const v = parseInt(col.getTag(mmcrTAGS.maxMonomerLength));
    maxMonomerLength = !isNaN(v) ? v : maxMonomerLength;
  }
  if (MmcrTemps.maxMonomerLength in col.temp) {
    const v = parseInt(col.temp[MmcrTemps.maxMonomerLength]);
    maxMonomerLength = !isNaN(v) ? v : maxMonomerLength;
  }
  const maxMonomerLengthInput = ui.input.int('Max Monomer Length', {
    value: maxMonomerLength!,
    nullable: true, min: 1, max: 50, step: 1,
    onValueChanged: (value) => {
      if (value == 0)
        setTimeout(() => { maxMonomerLengthInput.value = null!; }, 0);
      else {
        const newValue = value ?? '';
        const tagValue = newValue == null ? '' : newValue.toString();
        col.temp[MmcrTemps.maxMonomerLength] = tagValue;
        col.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
        col.dataFrame.fireValuesChanged();
      }
    },
    tooltipText: `The max length of monomer symbol displayed without shortening, empty to no limit`
  });

  let fontSize: number | null = getBioRendererDefaults().fontSize;
  if (MmcrTemps.fontSize in col.temp && !!col.temp[MmcrTemps.fontSize] && !isNaN(col.temp[MmcrTemps.fontSize]))
    fontSize = col.temp[MmcrTemps.fontSize];

  const fontSizeInput = ui.input.int('Font Size', {
    value: fontSize!,
    nullable: true, min: 1, max: 50, step: 1,
    onValueChanged: (value) => {
      if (value && value > 0) {
        const newValue = value ?? 12;
        col.temp[MmcrTemps.fontSize] = newValue;
        col.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
        col.dataFrame.fireValuesChanged();
      }
    },
    tooltipText: `The font size of monomer symbol in sequence renderer`
  });

  const gapLengthInput = ui.input.int('Monomer Margin', {
    value: col.temp[MmcrTemps.gapLength] ?? 0,
    onValueChanged: (value) => {
      col.temp[MmcrTemps.gapLength] = value;
      col.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'The size of margin between monomers (in pixels)'
  });

  const colorCodeInput = ui.input.bool('Color Code', {
    value: (col?.temp['color-code'] != null) ? col.temp['color-code'] : true,
    onValueChanged: (value) => {
      col.temp['color-code'] = value;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'Color code'
  });

  const referenceSequenceInput = ui.input.string('Reference Sequence', {
    value: (col?.temp['reference-sequence'] != null) ? col?.temp['reference-sequence'] : '',
    nullable: true,
    onValueChanged: (value) => {
      col.temp['reference-sequence'] = value;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'Reference sequence is not empty, then the sequence will be render ' + '\n' +
      'as a difference from the reference sequence'
  });

  const compareWithCurrentInput = ui.input.bool('Compare with current', {
    value: (col?.temp['compare-with-current'] != null) ? col.temp['compare-with-current'] : true,
    onValueChanged: (value) => {
      col.temp['compare-with-current'] = value;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'When on, all sequences get rendered in the "diff" mode'
  });

  const shouldShowMultilineToggle = (): boolean => {
    const units = col.meta.units;

    // Don't show for formats that have their own complex renderers (like Helm).
    if (units === NOTATION.HELM)
      return false;

    // For all other cases, including 'UN' (non-canonical), 'fasta', and 'separator' show the multiline toggle.
    return true;
  };

  let renderMultilineInput = null;
  if (shouldShowMultilineToggle()) {
    renderMultilineInput = ui.input.bool('Multiline Rendering', {
      value: col.getTag('renderMultiline') === 'true',
      onValueChanged: (value) => {
        col.tags['renderMultiline'] = value ? 'true' : 'false';
        col.dataFrame.fireValuesChanged();
      },
      tooltipText: 'Render sequences across multiple lines when they exceed cell width'
    });
  }

  // Rows with connections (BILN cyclizations, HELM bonds shown as a sequence) draw them as lines: the column's tag,
  // else the package default
  let showConnectionsInput: DG.InputBase<boolean> | null = null;
  if (canShowConnections(col)) {
    const tag = col.getTag(mmcrTAGS.showConnections);
    showConnectionsInput = ui.input.bool('Show Connections', {
      value: tag === 'true' || tag === 'false' ? tag === 'true' : getBioRendererDefaults().showConnections,
      onValueChanged: (value) => {
        col.setTag(mmcrTAGS.showConnections, value ? 'true' : 'false');
        col.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
        col.dataFrame.fireValuesChanged();
      },
      tooltipText: 'Draw the bonds of cyclized and multi-chain sequences as lines over the monomers when the rows ' +
        'are tall enough, instead of marks like C(1,3)',
    });
  }

  const inputsArray = [
    fontSizeInput,
    maxMonomerLengthInput,
    gapLengthInput,
    referenceSequenceInput,
    colorCodeInput,
    compareWithCurrentInput,
  ];

  if (renderMultilineInput)
    inputsArray.push(renderMultilineInput);
  if (showConnectionsInput)
    inputsArray.push(showConnectionsInput);


  const sequenceConfigInputs = ui.inputs(inputsArray);

  // --- Annotations section ---
  const annotations = getColumnAnnotations(col);
  const annotationsDiv = ui.divV([]);

  if (annotations.length > 0) {
    const scheme = col.getTag(bioTAGS.numberingScheme);
    const structAnnots = annotations.filter((a) => a.category === AnnotationCategory.Structure);
    const liabAnnots = annotations.filter((a) => a.category === AnnotationCategory.Liability);

    if (scheme)
      annotationsDiv.append(ui.divText(`Numbering: ${scheme}`, {style: {fontSize: '12px', marginBottom: '4px'}}));

    if (structAnnots.length > 0) {
      const regionNames = structAnnots.map((a) => a.name).join(', ');
      annotationsDiv.append(ui.divText(`Regions: ${regionNames}`, {style: {fontSize: '12px', marginBottom: '4px'}}));
    }

    if (liabAnnots.length > 0) {
      const totalHits = liabAnnots.reduce((sum, a) => {
        const match = a.description?.match(/\((\d+) hits\)/);
        return sum + (match ? parseInt(match[1]) : 0);
      }, 0);
      annotationsDiv.append(ui.divText(`Liabilities: ${liabAnnots.length} rules (${totalHits} hits)`, {style: {fontSize: '12px', marginBottom: '4px'}}));
    }

    const manageBtn = ui.button('Manage', () => {
      import('../utils/annotations/annotation-manager-ui').then((m) => m.showAnnotationManagerDialog());
    });
    const clearBtn = ui.button('Clear All', () => {
      clearAnnotations(col.dataFrame, col);
      col.dataFrame.fireValuesChanged();
    });
    annotationsDiv.append(ui.divH([manageBtn, clearBtn], {style: {gap: '4px', marginTop: '4px'}}));
  } else {
    annotationsDiv.append(ui.divText('No annotations', {style: {fontSize: '12px', color: '#888'}}));
    const scanBtn = ui.button('Scan Liabilities', () => {
      import('../utils/annotations/liability-scanner-ui').then((m) => m.showLiabilityScannerDialog());
    });
    const numberBtn = ui.button('Apply Numbering', () => {
      import('../utils/annotations/numbering-ui').then((m) => m.showNumberingSchemeDialog());
    });
    annotationsDiv.append(ui.divH([scanBtn, numberBtn], {style: {gap: '4px', marginTop: '4px'}}));
  }

  const acc = ui.accordion();
  acc.addPane('Renderer Settings', () => sequenceConfigInputs);
  acc.addPane('Annotations', () => annotationsDiv);

  return new DG.Widget(acc.root);
}
