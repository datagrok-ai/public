/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as mmcrTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';

import {MmcrTemps, rendererSettingsChangedState} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';

import {_package} from '../package';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

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

  let fontSize: number | null = (_package.properties ? _package.properties.fontSize : 12);
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
    if (units === NOTATION.HELM || units === NOTATION.CUSTOM)
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


  const sequenceConfigInputs = ui.inputs(inputsArray);
  return new DG.Widget(sequenceConfigInputs);
}
