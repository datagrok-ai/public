import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMolfilesFromSingleSeq} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {TAGS as mmcrTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';

import {
  Temps as mmcrTemps, rendererSettingsChangedState, Temps
} from '../utils/cell-renderer-consts';

import {_package} from '../package';
import {max} from 'rxjs/operators';


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

  let maxMonomerLength: number | null = (_package.properties ? _package.properties.MaxMonomerLength : 4);
  if (mmcrTAGS.maxMonomerLength in col.tags) {
    const v = parseInt(col.getTag(mmcrTAGS.maxMonomerLength));
    maxMonomerLength = !isNaN(v) ? v : maxMonomerLength;
  }
  if (Temps.maxMonomerLength in col.temp) {
    const v = parseInt(col.temp[Temps.maxMonomerLength]);
    maxMonomerLength = !isNaN(v) ? v : maxMonomerLength;
  }
  const maxMonomerLengthInput = ui.input.int('Max Monomer Length', {
    value: maxMonomerLength!,
    nullable: true, min: 1, max: 50, step: 1,
    onValueChanged: () => {
      if (maxMonomerLengthInput.value == 0)
        setTimeout(() => { maxMonomerLengthInput.value = null!; }, 0);
      else {
        const value = maxMonomerLengthInput.value ?? '';
        const tagValue = value == null ? '' : value.toString();
        col.temp[Temps.maxMonomerLength] = tagValue;
        col.temp[Temps.rendererSettingsChanged] = rendererSettingsChangedState.true;
        col.dataFrame.fireValuesChanged();
      }
    },
    tooltipText: `The max length of monomer symbol displayed without shortening, empty to no limit`
  });

  const gapLengthInput = ui.input.int('Monomer Margin', {
    value: col.temp[mmcrTemps.gapLength] ?? 0,
    onValueChanged: () => {
      col.temp[mmcrTemps.gapLength] = gapLengthInput.value;
      col.temp[mmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'The size of margin between monomers (in pixels)'
  });

  const colorCodeInput = ui.input.bool('Color Code', {
    value: (col?.temp['color-code'] != null) ? col.temp['color-code'] : true,
    onValueChanged: () => {
      col.temp['color-code'] = colorCodeInput.value;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'Color code'
  });

  const referenceSequenceInput = ui.input.string('Reference Sequence', {
    value: (col?.temp['reference-sequence'] != null) ? col?.temp['reference-sequence'] : '',
    nullable: true,
    onValueChanged: () => {
      col.temp['reference-sequence'] = referenceSequenceInput.value;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'Reference sequence is not empty, then the sequence will be render ' + '\n' +
      'as a difference from the reference sequence'
  });

  const compareWithCurrentInput = ui.input.bool('Compare with current', {
    value: (col?.temp['compare-with-current'] != null) ? col.temp['compare-with-current'] : true,
    onValueChanged: () => {
      col.temp['compare-with-current'] = compareWithCurrentInput.value;
      col.dataFrame.fireValuesChanged();
    },
    tooltipText: 'When on, all sequences get rendered in the "diff" mode'
  });

  const rdKitInputs = ui.inputs([
    maxMonomerLengthInput,
    gapLengthInput,
    referenceSequenceInput,
    colorCodeInput,
    compareWithCurrentInput,
  ]);

  return new DG.Widget(rdKitInputs);
}

/**
 * 3D representation widget of macromolecule.
 *
 * @export
 * @param {DG.Cell} macroMolecule macromolecule cell.
 * @param {any[]} monomersLibObject
 * @return {Promise<DG.Widget>} Widget.
 */
export async function representationsWidget(macroMolecule: DG.Cell, monomersLibObject: any[]): Promise<DG.Widget> {
  const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

  let widgetHost;
  let molBlock3D = '';
  try {
    try {
      const _atomicCodes = getMolfilesFromSingleSeq(macroMolecule, monomersLibObject);
      const result = '';//await getMacroMol(atomicCodes!);
      const molBlock2D = result[0];
      molBlock3D = (await grok.functions.call('Bio:Embed', {molBlock2D})) as unknown as string;
    } catch (e) {
      console.warn(e);
    }

    try {
      molBlock3D = molBlock3D.replaceAll('\\n', '\n');
      const stringBlob = new Blob([molBlock3D], {type: 'text/plain'});
      const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});

      //@ts-ignore
      const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'});
      //@ts-ignore
      stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(comp: NGL.StructureComponent) {
        stage.setSize(300, 300);
        comp.addRepresentation('ball+stick');
        comp.autoView();
      });
      const sketch = grok.chem.svgMol(molBlock3D);
      const panel = ui.divH([sketch]);

      widgetHost = ui.div([panel, nglHost]);
    } catch (e) {
      widgetHost = ui.divText('Couldn\'t get peptide structure');
    }
  } catch (e) {
    widgetHost = ui.divText('Couldn\'t get peptide structure');
  }
  pi.close();
  return new DG.Widget(widgetHost);
}
