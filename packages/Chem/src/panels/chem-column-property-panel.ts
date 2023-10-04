import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ALIGN_BY_SCAFFOLD_TAG, SCAFFOLD_COL, REGENERATE_COORDS, HIGHLIGHT_BY_SCAFFOLD_COL } from '../constants';

enum StructureFilterType {
  Sketch = 'Sketch',
  Categorical = 'Categorical'
}

/**
 * Ability to define the following:
 * - Scaffold column
 * - Substructure to highlight
 * - Filter type
 * - Show rendered structures or smiles
 *  */
export function getMolColumnPropertyPanel(col: DG.Column): DG.Widget {
  const NONE = 'None';
  const scaffoldColName = col.temp[SCAFFOLD_COL] ?? NONE;

  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns as any).filter(
    (c: any) => c.semType === DG.SEMTYPE.MOLECULE).map((c: any) => c.name);
  const columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  const scaffoldColumnChoice = ui.choiceInput('Scaffold column',
    scaffoldColName,
    [NONE].concat([...columnsSet].sort()),
    (s: string) => {
      col.temp[SCAFFOLD_COL] = s === NONE ? null : s;
      col.dataFrame.fireValuesChanged();
    });
  scaffoldColumnChoice.setTooltip('Align structures to a scaffold defined in another column');

  const highlightScaffoldsCheckbox = ui.boolInput('Highlight scaffold',
    col?.temp && col.temp[HIGHLIGHT_BY_SCAFFOLD_COL] === 'true',
    (v: any) => {
      col.temp[HIGHLIGHT_BY_SCAFFOLD_COL] = v.toString();
      col.dataFrame.fireValuesChanged();
    });
  highlightScaffoldsCheckbox.setTooltip('Highlight scaffold defined above');

  const regenerateCoordsCheckbox = ui.boolInput('Regen coords',
    col?.temp && col.temp[REGENERATE_COORDS] === 'true',
    (v: any) => {
      col.temp[REGENERATE_COORDS] = v.toString();
      col.dataFrame.fireValuesChanged();
    });
  regenerateCoordsCheckbox.setTooltip('Force regeneration of coordinates even for MOLBLOCKS');

  const moleculeFilteringChoice = ui.choiceInput('Filter type',
    col.tags[DG.TAGS.STRUCTURE_FILTER_TYPE] ?? StructureFilterType.Sketch,
    [StructureFilterType.Sketch, StructureFilterType.Categorical],
    (s: string) => {
      col.tags[DG.TAGS.STRUCTURE_FILTER_TYPE] = s;
      col.tags[DG.TAGS.IGNORE_CUSTOM_FILTER] = (s == StructureFilterType.Categorical).toString();
    },
  );
  moleculeFilteringChoice.setTooltip('Sketch a molecule, or use them as categories in a filter');

  const showStructures = ui.boolInput('Structures',
    col.tags['cell.renderer'] == DG.SEMTYPE.MOLECULE,
    (v: boolean) => col.tags['cell.renderer'] = v ? DG.SEMTYPE.MOLECULE : DG.TYPE.STRING);

  const rdKitInputs = ui.form([
    showStructures,
    scaffoldColumnChoice,
    highlightScaffoldsCheckbox,
    regenerateCoordsCheckbox,
    moleculeFilteringChoice,
  ]);
  const sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
  sketcher.syncCurrentObject = false;
  sketcher.setMolFile(col.tags[ALIGN_BY_SCAFFOLD_TAG]);
  sketcher.onChanged.subscribe((_: any) => {
    const molFile = sketcher.getMolFile();
    col.tags[ALIGN_BY_SCAFFOLD_TAG] = molFile;
    col.temp[ALIGN_BY_SCAFFOLD_TAG] = molFile;
    col.dataFrame.fireValuesChanged();
  });
  sketcher.root.classList.add('ui-input-editor');
  sketcher.root.style.marginTop = '3px';
  const scaffoldLabel = ui.label('Scaffold');
  scaffoldLabel.classList.add('ui-input-label');
  const scaffoldInput = ui.divH([
    scaffoldLabel,
    sketcher.root,
  ]);
  scaffoldInput.className = 'ui-input-root';
  rdKitInputs.append(scaffoldInput);
  return new DG.Widget(ui.box(rdKitInputs));
}
