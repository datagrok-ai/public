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

  const scaffoldColumnChoice = ui.input.choice('Scaffold column', {value: scaffoldColName, items: [NONE].concat([...columnsSet].sort()),
    onValueChanged: (value) => {
      col.temp[SCAFFOLD_COL] = value === NONE ? null : value;
      col.dataFrame?.fireValuesChanged();
    }});
  scaffoldColumnChoice.setTooltip('Align structures to a scaffold defined in another column');

  const highlightScaffoldsCheckbox = ui.input.bool('Highlight scaffold',
    {value: col?.temp && col.temp[HIGHLIGHT_BY_SCAFFOLD_COL] === 'true',
      onValueChanged: (value) => {
        col.temp[HIGHLIGHT_BY_SCAFFOLD_COL] = value.toString();
        col.dataFrame?.fireValuesChanged();
      }});
  highlightScaffoldsCheckbox.setTooltip('Highlight scaffold defined above');

  const regenerateCoordsCheckbox = ui.input.bool('Regen coords',
    {value: col?.temp && col.temp[REGENERATE_COORDS] === 'true',
      onValueChanged: (value) => {
        col.temp[REGENERATE_COORDS] = value.toString();
        col.dataFrame?.fireValuesChanged();
      }});
  regenerateCoordsCheckbox.setTooltip('Force regeneration of coordinates even for MOLBLOCKS');

  const moleculeFilteringChoice = ui.input.choice('Filter type',
    {value: col.tags[DG.TAGS.STRUCTURE_FILTER_TYPE] ?? StructureFilterType.Sketch, items: [StructureFilterType.Sketch, StructureFilterType.Categorical],
      onValueChanged: (value) => {
        col.tags[DG.TAGS.STRUCTURE_FILTER_TYPE] = value;
        col.tags[DG.TAGS.IGNORE_CUSTOM_FILTER] = (value == StructureFilterType.Categorical).toString();
      }});
  moleculeFilteringChoice.setTooltip('Sketch a molecule, or use them as categories in a filter');

  const showStructures = ui.input.bool('Structures',
    {value: col.tags[DG.TAGS.CELL_RENDERER] == DG.SEMTYPE.MOLECULE,
      onValueChanged: (value) => col.tags[DG.TAGS.CELL_RENDERER] = value ? DG.SEMTYPE.MOLECULE : DG.TYPE.STRING});

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
    col.dataFrame?.fireValuesChanged();
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
