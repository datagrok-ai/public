import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ALIGN_BY_SCAFFOLD_TAG, ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG, CHEM_ATOM_PICKER_TAG, SCAFFOLD_COL,
  SCAFFOLD_COL_SYNC, HIGHLIGHT_BY_SCAFFOLD_COL, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC, REGENERATE_COORDS,
  REGENERATE_COORDS_SYNC, getSyncTag, setSyncTag} from '../constants';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';

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
  const scaffoldColName = getSyncTag(col, SCAFFOLD_COL_SYNC, SCAFFOLD_COL) ?? NONE;

  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns as any).filter(
    (c: any) => c.semType === DG.SEMTYPE.MOLECULE).map((c: any) => c.name);
  const columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  const scaffoldColumnChoice = ui.input.choice('Scaffold column',
    {value: scaffoldColName, items: [NONE].concat([...columnsSet].sort()),
      onValueChanged: (value) => {
        setSyncTag(col, SCAFFOLD_COL_SYNC, SCAFFOLD_COL, value === NONE ? null : value);
        col.dataFrame?.fireValuesChanged();
      }});
  scaffoldColumnChoice.setTooltip('Align structures to a scaffold defined in another column');

  const highlightScaffoldInitValue = getSyncTag(col, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC, HIGHLIGHT_BY_SCAFFOLD_COL);
  const highlightScaffoldsCheckbox = ui.input.bool('Highlight scaffold',
    {value: highlightScaffoldInitValue === 'true',
      onValueChanged: (value) => {
        setSyncTag(col, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC, HIGHLIGHT_BY_SCAFFOLD_COL, value.toString());
        col.dataFrame?.fireValuesChanged();
      }});
  highlightScaffoldsCheckbox.setTooltip('Highlight scaffold defined above');

  const backupSubstructCol = col.temp[ChemTemps.SUBSTRUCT_BACKUP_COL];
  const substructCol = col.temp[ChemTemps.SUBSTRUCT_COL];
  const hasHighlightCol = !!substructCol;
  const highlightRGroupsCheckbox = ui.input.bool('Highlight R-Groups', {
    value: hasHighlightCol,
    onValueChanged: (enabled) => {
      if (enabled && backupSubstructCol)
        col.temp[ChemTemps.SUBSTRUCT_COL] = backupSubstructCol;
      else
        delete col.temp[ChemTemps.SUBSTRUCT_COL];
      col.dataFrame?.fireValuesChanged();
    },
  });

  const regenerateCoordsInitValue = getSyncTag(col, REGENERATE_COORDS_SYNC, REGENERATE_COORDS);
  const regenerateCoordsCheckbox = ui.input.bool('Regen coords',
    {value: regenerateCoordsInitValue === 'true',
      onValueChanged: (value) => {
        setSyncTag(col, REGENERATE_COORDS_SYNC, REGENERATE_COORDS, value.toString());
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

  // In-grid atom picker — drag / Alt+drag / Alt+click on molecule cells to
  // highlight parts of the structure. DISABLED by default on every
  // molecule column; checking writes `chem-atom-picker = 'true'` which
  // the cell renderer checks at mousedown to enable picker handling AND
  // which the ISubstructProvider checks at render time to show the picks.
  // Unchecking writes 'false' (or any other value); only the explicit
  // string 'true' enables the picker.
  //
  // We write via BOTH `col.tags[..]` assignment AND `col.setTag(..)`, and
  // read via both on the renderer side. Datagrok's two tag APIs don't
  // always agree about what was last written for a given tag name; using
  // both leaves either read path correct.
  const atomPickerCheckbox = ui.input.bool('Interactive Atom Selection',
    {
      value: ((col.tags as any)[CHEM_ATOM_PICKER_TAG] ??
        col.getTag(CHEM_ATOM_PICKER_TAG)) === 'true',
      onValueChanged: (value) => {
        const v = value ? 'true' : 'false';
        (col.tags as any)[CHEM_ATOM_PICKER_TAG] = v;
        col.setTag(CHEM_ATOM_PICKER_TAG, v);
        // Force a grid repaint so highlights appear/disappear immediately
        // and the renderer sees the new tag value on the next drag.
        grok.shell.tv?.grid.invalidate();
      },
    });
  atomPickerCheckbox.setTooltip(
    'Drag on a cell to select atoms. Alt+drag to add to the existing selection. ' +
    'Alt+click on an atom to toggle it.');

  const rdKitInputs = ui.form([
    showStructures,
    scaffoldColumnChoice,
    highlightScaffoldsCheckbox,
    ...(backupSubstructCol ? [highlightRGroupsCheckbox] : []),
    regenerateCoordsCheckbox,
    moleculeFilteringChoice,
    atomPickerCheckbox,
  ]);
  const sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
  sketcher.syncCurrentObject = false;
  sketcher.setMolFile(col.tags[ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG] || col.tags[ALIGN_BY_SCAFFOLD_TAG]);
  sketcher.onChanged.subscribe((_: any) => {
    const molFile = sketcher.getMolFile();
    if (col.tags[ALIGN_BY_SCAFFOLD_TAG])
      delete col.tags[ALIGN_BY_SCAFFOLD_TAG];
    col.tags[ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG] = molFile;
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
