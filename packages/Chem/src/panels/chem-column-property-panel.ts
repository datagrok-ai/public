import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

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
  let scaffoldColName = col.temp['scaffold-col'] ?? NONE;

  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns as any).filter(
    (c: any) => c.semType === DG.SEMTYPE.MOLECULE).map((c: any) => c.name);
  const columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  const scaffoldColumnChoice = ui.choiceInput('Scaffold column',
    scaffoldColName,
    [NONE].concat([...columnsSet].sort()),
    (s: string) => {
      col.temp['scaffold-col'] = s === NONE ? null : s;
      col.dataFrame.fireValuesChanged();
    });
  scaffoldColumnChoice.setTooltip('Align structures to a scaffold defined in another column');

  const highlightScaffoldsCheckbox = ui.boolInput('Highlight from column',
    col?.temp && col.temp['highlight-scaffold'] === 'true',
    (v: any) => {
      col.temp['highlight-scaffold'] = v.toString();
      col.dataFrame.fireValuesChanged();
    });
  highlightScaffoldsCheckbox.setTooltip('Highlight scaffold defined above');

  const regenerateCoordsCheckbox = ui.boolInput('Regenerate coords',
    col?.temp && col.temp['regenerate-coords'] === 'true',
    (v: any) => {
      col.temp['regenerate-coords'] = v.toString();
      col.dataFrame.fireValuesChanged();
    });
  regenerateCoordsCheckbox.setTooltip('Force regeneration of coordinates even for MOLBLOCKS');

  const moleculeFilteringChoice = ui.choiceInput('Filter type',
    col.tags['.structure-filter-type'] ?? StructureFilterType.Sketch,
    [StructureFilterType.Sketch, StructureFilterType.Categorical],
    (s: string) => {
      col.tags['.structure-filter-type'] = s;
      col.tags['.ignore-custom-filter'] = (s == StructureFilterType.Categorical).toString();
    }
  );
  moleculeFilteringChoice.setTooltip('Sketch a molecule, or use them as categories in a filter');

  const showStructures = ui.boolInput('Show structures',
    col.tags['cell.renderer'] == DG.SEMTYPE.MOLECULE,
    (v: boolean) => col.tags['cell.renderer'] = v ? DG.SEMTYPE.MOLECULE : DG.TYPE.STRING);

  const rdKitPane = ui.div([
    ui.inputs([
      showStructures,
      scaffoldColumnChoice,
      highlightScaffoldsCheckbox,
      regenerateCoordsCheckbox,
      moleculeFilteringChoice,
    ]),
  ]);

  const panes = ui.accordion('chem-settings');
  panes.addPane('RDKit', () => rdKitPane);
  panes.addPane('Scaffold', () => {
    const sketcher = new DG.chem.chem.Sketcher();
    sketcher.syncCurrentObject = false;
    sketcher.setMolFile(col.tags['chem-scaffold']);
    sketcher.onChanged.subscribe((_) => col.tags['chem-scaffold'] = sketcher.getMolFile());
    return sketcher.root;
  });

  return new DG.Widget(panes.root);
}
