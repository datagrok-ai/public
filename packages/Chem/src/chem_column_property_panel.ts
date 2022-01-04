import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let subscr: any = null;

export function getMolColumnPropertyPanel(col: DG.Column) {
  const NONE = 'None';
  let scaffoldColName = null;
  if (col?.temp && col.temp['scaffold-col']) {
    scaffoldColName = col.temp['scaffold-col'];
  } else {
    scaffoldColName = NONE;
  }
  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns).filter(
    (c: any) => c.semType === DG.SEMTYPE.MOLECULE).map((c: any) => c.name);
  const columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  const scaffoldColumnChoice = ui.choiceInput(
    'Scaffold column',
    scaffoldColName,
    [NONE].concat([...columnsSet].sort()));
  scaffoldColumnChoice.onChanged((_: any) => {
    const scaffoldColName = scaffoldColumnChoice.stringValue;
    col.temp['scaffold-col'] = scaffoldColName === NONE ? null : scaffoldColName;
    col.dataFrame.fireValuesChanged();
  });
  const highlightScaffoldsCheckbox = ui.boolInput(
    'Highlight from column', col?.temp && col.temp['highlight-scaffold'] === 'true',
    (v: any) => {
      col.temp['highlight-scaffold'] = v.toString();
      col.dataFrame.fireValuesChanged();
    });
  const regenerateCoordsCheckbox = ui.boolInput(
    'Regenerate coords', col?.temp && col.temp['regenerate-coords'] === 'true',
    (v: any) => {
      col.temp['regenerate-coords'] = v.toString();
      col.dataFrame.fireValuesChanged();
    });

  const matchMoleculeFilteringToDropdown = (v: string) => {
    if (v === 'categorical') return 'Categorical';
    if (v === 'sketching') return 'Sketching';
    return 'Dynamic';
  };

  const matchDropdownToMoleculeFiltering = (v: string) => {
    if (v === 'Categorical') {
      col.temp['.molecule-filtering'] = 'categorical';
    } else if (v === 'Sketching') {
      col.temp['.molecule-filtering'] = 'sketching';
    } else {
      col.temp.delete('.molecule-filtering');
    }
  };

  const moleculeFilteringChoice = ui.choiceInput('Filter Type',
    matchMoleculeFilteringToDropdown(col.temp['.molecule-filtering']),
    ['Dynamic', 'Categorical', 'Sketching']);
  moleculeFilteringChoice.onChanged((_: any) => {
    const v = moleculeFilteringChoice.stringValue;
    matchDropdownToMoleculeFiltering(v);
  });

  subscr?.unsubscribe();
  subscr = col.dataFrame.onMetadataChanged.subscribe((a: any) => {
    // Handling scaffold column
    const scaffoldColumnChoiceValue = scaffoldColumnChoice.stringValue;
    const scaffoldColumnTag = col.temp && col.temp['scaffold-col'] ? col.temp['scaffold-col'] : NONE;
    if (scaffoldColumnChoiceValue !== scaffoldColumnTag) {
      const htmlSelectElement = scaffoldColumnChoice.root.children[1] as HTMLSelectElement;
      if (scaffoldColumnTag === NONE) {
        htmlSelectElement.value = NONE;
      } else if (columnsSet.has(scaffoldColumnTag)) {
        htmlSelectElement.value = scaffoldColumnTag;
      } else {
        // TODO: handle a selection of a non-molecule column
      }
    }
    // handling highlight scaffolds selection
    const highlightScaffoldsCheckboxValue = highlightScaffoldsCheckbox.value;
    const highlightScaffoldsTagPresent = col.temp && col.temp['highlight-scaffold'] === 'true';
    const highlightScaffoldsCheckboxElement = highlightScaffoldsCheckbox.root.children[1] as HTMLInputElement;
    if (highlightScaffoldsCheckboxValue != highlightScaffoldsTagPresent) {
      highlightScaffoldsCheckboxElement.checked = highlightScaffoldsTagPresent;
    }
    // handling regenerate coords selection
    const regenerateCoordsCheckboxValue = regenerateCoordsCheckbox.value;
    const regenerateCoordsTagPresent = col.temp && col.temp['regenerate-coords'] === 'true';
    const regenerateCoordsCheckboxElement = regenerateCoordsCheckbox.root.children[1] as HTMLInputElement;
    if (regenerateCoordsCheckboxValue != regenerateCoordsTagPresent) {
      regenerateCoordsCheckboxElement.checked = regenerateCoordsTagPresent;
    }
    // handling molecule filtering choice value
    const moleculeFilteringChoiceValue = moleculeFilteringChoice.stringValue;
    const moleculeFilteringTag = matchMoleculeFilteringToDropdown(col?.temp['.molecule-filtering']);
    const moleculeFilteringChoiceElement = moleculeFilteringChoice.root.children[1] as HTMLSelectElement;
    if (moleculeFilteringChoiceValue != moleculeFilteringTag) {
      moleculeFilteringChoiceElement.value = moleculeFilteringTag;
    }
  });

  const widget = new DG.Widget(ui.div([
    ui.inputs([
      scaffoldColumnChoice,
      highlightScaffoldsCheckbox,
      regenerateCoordsCheckbox,
      moleculeFilteringChoice,
    ]),
  ]));

  widget.subs.push(subscr);

  return widget;
}
