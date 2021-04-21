let subscr = null;

function getMolColumnPropertyPanel(col) {

  const NONE = 'None';
  let scaffoldColName = null;
  if (col?.tags && col.tags['scaffold-col']) {
    scaffoldColName = col.tags['scaffold-col'];
  } else {
    scaffoldColName = NONE;
  }
  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns).filter(c => c.semType === DG.SEMTYPE.MOLECULE).map(c => c.name);
  let columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  let scaffoldColumnChoice = ui.choiceInput(
    'Scaffold column',
     scaffoldColName,
    [NONE].concat([...columnsSet].sort()));
  scaffoldColumnChoice.onChanged(_ => {
    const scaffoldColName = scaffoldColumnChoice.stringValue;
    col.tags['scaffold-col'] = scaffoldColName === NONE ? null : scaffoldColName;
  });
  let highlightScaffoldsCheckbox = ui.boolInput(
    'Highlight from column', col?.tags && col.tags['highlight-scaffold'] === 'true',
    v => { col.tags['highlight-scaffold'] = v.toString(); });
  let regenerateCoordsCheckbox = ui.boolInput(
    'Regenerate coords', col?.tags && col.tags['regenerate-coords'] === 'true',
    v => { col.tags['regenerate-coords'] = v.toString(); });
    
    
    const matchMoleculeFilteringToDropdown = (v) => {
    if (v === 'categorical') return 'Categorical';
    if (v === 'sketching') return 'Sketching';
    return 'Dynamic';
  };
  
  const matchDropdownToMoleculeFiltering = (v) => {
    if (v === 'Categorical')
      col.tags['.molecule-filtering'] = 'categorical';
    else if (v === 'Sketching')
      col.tags['.molecule-filtering'] = 'sketching';
    else
      col.tags.delete('.molecule-filtering');
  };
  
  let moleculeFilteringChoice = ui.choiceInput('Filter Type',
    matchMoleculeFilteringToDropdown(col.tags['.molecule-filtering']),
    ['Dynamic', 'Categorical', 'Sketching']);
  moleculeFilteringChoice.onChanged(_ => {
    const v = moleculeFilteringChoice.stringValue;
    matchDropdownToMoleculeFiltering(v);
  });

  subscr?.unsubscribe();
  subscr = col.dataFrame.onMetadataChanged.subscribe((a) => {
    // Handling scaffold column
    let scaffoldColumnChoiceValue = scaffoldColumnChoice.stringValue;
    const scaffoldColumnTag = col.tags && col.tags['scaffold-col'] ? col.tags['scaffold-col'] : NONE;
    if (scaffoldColumnChoiceValue !== scaffoldColumnTag) {
      if (scaffoldColumnTag === NONE) {
        scaffoldColumnChoice.root.children[1].value = NONE;
      } else if (columnsSet.has(scaffoldColumnTag)) {
        scaffoldColumnChoice.root.children[1].value = scaffoldColumnTag;
      } else {
        // TODO: handle a selection of a non-molecule column
      }
    }
    // handling highlight scaffolds selection
    const highlightScaffoldsCheckboxValue = highlightScaffoldsCheckbox.value;
    const highlightScaffoldsTagPresent = col.tags && col.tags['highlight-scaffold'] === 'true';
    if (highlightScaffoldsCheckboxValue != highlightScaffoldsTagPresent) {
      highlightScaffoldsCheckbox.root.children[1].checked = highlightScaffoldsTagPresent;
    }
    // handling regenerate coords selection
    const regenerateCoordsCheckboxValue = regenerateCoordsCheckbox.value;
    const regenerateCoordsTagPresent = col.tags && col.tags['regenerate-coords'] === 'true';
    if (regenerateCoordsCheckboxValue != regenerateCoordsTagPresent) {
      regenerateCoordsCheckbox.root.children[1].checked = regenerateCoordsTagPresent;
    }
    // handling molecule filtering choice value
    const moleculeFilteringChoiceValue = moleculeFilteringChoice.stringValue;
    const moleculeFilteringTag = matchMoleculeFilteringToDropdown(col?.tags['.molecule-filtering']);
    if (moleculeFilteringChoiceValue != moleculeFilteringTag) { 
      moleculeFilteringChoice.root.children[1].value = moleculeFilteringTag;
    }
  });

  let widget = new DG.Widget(ui.div([
    ui.inputs([
      scaffoldColumnChoice,
      highlightScaffoldsCheckbox,
      regenerateCoordsCheckbox,
      moleculeFilteringChoice
    ])
  ]));

  widget.subs.push(subscr);

  return widget;
}