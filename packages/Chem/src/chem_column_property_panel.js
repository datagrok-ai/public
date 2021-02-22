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

  let choiceScaffoldColumn = ui.choiceInput(
    'Scaffold column',
     scaffoldColName,
    [NONE].concat([...columnsSet].sort()));
  choiceScaffoldColumn.onChanged(_ => {
    const scaffoldColName = choiceScaffoldColumn.stringValue;
    col.tags['scaffold-col'] = scaffoldColName === NONE ? null : scaffoldColName;
  });
  return ui.div([
    ui.inputs([
      choiceScaffoldColumn,
      ui.boolInput('Highlight scaffolds', col?.tags && col.tags['highlight-scaffold'] === 'true',
        v => { col.tags['highlight-scaffold'] = v.toString(); }),
      ui.boolInput('Regenerate coords', col?.tags && col.tags['regenerate-coords'] === 'true',
        v => { col.tags['regenerate-coords'] = v.toString(); })
    ])
  ]);

  // TODO: react to changing the tags, GROK-8093

}