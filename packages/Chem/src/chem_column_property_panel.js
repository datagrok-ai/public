function getMolColumnPropertyPanel(col) {

  const names = col.dataFrame.columns.names();
  return ui.div([
    ui.inputs([
      ui.choiceInput('Scaffold column', 'None', ['None'].concat(names),
        v => {}),
      ui.boolInput('Highlight scaffolds', col?.tags && col.tags['highlight-scaffold'] === 'true',
        v => {}),
      ui.boolInput('Regenerate coords', col?.tags && col.tags['regenerate-coords'] === 'true',
        v => {})
    ])
  ]);

}