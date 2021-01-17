class WidgetsPackage extends DG.Package {

  //name: RadioButtonFilter
  //description: A filter that lets you select exactly one category
  //tags: filter
  //output: filter result
  radioButtonFilter() {
    return new RadioButtonFilter();
  }

  //name: TimeWidget
  //description: Shows current time
  //output: widget result
  timeWidget() {
    return new TimeWidget();
  }

  //name: SmilesLengthWidgetPanel
  //input: string smiles = CN1C=NC2=C1C(=O)N(C(=O)N2C)C {semType: Molecule}
  //tags: panel
  //output: widget result
  smilesLengthWidgetPanel(smiles) {
    return new SmilesLengthWidget().apply({smiles: smiles});
  }

  //name: SmilesLengthWidget
  //output: widget result
  smilesLengthWidget() {
    return new SmilesLengthWidget();
  }
}