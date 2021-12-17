class WidgetsPackage extends DG.Package {

  //name: Single Choice
  //description: A filter that lets you select exactly one category
  //tags: filter
  //output: filter result
  radioButtonFilter() {
    return new RadioButtonFilter();
  }

  //name: Multi Choice
  //description: A filter that works with columns of multi-value cells (such as lists of identifiers)
  //tags: filter
  //output: filter result
  multiValueFilter() {
    return new MultiValueFilter();
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