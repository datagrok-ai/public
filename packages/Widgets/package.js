class WidgetsPackage extends DG.Package {

  //name: RadioButtonFilter
  //description: A filter that lets you select exactly one category
  //tags: filter
  //output: filter result
  radioButtonFilter() {
    return new RadioButtonFilter();
  }

  //name: TimeWidget
  //output: widget result
  timeWidget() {
    return new TimeWidget();
  }
}