

@grok.decorators.viewer({
  name: '#{NAME}',
  description: 'Creates #{NAME} viewer'
})
export class #{NAME}Viewer extends DG.JsViewer {
  constructor() {
    super();
  }

  // Additional chart settings
  init() {
  }
  
  // Stream subscriptions
  onTableAttached() {
    this.init();
    this.render();
  }

  // Cancel subscriptions when the viewer is detached
  detach() {
  }
  
  // Override to handle property changes
  onPropertyChanged(property : DG.Property | null) {
    super.onPropertyChanged(property);
    this.render();
  }
  
  render(computeData = true) {
  }
}
