/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class should follow the `<PackageName>Detectors` notation, otherwise it won't be loaded.
 */
class HelmPackageDetectors extends DG.Package {

  //name: autostart
  //tags: autostart
  //description: Helm bootstrap
  autostart() {
    this.logger.debug('Helm: detectors.js: autostart()');

    this.autostartInputs();
  }

  autostartInputs() {
    const logPrefix = 'Helm: detectors.js: autostartInputs()';

    if (!!ui.input.helmAsync)
      this.logger.warning(`${logPrefix}, ui.input.helmAsync is already defined`);
    else
      ui.input.helmAsync = this.helmAsync;
  }

  //name: helmAsync
  //input: string name
  //input: object options
  //output: object result
  helmAsync(name, options) {
    return grok.functions.call('Helm:helmInput', {name, options});
  }
}
