/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class should follow the `<PackageName>Detectors` notation, otherwise it won't be loaded.
 */
class HelmPackageDetectors extends DG.Package {

  // -- autostart --

  //name: autostart
  //tags: autostart
  //description: Helm bootstrap
  autostart() {
    this.logger.debug('Helm: detectors.js: autostart()');
    this.autostartContextMenu();
  }

  autostartContextMenu() {
    grok.events.onContextMenu.subscribe((event) => {
      if (!event.args.item) return false;
      const item = event.args.item;
      if (item instanceof DG.FileInfo || item.constructor.name === 'FileInfo') {
        const fi = item;
        if (fi.path.endsWith('.json')) {
          grok.functions.call('Helm:addContextMenuForFileInfoJson', {event})
            .catch((err) => {
              this.logger.error('addContextMenu... must not throw any exception', err.toString());
              this.logger.error(error);
            });
          event.preventDefault();
          return true;
        }
      }
    });
  }
}
