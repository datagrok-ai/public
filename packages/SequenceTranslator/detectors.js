/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */

class SequenceTranslatorPackageDetectors extends DG.Package {
  //name: autostart
  //tags: autostart
  //description: SequenceTranslator bootstrap
  autostart() {
    this.logger.debug('ST: detectors.js: autostart()');

    this.autostartContextMenu();
  }

  autostartContextMenu() {
    const logPrefix = `ST: detectors.js: autostartContextMenu()`;
    this.logger.debug(`${logPrefix}, start`);
    grok.events.onContextMenu.subscribe((event) => {
      this.logger.debug(`${logPrefix}, onContextMenu, start`);
      const item = event.args.item;
      if (item) {
        if (item && (
          (item instanceof DG.GridCell || item.constructor.name === 'GridCell')
        )) {
          const packageName = this.name;
          grok.functions.call(`${packageName}:addContextMenu`, {event: event})
            .catch((err) => {
              this.logger.error('addContextMenu must not throw any exception');
              this.logger.error(err.toString());
            });
        }
      }
      this.logger.debug(`${logPrefix}, onContextMenu, end`);
    });
  }
}
