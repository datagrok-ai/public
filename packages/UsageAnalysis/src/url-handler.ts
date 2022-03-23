import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class UrlHandler {
  public constructor() {
    let arr = ['Events', 'Errors', 'Users', 'Overview', 'Data'];
    //
    // grok.events.onEvent('d4-current-view-changed').subscribe(
    //     () => {
    //       grok.shell.info(grok.shell.v.name);
    //       arr.indexOf(grok.shell.v.name);
    //       grok.shell.v.path = `/apps/UsageAnalysis/${grok.shell.v.name}`;
    //     });
  }
}