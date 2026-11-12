import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Widgets: containers', () => {
  test('empty takes out text and comments, and kills the elements', async () => {
    const host = ui.div();
    const child = ui.div();
    let killed = false;
    DG.Widget.registerCleanup(child, () => killed = true);
    host.append('plain text', child, document.createComment('note'));
    expect(ui.empty(host), host);
    expect(host.childNodes.length, 0);
    expect(killed, true);
  });

  test('divV and divH take null for no children', async () => {
    const v = ui.divV(null as unknown as HTMLElement[]);
    const h = ui.divH(null as unknown as HTMLElement[]);
    expect(v.childElementCount, 0);
    expect(h.childElementCount, 0);
    expect(v.classList.contains('d4-flex-col'), true);
    expect(h.classList.contains('d4-flex-row'), true);
  });
});
