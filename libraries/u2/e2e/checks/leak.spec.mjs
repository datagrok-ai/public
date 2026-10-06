/* P4.5 — the host's spec is a constant. `u2DesignerApp` hands the same `DESIGNER_SPEC` object to
   every designer of the session, and the editor patches its document in place (DD10): a designer
   that adopted it left its edits in the next one. Reopens the app the way a user does — closing
   every view, which is why this file runs last. */
import {ok, shot} from '../local.mjs';
import {dropControl, named, nodeCount, reopenApp, statusText} from '../lib.mjs';

export async function fixture(page) {
  await reopenApp(page);
}

async function checkSpecNotAdopted(page) {
  const fresh = await nodeCount(page);
  ok('leak/1a/the-app-opens-on-the-spec-its-host-wrote', fresh === 13 &&
    !(await named(page)).includes('textInput1'), `${fresh} nodes ${JSON.stringify(await named(page))}`);

  await dropControl(page, 'u2-text-input', 'text-input', 'details');
  const edited = await nodeCount(page);
  ok('leak/1b/an-insert-lands-in-this-designers-document', edited === fresh + 1 &&
    (await named(page)).includes('textInput1'), `${edited} nodes, status="${await statusText(page)}"`);

  await reopenApp(page);
  const again = await nodeCount(page);
  await shot(page, 'leak-1-reopened');
  ok('leak/1c/the-next-designer-of-the-session-opens-on-the-same-13-nodes', again === fresh &&
    !(await named(page)).includes('textInput1'),
  `${again} nodes ${JSON.stringify(await named(page))}`);
}

export const checks = [
  {id: 'leak/1 the host spec is a constant: an edit never reaches the next designer of the session',
    run: checkSpecNotAdopted},
];
