import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {ChoiceInput, DateInput, Scope} from '@datagrok-libraries/u2';
import {columnInput, columnsInput, fileInput, filesInput, groupPick, moleculeInput, tableInput, userPick}
  from '@datagrok-libraries/u2/dg';

interface Twin {root: HTMLElement; dispose(): void}

interface Pair {
  name: string;
  /** The platform input; or, for one JS cannot build, the affordances xamgle gives it. */
  platform: (() => DG.InputBase) | string[];
  /** Built under the test's `owner`; one that awaits first builds there itself. */
  twin: (owner: Scope) => Promise<Twin> | Twin;
  /** What the platform is known to offer, so that a platform input that stops naming its icons fails here
   * instead of making every twin pass. */
  floor?: string[];
}

function shown(root: HTMLElement, selector: string): string[] {
  return [...root.querySelectorAll(selector)]
    .filter((el) => el.getClientRects().length > 0)
    .map((el) => el.getAttribute('aria-label')!);
}

/** The affordances shown while `root` is in the page (a control hidden by its own styles has no box). */
function affordances(root: HTMLElement, selector: string): string[] {
  document.body.append(root);
  try {
    return shown(root, selector);
  } finally {
    root.remove();
  }
}

const demog = () => grok.data.demo.demog(10);

/** Every platform input `ui.input.X` with a u2/dg twin. The platform adds some affordances outside the input itself
 * (xamgle `inputs.dart` `_initInputCustomizations` hooks every Dart TableInput, FileInput and FilesInput), which u2
 * inputs never get, so a port that forgets one drops the feature silently (SEAM S-20). An affordance is a displayed
 * control with an accessible name: the platform's icons carry their tooltip as `aria-label`, u2's icon buttons too. */
const PAIRS: Pair[] = [
  {name: 'file', platform: () => ui.input.file('File'), twin: () => fileInput('File'),
    floor: ['Open file', 'Add file from Files']},
  {name: 'files', platform: () => ui.input.files('Files'), twin: () => filesInput('Files'),
    floor: ['Open files', 'Add file from Files']},
  {name: 'folder', platform: () => ui.input.folder('Folder'), twin: () => fileInput('Folder', {mode: 'folder'}),
    floor: ['Add file from Files']},
  // the platform builds a FileInput on a connection only in Dart (`data_connection_meta.dart`)
  {name: 'file on a connection', platform: ['Choose directory or file'], twin: async (owner) => {
    const connection = await grok.dapi.connections.filter('shortName = "DemoFiles"').first();
    return owner.build(() => fileInput('File', {connection}));
  }},
  {name: 'table', platform: () => ui.input.table('Table'), twin: () => tableInput('Table'),
    floor: ['Open file', 'Add file from Files', 'Query database']},
  {name: 'column', platform: () => ui.input.column('Column', {table: demog()}),
    twin: () => columnInput('Column', demog())},
  {name: 'columns', platform: () => ui.input.columns('Columns', {table: demog()}),
    twin: () => columnsInput('Columns', demog())},
  {name: 'user', platform: () => ui.input.user('User'), twin: () => userPick({label: 'User'})},
  {name: 'group', platform: () => ui.input.userGroups('Groups'), twin: () => groupPick({label: 'Groups'})},
  {name: 'molecule', platform: () => ui.input.molecule('Molecule'), twin: () => moleculeInput('Molecule')},
  {name: 'choice', platform: () => ui.input.choice('Choice', {items: ['a', 'b'], value: 'a'}),
    twin: () => new ChoiceInput({label: 'Choice', items: ['a', 'b'], value: 'a'})},
  {name: 'date', platform: () => ui.input.date('Date'), twin: () => new DateInput({label: 'Date'}),
    floor: ['Pick a date']},
];

category('u2/dg: affordance parity', () => {
  for (const pair of PAIRS) {
    test(`${pair.name}: the u2/dg twin offers what the platform input offers`, async () => {
      const platform = Array.isArray(pair.platform) ? pair.platform :
        affordances(pair.platform().root, '.grok-icon[aria-label]');
      for (const name of pair.floor ?? [])
        expect(platform.includes(name), true, `the platform ${pair.name} input lost "${name}" (${platform})`);
      const owner = Scope.root('affordance parity');
      try {
        const twin = await owner.build(() => pair.twin(owner));
        const offered = affordances(twin.root, '[aria-label]');
        const missing = platform.filter((name) => !offered.includes(name));
        expect(missing.length === 0, true,
          `the u2/dg ${pair.name} input lacks ${missing.map((m) => `"${m}"`).join(', ')} (it offers ${offered})`);
      } finally {
        owner.dispose();
      }
    });
  }
});
