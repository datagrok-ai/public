/** What holds a rule: `lint:u2/<rule>` (eslint-plugin-u2, run by `check-consumer` and `grok check`), `usage:<code>`
 * (the usage error the engine raises), `probe:<rule>` (the operability probe), `type:<Type.member>` (a compile error),
 * `ratchet:<count>` (`consumers-baseline.json`), `api:<rule>` and `choke:<id>` (u2's own declarations and sources)
 * and `debug:<test>` (a warning with u2 debug on, pinned by `tests/<test>.test.js`). */
export type U2Guard = `lint:u2/${string}` | `usage:${string}` | `probe:${string}` | `type:${string}` |
  `ratchet:${string}` | `api:${string}` | `choke:${string}` | `debug:${string}`;

/** What a rule says about a kind; the manifest's `kinds` carry it. */
export interface U2KindNote {
  /** One line: when this kind is the one to use rather than a neighbour. */
  whenToUse?: string;
  /** One line: a job this kind looks fit for, and what does it instead. */
  notFor?: string;
}

/** One of u2's usage rules: the README's "Rules that matter", the digest in `AGENTS.md` and the manifest's kind notes
 * are generated from these (`npm run manifest`). */
export interface U2Rule {
  /** Unique, kebab-case. */
  id: string;
  /** The rule, one sentence. */
  rule: string;
  /** Why it matters, one sentence. */
  why: string;
  /** The README's text under the rule: how, naming the API that makes it easy. */
  detail: string;
  /** One line of code that breaks it. */
  wrong: string;
  /** The same job done right, one line. */
  right: string;
  /** What catches a slip; none where nothing does yet. */
  guards: U2Guard[];
  /** Kind → what the rule says about it. */
  kinds?: Record<string, U2KindNote>;
}

const RULES: readonly U2Rule[] = [
  {
    id: 'entry-points',
    rule: 'Import from the entry points only.',
    why: 'A path inside the package is no API: its files move and change between releases.',
    detail: '`@datagrok-libraries/u2` (platform-free), `@datagrok-libraries/u2/dg` (the platform layer), ' +
      '`/filter` and `/dg/filter` (the filter model and its platform side) and `/testing` (the testing kit); ' +
      'never a path inside the package (`…/u2/src/…`, `…/u2/dist/…`). There is no CSS to import: the first u2 ' +
      'control a page loads installs u2\'s stylesheet.',
    wrong: 'import {TextInput} from \'@datagrok-libraries/u2/dist/src/index.js\';',
    right: 'import {TextInput} from \'@datagrok-libraries/u2\';',
    guards: ['lint:u2/import-entry', 'ratchet:deepImports', 'ratchet:deepDgIndex'],
  },
  {
    id: 'owner-function',
    rule: 'Build inside an owner, as a function.',
    why: 'What an owner builds is disposed with it; what nothing owns outlives the screen that showed it.',
    detail: 'Element factories and `each()` need an owner building. A piece of UI with a life of its own (shown ' +
      'and closed, rebuilt, in a dialog, addressed by tests) is a function component, `component((self) => { …; ' +
      'return {root, …api}; }, {name})`: its builder owns everything it makes, and all of it is disposed with the ' +
      'control; `self` has the owner verbs (`self.effect`, `self.subscribe`, `self.setInterval`, `self.task`). A ' +
      'piece that is only part of its parent is a plain function that builds it, called inside the parent\'s ' +
      'builder — it is the parent\'s. After an `await` nothing is building: `self.task(async (t) => … ' +
      't.build(…))`. A tree that nothing else owns (an app\'s home that outlives its views) starts on purpose as ' +
      '`Scope.root(\'<what>\')`, and whoever starts it disposes it. Consumer code does not subclass u2\'s ' +
      'controls. A control\'s own members (`name`, `root`, `scope`, `visible`, `effect`, `own`, `on`, `dispose`, ' +
      '…) are not names for a part: name it for what it holds (`nameInput`); u2 reports the collision ' +
      '(`reserved-key`). The first sample: docs/samples/first-control.ts.',
    wrong: 'const nameInput = new TextInput({label: \'Name\'});  // at module load or after an await: nothing ' +
      'owns it',
    right: 'component((self) => { const nameInput = new TextInput({label: \'Name\'}); return {root: [nameInput], ' +
      'nameInput}; })',
    guards: ['lint:u2/build-in-owner', 'lint:u2/no-ui-in-class', 'lint:u2/no-subclass', 'lint:u2/no-unowned-root',
      'usage:owner-missing', 'usage:owner-disposed', 'usage:reserved-key', 'ratchet:subclassU2',
      'ratchet:classBuildsUi', 'ratchet:newScope', 'ratchet:runWithUndefined'],
  },
  {
    id: 'input-value',
    rule: 'An input\'s value is `input.value`.',
    why: 'A slot u2 keeps live needs the signal: the plain value is read once there and never changes.',
    detail: 'Read inside an effect or a `computed` it is followed, and assigning it writes it; `input.peek()` reads ' +
      'it untracked, and `input.signal` is the signal for wiring (`bind: input.signal`, `span(input.signal)` — ' +
      '`span(input.value)` shows a snapshot). A validator may read other inputs: it re-judges whenever they ' +
      'change.',
    wrong: 'span(nameInput.value)',
    right: 'span(nameInput.signal)',
    guards: ['lint:u2/no-snapshot-in-live-slot', 'api:value-signal', 'ratchet:valueValue', 'ratchet:valuePeek',
      'usage:bind-readonly'],
  },
  {
    id: 'register-where-it-lives',
    rule: 'Register on the thing that lives as long.',
    why: 'An effect put on a long-lived control from a row is added again by every rebuild of the row.',
    detail: 'Inside an `each()` row or an effect run, put effects and cleanups on a control built in that row ' +
      '(`input.effect(...)`), not on the long-lived control around it — that one would gain one more per ' +
      'rebuild. With u2 debug on (`localStorage[\'u2.debug\'] = \'1\'`, and in every bdd run with the leak ' +
      'audit) u2 warns when the same call site registers on a long-lived control from two different rows or ' +
      'runs.',
    wrong: 'self.effect(() => paint(row));  // inside an each() row',
    right: 'rowInput.effect(() => paint(row));  // on a control built in that row',
    guards: ['debug:rebuilt-scope-warning'],
  },
  {
    id: 'on-input',
    rule: 'Act on the user with `onInput`.',
    why: '`onChanged` also hears writes from code, so the work runs again whenever code sets the value.',
    detail: '`onChanged` also hears writes from code (it matches the platform\'s `onChanged`). "The user picked a ' +
      'library, load it" is `onInput`.',
    wrong: 'new ChoiceInput({label: \'Library\', items, onChanged: (name) => load(name)})',
    right: 'new ChoiceInput({label: \'Library\', items, onInput: (name) => load(name)})',
    guards: [],
  },
  {
    id: 'reset-follow',
    rule: 'Load a record with `reset(value)`; follow a model with `follow(read)`.',
    why: 'A plain write fires `onChanged` and keeps the touched state of the record before, so a loaded record ' +
      'looks edited.',
    detail: 'Neither fires `onChanged`; `reset` also clears the touched state and any unreadable text. ' +
      '`Input.system()` is the input\'s own pruning, and it is *not* silent.',
    wrong: 'nameInput.value = record.name;  // loading a record',
    right: 'nameInput.reset(record.name);',
    guards: [],
  },
  {
    id: 'two-verdicts',
    rule: 'Save goes by `validity`; what a field shows goes by `error`.',
    why: '`error` stays empty until a field is edited, left or validated, so a Save that follows it is on for an ' +
      'empty form.',
    detail: 'An input has two verdicts. `validity` judges the value from the start (what Save goes by — disable ' +
      'it and say why in its tooltip); `error` is what the input shows, only once it was edited, left or ' +
      'validated (what a tab marker or a status line goes by). A `Form` has both.',
    wrong: 'button(\'Save\', save, {disabled: nameInput.error})',
    right: 'button(\'Save\', save, {disabled: form.validity})',
    guards: [],
  },
  {
    id: 'tabs-mark',
    rule: 'A tab marks itself: give it a `validity`, or make its content a `Form`.',
    why: 'A marker set by hand drifts from the fields it stands for.',
    detail: 'Before reaching for tabs in a narrow panel, consider showing what matters and folding the rest into a ' +
      '`Section` with a `hint` (what it holds), a `validity` (an error is never hidden by the fold) and a ' +
      '`persistKey` (the user\'s open or close is remembered).',
    wrong: 'tabs.setTabTitle(\'fields\', \'Fields (!)\');  // marked by hand',
    right: '{id: \'fields\', title: \'Fields\', content: () => fieldsForm, validity: fieldsForm.validity}',
    guards: [],
  },
  {
    id: 'say-why-disabled',
    rule: 'Say why something is disabled.',
    why: 'A control that is off with no reason leaves the user guessing what to do.',
    detail: 'Every off takes the reason the user reads, a string or a signal of one (`null` while it is on): ' +
      '`disabled: \'Pick a file first\'` on an input, a button, an action, a menu item, a ribbon button; ' +
      '`dialog.okDisabled(reason)`. There is no `enabled` and no `disabled: true`. The reason is the tooltip; a ' +
      'dialog\'s OK also names the form\'s first failing field. A disabled u2 button is `aria-disabled` and keeps ' +
      'its tab stop, so keyboard focus shows the reason too; a press on it does nothing.',
    wrong: 'button(\'Save\', save, {disabled: true})',
    right: 'button(\'Save\', save, {disabled: \'Pick a file first\'})',
    guards: ['type:ButtonOptions.disabled', 'type:InputOptions.disabled', 'api:no-enabled', 'api:disabled-reason',
      'lint:u2/reason-text', 'probe:disabled-without-reason', 'ratchet:enabledOpt', 'ratchet:okEnabled',
      'usage:reason-boolean', 'usage:reason-missing'],
  },
  {
    id: 'run-returns',
    rule: 'An action\'s `run` returns a promise or nothing.',
    why: 'u2 waits for what `run` answers: a value handed back by mistake, such as a `DG.Row`, is no work to ' +
      'wait for.',
    detail: '`run: () => save()` hands back what `save` answers; a value handed back by mistake (a `DG.Row`) is a ' +
      'type error, so wrap a call that answers one in braces. A `run` or a click handler that throws or rejects ' +
      'goes to its owner\'s `onError`, then up its owners, then to `Scope.onUnhandledError`, which logs it to the ' +
      'console: the user sees nothing. Say what failed where the work is (`notify.error(…)`), or once for a whole ' +
      'component (`self.onError = (e) => notify.error(…)`); `confirmDialog` and a dialog\'s `onOK` keep their ' +
      'own failure in the dialog.',
    wrong: '{text: \'Add row\', run: () => df.rows.addNew()}',
    right: '{text: \'Add row\', run: () => { df.rows.addNew(); }}',
    guards: ['type:Action.run'],
  },
  {
    id: 'label-column',
    rule: 'Forms that must line up share one `LabelColumn`.',
    why: 'Each form aligns its own captions, so two forms one above the other start their editors at different ' +
      'places.',
    detail: 'The tabs of one editor, a header over a body: each `Form` takes the same `labelColumn`; `labelWidth` ' +
      'fixes a column (one-letter captions).',
    wrong: '[new Form(), new Form()]  // two tabs of one editor',
    right: 'const labelColumn = new LabelColumn(); const tabs = [new Form({labelColumn}), ' +
      'new Form({labelColumn})];',
    guards: [],
  },
  {
    id: 'widths-options',
    rule: 'Widths are options, not CSS.',
    why: 'A width rule over u2\'s classes fights u2\'s own sizing (hence the `!important`) and breaks when they ' +
      'change.',
    detail: '`fill: true` stretches an editor across its row; `width` (px or a CSS length) sets it outright. No ' +
      '`!important` needed. A form is 700px at most: `new Form({maxWidth: 820})` sets its own. Several inputs ' +
      'under one caption (Colors: Line, Background, Text) are `form.addGroup(caption, inputs)`. A form\'s button ' +
      'row is on the right (`form.addButtons(build, {align: \'start\' | \'editors\'})` moves it). Spacing too: ' +
      '`divH([cancel, save], {gap: \'m\'})` on the token scale (`xs` 2px … `xxl` 24px), never `style.gap`.',
    wrong: '.my-panel .u2-input-editor { width: 100% !important; }',
    right: 'new TextInput({label: \'Name\', fill: true})',
    guards: ['ratchet:cssImportant', 'ratchet:cssInternalLines'],
  },
  {
    id: 'page-layout',
    rule: 'Lay a view out with `Page`; a list beside its details is a `MasterDetail`, past a few hundred items a ' +
      '`VirtualList`.',
    why: 'A layout written in package CSS breaks as the view narrows, and drifts from every other view\'s.',
    detail: '`new Page({header, sections, footer, maxWidth})`: the header stays put, the body alone scrolls and the ' +
      'action row stays pinned under it, its buttons on the right and a `statusLine` on their left; the forms ' +
      'on the page share one label column. A section is a heading with its actions over its content ' +
      '(`{title, actions, content, maxWidth}`). `new MasterDetail({items, keyOf, row, detail, disabled})` puts ' +
      'the rows — a check, a hint — beside the current one\'s detail, built when first shown and kept, off with ' +
      'its reason; the arrow keys move the current row and Space checks it. It renders every row, so it is for ' +
      'lists up to a few hundred items; a longer one is a `VirtualList` (`new VirtualList({keyOf, render})`, ' +
      '`setItems(items)`, its `selectedIndex` the current row) beside a detail of your own, in a `Splitter`. ' +
      'Spacing and heights are the `--u2-page-*` and `--u2-master-detail-height` properties, set on an element ' +
      'of your own.',
    wrong: '.my-editor { display: flex; flex-direction: column; height: 100%; }  /* and a pinned footer by hand */',
    right: 'new Page({content: form, footer: [statusLine(form.shownProblem), saveButton], maxWidth: 820})',
    guards: ['ratchet:cssLines', 'ratchet:cssInternalLines'],
    kinds: {
      'page': {whenToUse: 'the content of a view or a dialog-sized editor: a header, a body that alone scrolls, ' +
        'a pinned action row, sections with their actions'},
      'master-detail': {whenToUse: 'a short list (up to a few hundred items: it renders every row) beside the ' +
        'detail of its current item, each item with a check and a hint, its detail kept while it is not shown',
      notFor: 'many records to sort and group: `CatalogTable`, its details in the context panel; a longer list ' +
        'beside a detail: `VirtualList` and a detail of your own in a `Splitter`'},
      'splitter': {notFor: 'a list beside the details of its current item: `MasterDetail`'},
    },
  },
  {
    id: 'own-styling',
    rule: 'Style your own classes, never u2\'s.',
    why: 'u2\'s classes and attributes are its internals: they change without notice, and a rule over them ' +
      'fights u2\'s own.',
    detail: 'A package\'s CSS styles its own classes and sets u2\'s `--u2-*` custom properties on its own ' +
      'elements — the ones each kind lists as `styling` in manifest.json (README › Styling); it never targets ' +
      '`.u2-*` classes or `[data-u2…]` attributes, and code never looks u2\'s parts up by class name. Tests reach ' +
      'a control through the testing kit (`byName`, `byKind`, `state`).',
    wrong: '.my-app .u2-table-row:hover { background: #eef; }',
    right: '.my-app { --u2-row-hover-bg: var(--blue-1); }',
    guards: ['lint:u2/no-internal-selectors', 'ratchet:tsInternalSelectors', 'ratchet:cssInternalLines'],
  },
  {
    id: 'element-value',
    rule: 'A value shown as an element is a `DynamicInput` row; the rail takes buttons.',
    why: 'The rail beside an editor is one row high: a drawing put there spills over the rows under it.',
    detail: 'Inside a form, a value drawn as an element (a structure drawing, a chart) is a `DynamicInput`: a row ' +
      'of the form, in its label column, that grows to the element. The rail beside an editor takes buttons and ' +
      'links (`addOptions(iconButton(…))`). Read-only details in a pane of their own are a `tableFromMap`: muted ' +
      'keys, wrapping values (docs/recipes/forms-and-details.md).',
    wrong: 'smiles.addOptions(grok.chem.drawMolecule(value, 160, 100))',
    right: 'new DynamicInput({label: \'Structure\', value, render: (s) => s === null ? null : ' +
      'grok.chem.drawMolecule(s, 160, 100)})',
    guards: ['type:Input.addOptions', 'probe:oversized-option'],
    kinds: {
      'dynamic-input': {
        whenToUse: 'a value shown as an element (a structure drawing, a chart): its row grows to fit it, where an ' +
          'options-rail child would spill',
        notFor: 'read-only details in a pane of their own: `tableFromMap(map)`',
      },
      'table': {
        whenToUse: 'read-only details of an item in a pane of their own: `tableFromMap(map)`, muted keys and ' +
          'wrapping values',
        notFor: 'a value drawn as an element among a form\'s inputs: `DynamicInput`',
      },
    },
  },
  {
    id: 'listen-follow',
    rule: 'Listen with `on(el, type, handler)`, follow with `effect(fn)`, on the elements a row builds.',
    why: 'Registered on the long-lived control instead, the handlers and the element stay with it after the row ' +
      'is gone.',
    detail: 'Both go with the row or the card that built the element. `self.effect` on the component around the ' +
      'row would tie them, and the element, to the long-lived control.',
    wrong: 'self.listen(card, \'click\', open);  // inside an each() row',
    right: 'on(card, \'click\', open);',
    guards: ['choke:handlers', 'debug:rebuilt-scope-warning'],
  },
  {
    id: 'docked-editor',
    rule: 'A docked editor panel is laid out by its form: columns, a status line and folds.',
    why: 'A panel laid out by package CSS breaks when the dock narrows it.',
    detail: '`form.addColumns(inputs, {minWidth})` puts short fields two to a row, each editor filling its ' +
      'column, one to a row when narrow. `statusLine(problem, otherwise)` beside Save says what it does or why it ' +
      'cannot yet; a flush `Accordion` pane with a `hint` folds the rest; a `ChoiceChip` (each item\'s ' +
      '`description` says what it stands for) changes a small choice in two clicks. The recipe: ' +
      'docs/recipes/forms-and-details.md. `new Form({quiet: true})` drops the underlines at rest: for an ' +
      'inspector, not a form people fill in.',
    wrong: '.my-editor .u2-form { display: grid; grid-template-columns: 1fr 1fr; }',
    right: 'form.addColumns([symbol, name], {minWidth: 160})',
    guards: [],
    kinds: {
      'choice-chip': {
        whenToUse: 'a small choice inline among others, its label inside it ("R1 OH ▾", "Sort Name ▾"): a menu of ' +
          'the choices, each with what it stands for, the current one checked',
        notFor: 'a field of a form, in its label column: `ChoiceInput`',
      },
      'choice-input': {notFor: 'a small choice inline among others, its label inside it: `ChoiceChip`'},
      'dropdown-button': {notFor: 'a choice that keeps a current value: `ChoiceChip`'},
    },
  },
  {
    id: 'dock',
    rule: 'Dock with `ViewHost.of(view).dock(control, place)`.',
    why: 'Docked by hand, a control outlives its panel and its view.',
    detail: 'Closing the panel — its ✕, or `dockManager.close` from code — disposes the control, and so does the ' +
      'view\'s close; disposing the control closes the panel. Chrome that is not a view (a context-panel widget, ' +
      'a pane) takes `host(control)`, which the platform\'s close of that chrome disposes.',
    wrong: 'grok.shell.dockManager.dock(panel.root, \'right\')',
    right: 'ViewHost.of(view).dock(panel, \'right\', {title: \'Details\'})',
    guards: [],
  },
  {
    id: 'find-from-root',
    rule: 'Find a view\'s control from its root.',
    why: 'The platform hands out a new wrapper object for a view each time, so nothing stored on one is there the ' +
      'next time.',
    detail: 'The platform hands out a new wrapper object for a view each time (`grok.shell.views`), so nothing ' +
      'stored on one is there the next time; reach the content with ' +
      '`Control.forElement(view.root.querySelector(\'.my-view\'))`.',
    wrong: '(view as any).editor = editor;  // read back from grok.shell.v later',
    right: 'Control.forElement(view.root.querySelector(\'.my-view\'))',
    guards: [],
  },
  {
    id: 'dialog-cancel',
    rule: 'A dialog comes with CANCEL; one whose work applies as it is done is `closeOnly`.',
    why: 'A dialog that applies its work as it goes has nothing to cancel: one CLOSE says so, and ✕ and Escape ' +
      'agree.',
    detail: 'The first `onOK` / `addButton` adds CANCEL, as in every platform dialog. A dialog whose work applies ' +
      'as it is done (a tracker, a viewer of progress) is `new Dialog(title, {closeOnly: true})`: one CLOSE, ' +
      'never a CANCEL; ✕ and Escape close it the same way.',
    wrong: 'dialog.getButton(\'CANCEL\')?.remove();',
    right: 'new Dialog(\'Progress\', {closeOnly: true})',
    guards: ['usage:close-only-cancel'],
  },
  {
    id: 'settings-scope',
    rule: 'Personal settings apply at once; shared and admin settings stage behind SAVE N CHANGES.',
    why: 'The user\'s own preference should take effect as it is set, while a change other people feel waits until ' +
      'it is saved on purpose.',
    detail: 'A form of the user\'s own settings (`grok.userSettings`, a view\'s look) writes and keeps each valid ' +
      'edit as it is made, with no Save button: `propertyForm(props, settings, {onChanged: save})` ' +
      '(docs/samples/property-form.ts). Settings shared with a team or set by an admin keep the edits in the form ' +
      'and save them together with one button that counts them ("SAVE 3 CHANGES"), off with its reason while ' +
      'there is nothing to save or a field is invalid.',
    wrong: 'form.addButtons((row) => row.append(button(\'Save\', save)))  // the user\'s own settings',
    right: 'propertyForm(props, settings, {onChanged: () => save(settings)})  // the user\'s own settings',
    guards: [],
  },
  {
    id: 'dialog-outlives-row',
    rule: 'A dialog opened from a row outlives the row.',
    why: 'A rebuild of the row does not close the dialog, so content the row built would lose its handlers under ' +
      'it.',
    detail: 'Shown, a dialog or a `confirmDialog` leaves the row or effect run that built it, so a rebuild does ' +
      'not close it; the content\'s own handlers and live text stay with whoever built them, so build that ' +
      'content with `dialog.add(() => …)`. A `keepAlive` dialog belongs to its builder: make it in the control ' +
      'that reopens it, not in a row.',
    wrong: 'dialog.add(detailsOf(row));  // built by the row',
    right: 'dialog.add(() => detailsOf(row));',
    guards: [],
  },
  {
    id: 'balloon-content',
    rule: 'A balloon\'s content is the balloon\'s.',
    why: 'Content built in a row\'s handler keeps its handlers on the row, so its Undo dies with a rebuild.',
    detail: '`notify.info(() => div([\'Deleted. \', link(\'Undo\', undo)]), {onClose})` builds the content owned ' +
      'by the balloon and ends it when the balloon goes; an element built in a row\'s handler would keep its ' +
      'handlers on the row.',
    wrong: 'notify.info(div([\'Deleted. \', link(\'Undo\', undo)]))',
    right: 'notify.info(() => div([\'Deleted. \', link(\'Undo\', undo)]), {onClose})',
    guards: [],
  },
  {
    id: 'app-header',
    rule: 'An app home starts with `AppHeader`.',
    why: 'The same header in every app tells the user what the app does; one built by hand drifts from the rest.',
    detail: 'The icon, what the app does as bullets (text, an element or a link each) and Learn more: the same ' +
      'header in every app, built in one call (docs/recipes/shell-integration.md). It follows its own width, so ' +
      'build it under the view\'s host — in the `ViewHost.create` builder, or in `ViewHost.of(view).build(…)` on ' +
      'a view you made — and it goes with the view.',
    wrong: 'div([h1(\'My App\'), span(\'Design molecules with your team\')])',
    right: 'ViewHost.create({name: \'My App\', content: () => [header(), campaigns()]})  // header() builds the ' +
      'AppHeader',
    guards: [],
    kinds: {'app-header': {whenToUse: 'the top of an app\'s home: its icon, what it does and Learn more'}},
  },
  {
    id: 'owner-subscribe',
    rule: 'Own a platform subscription with the owner\'s `subscribe`.',
    why: 'A raw subscription is never ended: its handler keeps running after the view or the control is gone.',
    detail: '`self.subscribe(grok.events.onViewRemoved, handler)` — or `host.subscribe(…)` on a view\'s host — ' +
      'ends with that owner, runs the handler batched and contained, and shows in the leak audit. The free ' +
      '`subscribe` from u2 takes whoever is building: an effect run or an `each()` row, so a subscription made ' +
      'there goes when the run or the row does. After an `await` nothing is building: subscribe on the owner, ' +
      'and start timers on it (`self.setInterval`).',
    wrong: 'grok.events.onViewRemoved.subscribe(handler)',
    right: 'self.subscribe(grok.events.onViewRemoved, handler)',
    guards: ['lint:u2/no-raw-subscribe', 'usage:source-missing', 'api:owner-positional', 'ratchet:rawSubscribe',
      'ratchet:freeSubscribe'],
  },
  {
    id: 'view-host',
    rule: 'Code living in a platform view is owned by its `ViewHost`.',
    why: 'The view is what the user closes: what its host owns ends then, whichever wrapper of the view the code ' +
      'held.',
    detail: '`ViewHost.of(tv)` is one host per view, whichever wrapper you hold, disposed when the shell removes ' +
      'the view: build under it (`host.build(…)`; after an `await`, `host.task(async (t) => … t.build(…))`), ' +
      'subscribe with `host.subscribe(stream, handler)`, read `host.signals` (the current row, selection, ' +
      'columns and the `grid` as signals; an effect over `grid` re-subscribes per grid, no timers), put the ' +
      'ribbon group up with `host.ribbon({actions, buttons})`, show its dialogs with `host.show(dialog)` and dock ' +
      'with `host.dock`. A new view is `ViewHost.create({name, content})` (docs/recipes/shell-integration.md).',
    wrong: 'grok.events.onViewRemoved.subscribe((v) => { if (v.name === tv.name) panel.dispose(); })',
    right: 'ViewHost.of(tv).build(() => detailsPanel())',
    guards: ['usage:view-not-ready', 'ratchet:viewZoo'],
    kinds: {'ribbon-group': {whenToUse: 'a package\'s group on a platform view\'s ribbon: ' +
      '`ViewHost.of(view).ribbon({actions, buttons})`'}},
  },
  {
    id: 'confirm',
    rule: 'Ask with `confirmDialog()` before what cannot be undone.',
    why: 'The platform asks the same way: CANCEL first and the action named as a verb, busy until the work is ' +
      'done.',
    detail: 'CANCEL and the action as a verb; `onConfirm` keeps the dialog busy until the work is done.',
    wrong: 'if (window.confirm(\'Delete?\')) remove();',
    right: 'confirmDialog({title: \'Delete campaign\', message: \'Delete KHL-3?\', confirmText: \'Delete\', ' +
      'onConfirm: () => remove()})',
    guards: [],
  },
  {
    id: 'records-checklist',
    rule: 'Rows the user edits are a `RecordsInput`; a pick with an order is a `ChecklistInput`.',
    why: 'Both are inputs: they validate, mark their form and read back as one value.',
    detail: 'A template\'s fields or stages, the functions to run in order — never an `ItemsGrid` or a hidden ' +
      'columns input (docs/recipes/forms-and-details.md).',
    wrong: 'new MultiChoiceInput({label: \'Run\', items: functionNames})  // the order is lost',
    right: 'new ChecklistInput({label: \'Run\', items: functionNames, orderable: true})',
    guards: [],
    kinds: {
      'records-input': {whenToUse: 'rows the user edits in place: a template\'s fields or stages'},
      'checklist-input': {whenToUse: 'a pick where the order matters too: the functions to run, in order'},
      'multi-choice-input': {notFor: 'a pick where the order matters too: `ChecklistInput`'},
    },
  },
  {
    id: 'catalog-table',
    rule: 'A list of records is a `CatalogTable`.',
    why: 'It owns the order, the groups and the chosen columns, and keeps them where the app keeps them.',
    detail: 'Columns (`value` to sort and group by, `render` to draw), groupings, row `actions` and a `persist` ' +
      'adapter; it owns the order, the groups and the chosen columns, and `chrome()` gives the heading its ' +
      'columns, group and refresh icons. Never a `ui.table` with header listeners and injected rows ' +
      '(docs/recipes/catalog-table.md).',
    wrong: 'ui.table(campaigns, (c) => [c.name, c.status])  // header listeners, rows put in by hand',
    right: 'new CatalogTable({items, keyOf: (c) => c.code, columns, groupings, persist: \'my-app.campaigns\'})',
    guards: [],
    kinds: {
      'catalog-table': {whenToUse: 'a list of records to sort, group, choose the columns of and act on'},
      'table': {notFor: 'records to sort, group and act on: `CatalogTable`'},
    },
  },
  {
    id: 'column-choice',
    rule: 'A column choice is `columnInput`.',
    why: 'A list built from `df.columns` offers the hidden `~` columns and nameless ones the platform\'s pickers ' +
      'leave out.',
    detail: 'Picking a column of a table is u2/dg\'s `columnInput(label, df, …)` (or `Platform.listedColumns(df, ' +
      'filter)` for a list of your own): they list what the platform\'s pickers list, hidden `~` columns left ' +
      'out. Never a `ChoiceInput` over `df.columns`.',
    wrong: 'new ChoiceInput({label: \'Column\', items: df.columns.names()})',
    right: 'columnInput(\'Column\', df, {filter: (c) => c.semType === DG.SEMTYPE.MOLECULE})',
    guards: ['lint:u2/no-raw-column-list', 'choke:column listing'],
    kinds: {
      'choice-input': {notFor: 'a column choice: that is `columnInput(label, df)` (u2/dg), never a list built ' +
        'from `df.columns`'},
      'column-combo': {whenToUse: 'a column choice in the platform\'s column combo: ' +
        '`columnInput(label, df, {grid: true})`, never a list built from `df.columns`'},
      'column-picker': {whenToUse: 'a column choice with a search over the columns: ' +
        '`columnInput(label, df, {rich: true})`, never a list built from `df.columns`'},
    },
  },
  {
    id: 'table-choice',
    rule: 'A table choice is `tableInput`.',
    why: 'A list of the shell\'s tables offers a table with no name, which no picker can hold.',
    detail: '`tableInput(label)` (u2/dg) offers the open tables the platform lists and holds one of them; ' +
      '`Platform.tableNames()` is that list for a picker of your own. A table with no name is never offered: a ' +
      'choice list holding one reports `item-null`, and a table input given a table it does not offer reports ' +
      '`table-not-offered`.',
    wrong: 'new ChoiceInput({label: \'Table\', items: grok.shell.tables.map((t) => t.name)})',
    right: 'tableInput(\'Table\')',
    guards: ['usage:item-null', 'usage:table-not-offered', 'choke:table names'],
    kinds: {'choice-input': {notFor: 'a table choice: that is `tableInput(label)` (u2/dg)'}},
  },
  {
    id: 'format',
    rule: 'One text for a value: `Format`.',
    why: 'A value reads the same everywhere u2 shows it: a date in the user\'s platform format, in one zone.',
    detail: '`Format.value(v)` is the text u2 shows for a value — a date in the user\'s platform date format, read ' +
      'in the browser\'s zone; `Format.date(d, {format})` takes a platform pattern; a Date on midnight UTC that ' +
      'stands for a day reads in UTC with `{utcDates: true}` (docs/recipes/forms-and-details.md).',
    wrong: 'span(String(row.created))',
    right: 'span(Format.value(row.created))',
    guards: ['choke:formatter'],
  },
  {
    id: 'comments',
    rule: 'Comment only what the code cannot say.',
    why: 'A narrative goes stale the day the code changes, and the code, the README and the history already say ' +
      'it.',
    detail: 'Doc blocks go on the API a caller reads; no file-header narratives and no doc blocks on private or ' +
      'protected members; a one-line `//` says why where the code cannot.',
    wrong: '/* The campaign editor: fields, stages and the save flow, reworked for the new templates. */',
    right: '// a rename fires onColumnNameChanged only, so the held name is remapped here',
    guards: ['lint:u2/comments', 'ratchet:fileHeaders', 'ratchet:privateDocs'],
  },
];

/** The rules, each one nothing checks yet saying so after its detail. */
export const U2_RULES: readonly U2Rule[] = RULES.map((r) => r.guards.length > 0 ? r :
  {...r, detail: `${r.detail} Docs only: no check catches a slip yet.`});
