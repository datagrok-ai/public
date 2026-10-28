'use strict';

const tester = require('./tester');
const rule = require('../rules/no-unowned-root');

const U2 = `import {AppHeader, Form, TextInput, Dialog, component, Scope} from '@datagrok-libraries/u2';
import {ViewHost, funcForm} from '@datagrok-libraries/u2/dg';
`;
const options = [{controls: ['Control', 'AppHeader', 'Form', 'TextInput', 'Dialog', 'ViewHost'],
  roots: ['Dialog', 'ViewHost'], factories: ['funcForm']}];
const ok = (code) => ({code: U2 + code, options});
const unowned = (code, what, messageId = 'unownedControl') => ({code: U2 + code, options,
  errors: [{messageId, data: {what}}]});

tester.run('no-unowned-root', rule, {
  valid: [
    `const built = Scope.runWith(this.scope, () => new HitHome(options));`,
    `const scope = new Scope(Ownership.scopeOf(owner, 'subscribe()'), {label: 'subscription'});`,
    `const home = Scope.root('campaign landing').build(() => new HitHome(options));`,
    `const s = new DG.U2.Scope(parent);`,
    `runWith(undefined, fn);`,
    ok(`export function app() {
      const view = DG.View.create();
      view.root.append(ViewHost.of(view).build(() => new AppHeader({icon: 'flask'})).root);
      return view;
    }`),
    ok(`export const home = () => ViewHost.create({name: 'Home', content: () => new AppHeader({icon: 'flask'})}).view;`),
    ok(`export const editor = () => component(() => {
      const name = new TextInput({label: 'Name'});
      return {root: new Form({inputs: [name]}).root, name};
    });`),
    ok(`class MolTrackEms {
      static home(view) { return ViewHost.of(view).build(() => MolTrackEms._header()); }
      private static _header() { return new AppHeader({icon: 'flask'}).root; }
    }`),
    ok(`export function nameInput() { return new TextInput({label: 'Name'}); }
    export const fields = (names) => names.map((n) => new TextInput({label: n}));`),
    ok(`export function inputsPage(): HTMLElement {
      const name = new TextInput({label: 'Name'});
      return divV([name.root]);
    }`),
    ok(`export function open() { new Dialog('Rename').show(); }`),
    ok(`export function edit(dialog) { dialog.add(() => new TextInput({label: 'Name'})); }`),
    ok(`export function load(host) {
      host.task(async (t) => { const rows = await fetchRows(); t.build(() => new Form({rows})); });
    }`),
    ok(`export class Editor extends Control { refresh() { this.root.append(new Form().root); } }`),
    ok(`import {Form as U2Form} from '@datagrok-libraries/u2';
    export const x = () => { const f = new LocalForm(); f.show(); };`),
    {code: U2 + `export function app() { DG.View.fromRoot(new AppHeader({icon: 'flask'}).root); }`,
      options: [{controls: [], factories: []}]},
  ],
  invalid: [
    {code: `this.home = Scope.runWith(undefined, () => new HitHome(this.homeOptions(preset)));`,
      errors: [{messageId: 'runWith'}]},
    {code: `const editor = Scope.runWith(null, () => new TemplateEditor(app));`, errors: [{messageId: 'runWith'}]},
    {code: `DG.U2.Scope.runWith(void 0, () => build());`, errors: [{messageId: 'runWith'}]},
    {code: `const scope = new Scope();`, errors: [{messageId: 'newScope'}]},
    {code: `const scope = new Scope(undefined, {label: 'dialog'});`, errors: [{messageId: 'newScope'}]},
    {code: `const scope = new DG.U2.Scope();`, errors: [{messageId: 'newScope'}]},
    unowned(`export class PackageFunctions {
      static async benchlingLinkApp() {
        const appHeader = new AppHeader({icon: 'benchling.png', items: ['Integrate with your Benchling account.']});
        const view = DG.View.fromRoot(appHeader.root);
        view.name = 'Benchling';
        return view;
      }
    }`, 'new AppHeader(…)'),
    unowned(`function getAppHeader() {
      const appHeader = new AppHeader({icon: 'revvity.png', title: 'RevvityLink'});
      return appHeader.root;
    }
    export function revvityApp() { return DG.View.fromRoot(ui.divV([getAppHeader()])); }`, 'new AppHeader(…)'),
    unowned(`export class Boltz1AppView {
      constructor() { this.divV = ui.divV([]); this.createHeader(); }
      createHeader() { this.divV.append(new AppHeader({icon: 'boltz.png'}).root); }
    }`, 'new AppHeader(…)'),
    unowned(`export async function openApp(standard) {
      const studies = await loadStudies(standard);
      const header = new AppHeader({icon: 'clin.png'});
      return ui.divV([header.root, studies]);
    }`, 'new AppHeader(…)', 'afterAwait'),
    unowned(`export function load(host) {
      host.task(async () => { const rows = await fetchRows(); host.root.append(new Form({rows}).root); });
    }`, 'new Form(…)', 'afterAwait'),
    unowned(`export const pane = (func) => { const form = funcForm(func); grok.shell.o = form.root; };`, 'funcForm(…)'),
    unowned(`export function wire(button) {
      button.addEventListener('click', () => { document.body.append(new TextInput({label: 'Name'}).root); });
    }`, 'new TextInput(…)'),
    unowned(`const header = new DG.U2.Control();
    document.body.append(header.root);`, 'new Control(…)'),
    {code: `import * as u2 from '@datagrok-libraries/u2';
    export function app() { return DG.View.fromRoot(new u2.AppHeader({icon: 'flask'}).root); }`,
    errors: [{messageId: 'unownedControl', data: {what: 'new AppHeader(…)'}}]},
  ],
});
