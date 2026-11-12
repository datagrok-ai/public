'use strict';

const tester = require('./tester');
const rule = require('../rules/build-in-owner');

const options = [{controls: ['Control', 'Form', 'TextInput', 'Dialog', 'NameEditor'], factories: ['textInput']}];
const outside = (code, what, cls, where, opts = options) => ({code, options: opts,
  errors: [{messageId: 'outside', data: {what, cls, where}}]});

tester.run('build-in-owner', rule, {
  valid: [
    {code: `class NameEditor extends Control {
      readonly nameInput = this.runInScope(() => new TextInput({label: 'Name'}));
    }`, options},
    {code: `class NameEditor extends Control {
      constructor() {
        super(div());
        this.runInScope(() => { this.form = new Form(); this.form.add(textInput('Name')); });
      }
    }`, options},
    {code: `class NameEditor extends Control {
      readonly nameInput = this.build(() => new TextInput({label: 'Name'}));
    }`, options},
    {code: `class NameEditor extends Control {
      constructor() { super(div()); this.on(this.root, 'click', () => new Dialog('Rename').show()); }
      open() { return new Dialog('Rename'); }
    }`, options},
    {code: `class Model { readonly nameInput = new TextInput({label: 'Name'}); }`, options},
    {code: `class NameEditor extends Control {
      readonly box = div();
      readonly open = () => new Dialog('x');
    }`, options},
    {code: `class HitHome extends Control {
      refresh() { this.root.replaceChildren(this.build(() => new Form())); }
      onClick() { this.on(this.root, 'click', () => new Dialog('Rename').show()); }
    }`, options: [{...options[0], methods: true}]},
  ],
  invalid: [
    outside(`export class NameEditor extends Control {
      readonly nameInput = new TextInput({label: 'Name'});
    }`, '`new TextInput(…)`', 'NameEditor', 'field initializer'),
    outside(`class Wizard extends Control {
      constructor() { super(div()); this.form = new Form(); }
    }`, '`new Form(…)`', 'Wizard', 'constructor'),
    outside(`class Wizard extends NameEditor {
      readonly name = textInput('Name');
    }`, '`textInput(…)`', 'Wizard', 'field initializer'),
    outside(`class Wizard extends Control {
      readonly body = Control.build(() => [div()]);
    }`, '`Control.build(…)`', 'Wizard', 'field initializer', [{controls: ['Control'], factories: []}]),
    outside(`class Picker extends TextInput { readonly menu = new Dialog('Pick'); }`,
      '`new Dialog(…)`', 'Picker', 'field initializer', []),
    outside(`class HitHome extends Control {
      refresh() { this.root.replaceChildren(new Form()); }
      rebuild() { this.root.replaceChildren(this.build(() => new Form())); }
      static of() { return new Form(); }
    }`, '`new Form(…)`', 'HitHome', 'method refresh', [{...options[0], methods: true}]),
    outside(`class HitHome extends Control {
      get editor() { return textInput('Name'); }
    }`, '`textInput(…)`', 'HitHome', 'getter editor', [{...options[0], methods: true}]),
  ],
});
