import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, Subject} from 'rxjs';

import {BoolInput, ChoiceInput, ColorInput, Dialog, Form, LabelColumn, RecordsInput, Scope, Splitter, TextArea,
  TextInput, VirtualTree, component, computed, div, divV, each, h3, iconButton, notify, signal, statusLine,
  tableFromMap, untracked} from '@datagrok-libraries/u2';
import type {TreeNode} from '@datagrok-libraries/u2';
import {columnInput, inputForProperty, tableInput} from '@datagrok-libraries/u2/dg';
import {MpoDesirabilityLineEditor} from '@datagrok-libraries/statistics/src/mpo/editors/mpo-line-editor';

import {DEFAULT_TEMPLATE_FILE, Model, ModelColoring, ModelProperty, Subgroup, Template, TEMPLATES_FOLDER}
  from './constants';
import {getTemplates} from './admetica-utils';
import '../css/admetica.css';

type Item = Subgroup | Model;

const isGroup = (item: Item): item is Subgroup => 'models' in item;
const isMolecule = (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE;

export function admeticaFuncCallEditor(call: DG.FuncCall, table: DG.DataFrame): DG.FuncCallEditor {
  const scope = Scope.root('admetica-editor');
  const changes = new Subject<void>();
  const editor = scope.build(() => component((self) => {
    self.onError = (e) => notify.error(`Admetica templates: ${e instanceof Error ? e.message : String(e)}`);
    const template = signal<Template | null>(null);
    const saved = signal('');
    const edits = signal(0);
    const changed = () => edits.value++;
    const dirty = computed(() => {
      edits.value;
      const t = template.value;
      return t !== null && JSON.stringify(t) !== saved.value;
    });
    const invalidColors = signal(false);

    const tableIn = tableInput('Table', {value: table.name});
    const molecules = columnInput('Molecules', table, {grid: true, filter: isMolecule});
    self.effect(() => {
      const t = tableIn.table;
      untracked(() => {
        if (t !== molecules.table)
          molecules.changeTable(t, isMolecule);
        if (molecules.peek() === null)
          molecules.value = molecules.names()[0] ?? null;
      });
    });

    const templates = signal<string[]>([]);
    const profile = new ChoiceInput({label: 'Profile', items: templates, nullable: false});
    profile.addOptions(iconButton('save', () => saveTemplate(), {tooltip: 'Save as a new template',
      disabled: computed(() => dirty.value ? null : 'The template has no changes to save')}));
    profile.addOptions(iconButton('folder', () => manageTemplates(), {tooltip: 'Manage templates (rename/delete)'}));
    loadTemplates();

    self.effect(() => {
      const name = profile.value;
      if (name === null)
        return;
      untracked(() => {
        template.value = null;
        self.task(async () => {
          const loaded: Template = JSON.parse(await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${name}.json`));
          if (profile.peek() !== name)
            return;
          saved.value = JSON.stringify(loaded);
          template.value = loaded;
        });
      });
    });

    const tree = new VirtualTree<Item>({onCheck: (node, checked) => check(node.data!, checked)});
    tree.name = 'properties';
    tree.root.classList.add('admetica-editor-tree');
    self.effect(() => {
      edits.value;
      tree.setRoots(template.value?.subgroup.map(groupNode) ?? []);
    });
    self.effect(() => {
      const t = template.value;
      const selected = tree.selectedNode.peek()?.id;
      if (t === null)
        return;
      untracked(() => {
        tree.expanded.value = new Set(t.subgroup.filter((s) => s.expanded).map((s) => s.name));
        const group = t.subgroup.find((s) => s.name === selected || s.models.some((m) => m.name === selected)) ??
          t.subgroup[0];
        void tree.expandPath(group.name === selected ? [group.name] : [group.name, selected ?? group.models[0].name]);
      });
    });

    const detail = div([], 'admetica-editor-detail');
    each(detail, () => {
      const node = tree.selectedNode.value;
      return node?.data === undefined ? [] : [node.data];
    }, (item) => item.name, (item) => isGroup(item) ? groupDetail(item) : modelDetail(item));

    const addPiechart = new BoolInput({label: 'Add piechart', value: true});
    const addForm = new BoolInput({label: 'Add form', value: true});
    const labelColumn = new LabelColumn();

    const models = computed(() => {
      edits.value;
      return template.value?.subgroup.flatMap((s) => s.models).filter((m) => m.checked).map((m) => m.name) ?? [];
    });
    const problem = computed(() => {
      if (tableIn.table === null)
        return 'Pick a table';
      if (molecules.value === null)
        return 'The table has no molecule column';
      if (template.value === null)
        return 'The template is not loaded';
      if (invalidColors.value)
        return 'Fix the color ranges';
      return models.value.length === 0 ? 'Check the properties to calculate' : null;
    });

    self.effect(() => {
      const t = tableIn.table;
      const col = molecules.value;
      const checked = models.value;
      const piechart = addPiechart.value;
      const form = addForm.value;
      edits.value;
      untracked(() => {
        call.inputs['table'] = t;
        call.inputs['molecules'] = t === null || col === null ? null : t.col(col);
        call.inputs['template'] = JSON.stringify(template.value);
        call.inputs['models'] = checked;
        call.inputs['addPiechart'] = piechart;
        call.inputs['addForm'] = form;
        changes.next();
      });
    });

    const split = new Splitter([tree, detail], {direction: 'horizontal', sizes: [0.36, 0.64]});
    split.root.classList.add('admetica-editor-models');
    return {
      root: div([
        new Form({labelColumn}).addAll([tableIn, molecules, profile]),
        split,
        statusLine(problem),
        new Form({labelColumn}).addAll([addPiechart, addForm]),
      ], 'admetica-editor'),
      valid: computed(() => problem.value === null),
    };

    function groupNode(group: Subgroup): TreeNode<Item> {
      return {id: group.name, label: group.name, data: group, checked: group.checked,
        children: group.models.map((m) => ({id: m.name, label: m.name, tooltip: m.units, data: m,
          checked: m.checked}))};
    }

    function check(item: Item, checked: boolean): void {
      const group = isGroup(item) ? item : template.peek()!.subgroup.find((s) => s.models.includes(item))!;
      for (const m of isGroup(item) ? group.models : [item])
        m.checked = checked;
      group.checked = group.models.every((m) => m.checked);
      changed();
    }

    function groupDetail(group: Subgroup): HTMLElement {
      return divV([h3(group.name), div([group.description])], {gap: 'm'});
    }

    function modelDetail(model: Model): HTMLElement {
      return component((detailScope) => {
        const parts: HTMLElement[] = [h3(model.name)];
        for (const p of model.properties) {
          const part = propertyPart(p, model);
          if (part)
            parts.push(part);
        }
        if (model.line) {
          const line = new MpoDesirabilityLineEditor(model, 360, 120);
          detailScope.subscribe(line.onChanged, () => changed());
          parts.push(line.root);
        }
        if (model.coloring?.type === 'Conditional')
          parts.push(conditionalColors(model).root);
        else if (model.coloring?.type === DG.COLOR_CODING_TYPE.LINEAR)
          parts.push(colorScheme(model.coloring));
        return [divV(parts, {gap: 'm'})];
      }, {name: 'admetica-model'}).root;
    }

    function propertyPart(p: ModelProperty, model: Model): HTMLElement | null {
      const {name, inputType, enable, skip} = p.property;
      if (skip)
        return null;
      if (inputType === DG.InputType.Map)
        return tableFromMap(p.object[name]);
      if (!enable)
        return div([String(p.object[name] ?? '')]);
      const options = {value: p.object[name], onInput: (value: any) => {
        p.object[name] = value;
        if (name === 'direction' && model.coloring.type === DG.COLOR_CODING_TYPE.LINEAR)
          [model.coloring.min, model.coloring.max] = [model.coloring.max, model.coloring.min];
        changed();
      }};
      return (inputType === DG.InputType.TextArea ? new TextArea({...options, autoGrow: true, fill: true}) :
        inputForProperty(DG.Property.fromOptions(p.property), options)).root;
    }

    function colorScheme(coloring: ModelColoring): HTMLElement {
      const canvas = ui.schemeInput(JSON.parse(coloring.colors!)).querySelector('canvas')!;
      canvas.classList.add('admetica-color-scheme');
      return canvas;
    }

    function conditionalColors(model: Model): RecordsInput<{range: string, color: string}> {
      const {type, ...colors} = model.coloring as unknown as Record<string, string>;
      const ranges = model.properties.find((p) => p.property.name === 'ranges')?.object.ranges;
      let keys = Object.keys(colors);
      invalidColors.value = false;
      const problemOf = (rules: {range: string}[]) => rules.some((r) => r.range.trim() === '') ?
        'A range cannot be empty' :
        new Set(rules.map((r) => r.range.trim())).size < rules.length ? 'Each range must be different' : null;
      const records = new RecordsInput({name: 'Colors', allowAdd: false, allowRemove: false,
        value: keys.map((range) => ({range, color: colors[range]})),
        columns: [
          {name: 'range', header: 'Range', input: (range: string) => new TextInput({value: range, inline: true,
            commitOn: 'change', tooltipText: ranges?.[range] ?? ''})},
          {name: 'color', header: 'Color', input: (color: string) => new ColorInput({value: color, inline: true,
            swatchOnly: true})},
        ],
        onInput: (rules) => {
          invalidColors.value = problemOf(rules) !== null;
          if (invalidColors.value)
            return;
          const coloring: Record<string, string> = {type};
          rules.forEach((rule, i) => {
            const range = rule.range.trim();
            coloring[range] = rule.color;
            if (ranges && range !== keys[i]) {
              ranges[range] = ranges[keys[i]];
              delete ranges[keys[i]];
            }
          });
          model.coloring = coloring as unknown as ModelColoring;
          keys = rules.map((rule) => rule.range.trim());
          changed();
        }});
      records.addValidator(problemOf);
      return records;
    }

    function loadTemplates(select?: string): void {
      self.task(async () => {
        const names = await getTemplates();
        templates.value = names;
        const builtIn = DEFAULT_TEMPLATE_FILE.replace(/\.json$/, '');
        profile.value = [select, profile.peek(), builtIn].find((n) => n != null && names.includes(n)) ??
          names[0] ?? null;
      });
    }

    function saveTemplate(): void {
      self.build(() => {
        const name = new TextInput({label: 'Name', required: true});
        name.addValidator((value) => {
          const n = value.trim();
          if (/[./\\]/.test(n))
            return 'A template name cannot contain ".", "/" or "\\"';
          if (`${n}.json` === DEFAULT_TEMPLATE_FILE)
            return 'The built-in default template cannot be overwritten';
          return templates.value.includes(n) ? 'Template with this name already exists' : null;
        });
        new Dialog('Save template')
          .add(new Form().add(name))
          .onOK(async () => {
            const t = template.peek()!;
            for (const group of t.subgroup)
              group.expanded = tree.expanded.peek().has(group.name);
            await grok.dapi.files.writeAsText(`${TEMPLATES_FOLDER}/${name.peek().trim()}.json`, JSON.stringify(t));
            notify.info(`Template "${name.peek().trim()}" has been saved`);
            loadTemplates(name.peek().trim());
          }, {text: 'Save', primary: true})
          .show();
      });
    }

    function manageTemplates(): void {
      ui.dialog({title: 'Manage files'})
        .add(ui.fileBrowser({path: TEMPLATES_FOLDER}).root)
        .onOK(() => loadTemplates())
        .show();
    }
  }, {name: 'admetica-editor'}));
  DG.Widget.registerCleanup(editor.root, () => scope.dispose());
  return new class extends DG.FuncCallEditor {
    get isValid(): boolean {
      return editor.valid.peek();
    }

    get onInputChanged(): Observable<void> {
      return changes;
    }
  }(editor.root);
}
