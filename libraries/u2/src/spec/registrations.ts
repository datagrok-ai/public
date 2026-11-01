/* Every component a dg-ui/1 spec can build, with the prop metadata that validates specs and feeds
   manifest.json. Left out: components whose options are functions (list, tree, combobox, dialog,
   toolbar, async view). */
import {Signal, signal} from '../core/signals.js';
import {Scope} from '../core/scope.js';
import {Control} from '../core/component.js';
import {button, divH, divV, panel, span, type SpaceSize, type Text} from '../core/elements.js';
import {icon, IconVariant} from '../components/display/icon.js';
import {Input, InputOptions, LiveOption} from '../core/input-base.js';
import type {Reason} from '../core/reason.js';
import {TextInput, TextArea} from '../components/inputs/text-input.js';
import {BoolInput} from '../components/inputs/bool-input.js';
import {NumberInput} from '../components/inputs/number-input.js';
import {ChoiceInput, MultiChoiceInput} from '../components/inputs/choice-input.js';
import {SliderInput} from '../components/inputs/slider-input.js';
import {RangeSlider} from '../components/inputs/range-slider.js';
import {RadioInput} from '../components/inputs/radio-input.js';
import {ColorInput} from '../components/inputs/color-input.js';
import {ListInput} from '../components/inputs/list-input.js';
import {MapInput} from '../components/inputs/map-input.js';
import {RecordsInput} from '../components/inputs/records-input.js';
import {ChecklistInput} from '../components/inputs/checklist-input.js';
import {ChipsInput} from '../components/inputs/chips-input.js';
import {AccessGrid} from '../components/forms/access-grid.js';
import type {AccessRow} from '../components/forms/access-grid.js';
import {QNum} from '../core/qnum.js';
import {QNumInput} from '../components/inputs/qnum-input.js';
import {FontInput} from '../components/inputs/font-input.js';
import {ImageInput} from '../components/inputs/image-input.js';
import {IconInput} from '../components/inputs/icon-input.js';
import {Form} from '../components/forms/form.js';
import type {FormLayout} from '../components/forms/form.js';
import {Section} from '../components/containers/section.js';
import {Page} from '../components/containers/page.js';
import {ChoiceChip} from '../components/inputs/choice-chip.js';
import {Splitter} from '../components/containers/splitter.js';
import {Accordion} from '../components/containers/accordion.js';
import {Card} from '../components/containers/card.js';
import {StatCard} from '../components/display/stat-card.js';
import {AppHeader} from '../components/display/app-header.js';
import {DataTable} from '../components/collections/data-table.js';
import type {DataTableColumn} from '../components/collections/data-table.js';
import {CatalogTable} from '../components/collections/catalog-table.js';
import type {CatalogColumn} from '../components/collections/catalog-table.js';
import {TabStrip, TabStripOptions} from '../components/containers/tabs.js';
import {Wizard} from '../components/containers/wizard.js';
import {PropertyGrid, PropDescriptor} from '../components/forms/property-grid.js';
import {Breadcrumbs} from '../components/navigation/breadcrumbs.js';
import {FilterBuilder} from '../components/filter/filter-builder.js';
import type {FilterMode} from '../components/filter/filter-builder.js';
import {FilterQueryInput} from '../components/filter/filter-query-input.js';
import {Filters} from '../core/filter/index.js';
import type {FilterGroup, FilterSchema, FilterScalar, FilterTarget, FilterTemplate} from '../core/filter/index.js';
import type {IProperty} from '../core/property-like.js';
import {StateSource} from '../sources/state.js';
import {registerDataSources} from '../sources/registrations.js';
import {ComponentMeta, SpecPropMeta, Registry, registry as globalRegistry} from './registry.js';
import {childProp} from './spec.js';
import type {SpecNode} from './spec.js';

type Props = Record<string, unknown>;
type Child = Control | HTMLElement;

/** A bound `value` arrives as the context signal itself, so the input adopts it instead of copying:
 * no bridge, no echo, and the spec's two-way contract holds through every edit. */
function inputOptions<T>(props: Props): InputOptions<T> {
  return {
    label: props.label as LiveOption<string> | undefined,
    name: props.name as string | undefined,
    inline: props.inline as boolean | undefined,
    postfix: props.postfix as LiveOption<string> | undefined,
    tooltipText: props.tooltipText as LiveOption<string> | undefined,
    disabled: props.disabled as Reason | undefined,
    required: props.required as LiveOption<boolean> | undefined,
    ...Input.specValue<T>(props.value, props.readOnly as LiveOption<boolean> | undefined),
  };
}

/** A bound `items` list follows its signal instead of being copied once — how a data source feeds
 * a picker (`bind: {"items": "$.people.names"}`). One-way: the list is data, not an edit. */
function boundItems(input: Control & {setItems(items: string[]): void}, items: unknown): void {
  if (items instanceof Signal)
    input.effect(() => input.setItems(items.value as string[] ?? []));
}

/** `['name', 'cas:CAS number']` — the property, and after a colon the header it gets. */
function dataTableColumns(spec: unknown): DataTableColumn<Record<string, unknown>>[] {
  return itemList(spec).map((entry) => DataTable.column(entry));
}

function itemList(items: unknown): string[] {
  return (items instanceof Signal ? items.peek() : items) as string[] ?? [];
}

function inputProps(value: string, ...extra: SpecPropMeta[]): SpecPropMeta[] {
  // chrome is live-bindable one-way (UB-10); `inline` changes the layout and is structural
  return [
    {name: 'label', type: 'string', bindable: true},
    {name: 'name', type: 'string', description: 'Stable key for forms and dumps; defaults to the label.'},
    {name: 'value', type: value, bindable: true, twoWay: true},
    {name: 'inline', type: 'bool', structural: true, description: 'Compact one-row variant without a label.'},
    {name: 'postfix', type: 'string', bindable: true,
      description: 'Units or a short suffix shown after the editor.'},
    {name: 'tooltipText', type: 'string', bindable: true},
    {name: 'disabled', type: 'string', bindable: true,
      description: 'Why the input is off — greyed out, edits blocked, the reason as its tooltip; empty while on.'},
    {name: 'readOnly', type: 'bool', bindable: true,
      description: 'Shows the value, selectable, but blocks edits.'},
    {name: 'required', type: 'bool', bindable: true,
      description: 'An empty value is invalid once the input was edited, left or its form validated.'},
    ...extra,
  ];
}

/** A spec's `schema` is either a literal `{properties, values?}` or a bound `FilterSchema` (its
 * `values` is a function); the literal goes through `Filters.schema`. */
function filterSchema(x: unknown): FilterSchema {
  const raw = (x instanceof Signal ? x.peek() : x) as
    {properties?: IProperty[], values?: unknown} | undefined;
  if (!raw || !Array.isArray(raw.properties))
    return Filters.schema([]);
  return typeof raw.values === 'function' ? raw as FilterSchema :
    Filters.schema(raw.properties, raw.values as Record<string, FilterScalar[]> | undefined);
}

/** A spec's `template` is a literal whose `root` carries no node ids; the builder gets one with
 * them, so the value it seeds can be followed by `Filters.diff`. */
function filterTemplate(x: unknown): FilterTemplate | undefined {
  const t = x as FilterTemplate | undefined;
  if (!t?.root)
    return undefined;
  const root = FilterBuilder.withIds(t.root);
  return root === t.root ? t : {...t, root};
}

function element(child: Child): HTMLElement {
  return Control.is(child) ? child.root : child;
}

/** The title a container reads off its child's spec node — the one prop a parent owns. */
function childTitle(node: SpecNode | undefined, fallback: string): string {
  const title = childProp(node, 'title');
  return typeof title === 'string' ? title : fallback;
}

function childIcon(node: SpecNode | undefined): string | undefined {
  const icon = childProp(node, 'icon');
  return typeof icon === 'string' && icon !== '' ? icon : undefined;
}

function tabStrip(props: Props): TabStripOptions {
  return {
    orientation: props.orientation === 'vertical' ? 'vertical' : 'horizontal',
    variant: props.variant === 'document' ? 'document' : 'platform',
  };
}

function splitter(props: Props, children: Child[]): Splitter {
  return new Splitter(children.map(element), {
    direction: props.direction === 'vertical' ? 'vertical' : 'horizontal',
    sizes: props.sizes as number[] | undefined,
    minSize: props.minSize as number | undefined,
  });
}

function card(props: Props, children: Child[]): Card {
  return new Card({
    title: props.title as string | undefined,
    subtitle: props.subtitle as string | undefined,
    icon: props.icon as string | undefined,
    media: props.media as string | undefined,
    clickable: props.clickable as boolean | undefined,
    selectable: props.selectable as boolean | undefined,
    selected: props.selected as boolean | Signal<boolean> | undefined,
    children: children.map(element),
  });
}

/** Every spec child is a step; a wizard needs one, so an empty node gets a blank first step. */
function wizard(children: Child[], nodes: SpecNode[]): Wizard {
  const steps = children.map((child, i) => ({
    id: String(childProp(nodes[i], 'id') ?? `step-${i}`), title: childTitle(nodes[i], `Step ${i + 1}`),
    content: element(child), done: childProp(nodes[i], 'done') === true,
  }));
  return new Wizard({steps: steps.length > 0 ? steps : [{id: 'step-0', title: 'Step 1', content: divV([])}]});
}

function page(props: Props, children: Child[]): Page {
  return new Page({content: children, maxWidth: props.maxWidth as number | undefined});
}

function section(props: Props, children: Child[]): Section {
  const s = new Section({
    title: (props.title as Text | undefined) ?? '',
    expanded: props.expanded as boolean | Signal<boolean> | undefined,
    collapsible: props.collapsible as boolean | undefined,
  });
  return s.add(...children.map(element));
}

/** `u2-button`: text, primary, icon and disabled are written in place — the element stays. */
class SpecButton extends Control {
  private declare _text: string;
  private declare _icon: string;
  private declare _disabled: Signal<string | null>;

  constructor(props: Props) {
    const disabled = signal<string | null>((props.disabled as string | undefined) ?? null);
    // the button follows the reason under the control's own scope
    const scope = new Scope();
    super(scope.build(() => button('', () => {}, {primary: props.primary === true, disabled})), scope);
    this._disabled = disabled;
    this._text = (props.text as string | undefined) ?? '';
    this._icon = (props.icon as string | undefined) ?? '';
    this._paint();
  }

  /** Why the button is off — its tooltip, and what a screen reader says on focus; empty or null on. */
  get disabled(): string | null {
    return this._disabled.peek();
  }

  set disabled(x: string | null) {
    this._disabled.value = x || null;
  }

  get text(): string {
    return this._text;
  }

  set text(x: string) {
    this._text = x ?? '';
    this._paint();
  }

  get icon(): string {
    return this._icon;
  }

  set icon(x: string) {
    this._icon = x ?? '';
    this._paint();
  }

  get primary(): boolean {
    return this.root.classList.contains('u2-btn-primary');
  }

  set primary(x: boolean) {
    this.root.classList.toggle('u2-btn-primary', x === true);
  }

  /** What `button` and `buttonWithIcon` put inside: the text alone, or the icon and a label. */
  private _paint(): void {
    this.root.replaceChildren(...this._icon === '' ? [this._text] :
      [icon(this._icon), span(this._text, 'u2-btn-label')]);
  }
}

/** A layout container whose extra `cls` and `gap` follow their props in place. */
class SpecContainer extends Control {
  private declare _cls: string[];
  private declare _gap: SpaceSize | '';

  get cls(): string {
    return (this._cls ?? []).join(' ');
  }

  set cls(x: string) {
    this.root.classList.remove(...this._cls ?? []);
    this._cls = String(x ?? '').split(/\s+/).filter((c) => c);
    this.root.classList.add(...this._cls);
  }

  get gap(): SpaceSize | '' {
    return this._gap ?? '';
  }

  set gap(x: SpaceSize | '') {
    if (this._gap)
      this.root.classList.remove(`u2-gap-${this._gap}`);
    this._gap = x ?? '';
    if (this._gap)
      this.root.classList.add(`u2-gap-${this._gap}`);
  }
}

/** `u2-form-buttons`: a row of buttons, 8px apart; inside a u2-form it is the form's button row. */
class SpecFormButtons extends SpecContainer {
  constructor(readonly align: 'end' | 'editors' | 'start' | undefined) {
    super(divH([], {gap: 'm'}));
    // no kind of its own: to automation it is a horizontal stack, and features address its buttons, never the row
    this.root.dataset.u2 = 'div-h';
  }
}

/** `u2-icon`: variant and size follow their props in place; the glyph is read once (`name` is the
 * node's identity accessor). */
class SpecIcon extends Control {
  private declare _glyph: string;
  private declare _variant: IconVariant | undefined;
  private declare _size: IconSize | undefined;

  constructor(props: Props) {
    const glyph = (props.name as string | undefined) ?? 'question';
    const options = {variant: props.variant as IconVariant | undefined, size: props.size as IconSize | undefined};
    super(icon(glyph, options));
    this._glyph = glyph;
    this._variant = options.variant;
    this._size = options.size;
  }

  get variant(): IconVariant | undefined {
    return this._variant;
  }

  set variant(x: IconVariant | undefined) {
    this._variant = x;
    this._paint();
  }

  get size(): IconSize | undefined {
    return this._size;
  }

  set size(x: IconSize | undefined) {
    this._size = x;
    this._paint();
  }

  private _paint(): void {
    this.root.className = icon(this._glyph, {variant: this._variant, size: this._size}).className;
  }
}

type IconSize = 'small' | 'medium' | 'large';

function container(tag: string, create: () => HTMLElement, description: string): ComponentMeta {
  return {
    tag,
    category: 'Containers',
    create: (props) => {
      const built = new SpecContainer(create());
      built.cls = (props.cls as string | undefined) ?? '';
      built.gap = (props.gap as SpaceSize | undefined) ?? '';
      return built;
    },
    description,
    props: [{name: 'cls', type: 'string', description: 'Extra layout class.'},
      {name: 'gap', type: 'string', choices: ['xs', 's', 'm', 'l', 'xl', 'xxl'],
        description: 'Space between the children on the token scale: 2, 4, 8, 12, 16 or 24px.'}],
    acceptsChildren: true,
    example: {tag, children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
  };
}

const METAS: ComponentMeta[] = [
  {
    tag: 'u2-text-input',
    category: 'Inputs',
    create: (props) => new TextInput({
      ...inputOptions<string>(props),
      placeholder: props.placeholder as LiveOption<string> | undefined,
      search: props.search as boolean | undefined,
      password: props.password as boolean | undefined,
      autoResize: props.autoResize as boolean | undefined,
    }),
    description: 'Single-line text editor with validation.',
    usage: 'With `search: true` this is the filter box; in a platform view the main filter belongs ' +
      'in the view ribbon (the `ribbon` a `ViewHost.create` builder answers), not inside the content area. ' +
      'Every input: act on the user\'s edits with `onInput` (`onChanged` also hears writes from code); load ' +
      'another record with `reset(value)`, keep in step with a model with `follow(read)`; a status line or ' +
      'tab marker reads `error` (shown, touch-gated), Save reads `validity`.',
    props: inputProps('string',
      {name: 'placeholder', type: 'string', bindable: true},
      {name: 'search', type: 'bool', structural: true, description: 'Magnifier affordance plus a clear button.'},
      {name: 'password', type: 'bool', structural: true, description: 'Masked field with an eye toggle.'},
      {name: 'autoResize', type: 'bool', structural: true,
        description: 'Width follows the text between 100 and 300px.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-text-input', props: {label: 'Name', value: 'Aspirin'}},
  },
  {
    tag: 'u2-text-area',
    category: 'Inputs',
    create: (props) => new TextArea({
      ...inputOptions<string>(props),
      placeholder: props.placeholder as LiveOption<string> | undefined,
      autoGrow: props.autoGrow as boolean | undefined,
      rows: props.rows as number | undefined,
    }),
    description: 'Multi-line text editor.',
    props: inputProps('string',
      {name: 'placeholder', type: 'string', bindable: true},
      {name: 'autoGrow', type: 'bool', structural: true, description: 'Grows with the text instead of scrolling.'},
      {name: 'rows', type: 'int', description: 'Visible text lines.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-text-area', props: {label: 'Description', placeholder: 'Up to 40 characters'}},
  },
  {
    tag: 'u2-bool-input',
    category: 'Inputs',
    create: (props) => new BoolInput({
      ...inputOptions<boolean>(props),
      switch: props.switch as boolean | undefined,
      captionAfter: props.captionAfter as boolean | undefined,
    }),
    description: 'Checkbox or switch.',
    props: inputProps('bool',
      {name: 'switch', type: 'bool', structural: true, description: 'Compact toggle instead of a checkbox.'},
      {name: 'captionAfter', type: 'bool', structural: true, description: 'Checklist row: the box, then its caption.'}),
    events: ['change'],
    example: {tag: 'u2-bool-input', props: {label: 'Active', value: true, switch: true}},
  },
  {
    tag: 'u2-number-input',
    category: 'Inputs',
    create: (props) => new NumberInput({
      ...inputOptions<number | null>(props),
      mode: props.mode as 'int' | 'float' | undefined,
      min: props.min as number | undefined,
      max: props.max as number | undefined,
      step: props.step as number | undefined,
      slider: props.slider as boolean | undefined,
      clicker: props.clicker as boolean | undefined,
      spinner: props.spinner as boolean | undefined,
      placeholder: props.placeholder as LiveOption<string> | undefined,
    }),
    description: 'Numeric editor; text that does not parse or falls outside min/max stays on ' +
      'screen and only marks the input invalid.',
    props: inputProps('double',
      {name: 'mode', type: 'string', choices: ['int', 'float'], structural: true,
        description: '"int" or "float" (default).'},
      {name: 'min', type: 'double'},
      {name: 'max', type: 'double'},
      {name: 'step', type: 'double', description: 'Spinner, clicker and arrow-key increment; 1 by default.'},
      {name: 'slider', type: 'bool', structural: true,
        description: 'Range control under the box, revealed on hover; needs both bounds.'},
      {name: 'clicker', type: 'bool', structural: true,
        description: '− / + buttons on the options rail, revealed on hover.'},
      {name: 'spinner', type: 'bool', structural: true, description: 'Hover spinner inside the field.'},
      {name: 'placeholder', type: 'string', bindable: true}),
    events: ['input', 'change'],
    example: {tag: 'u2-number-input', props: {label: 'Amount', value: 10, min: 0, max: 100}},
  },
  {
    tag: 'u2-choice-input',
    category: 'Inputs',
    create: (props) => {
      const input = new ChoiceInput({
        ...inputOptions<string | null>(props),
        items: itemList(props.items),
        nullable: props.nullable as boolean | undefined,
      });
      boundItems(input, props.items);
      return input;
    },
    description: 'Single choice over a native select.',
    usage: 'In a form a choice is as wide as its longest item; `width` sets it outright (a short code ' +
      'list), `fill: true` stretches it across the row. Act on the user\'s pick with `onInput` — ' +
      '`onChanged` also hears the program\'s writes. `items` may be a signal, and the value may be set ' +
      'before its item arrives.',
    props: inputProps('string',
      {name: 'items', type: 'string_list', bindable: true},
      {name: 'nullable', type: 'bool', structural: true, description: 'Offers an empty option; true by default.'}),
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    events: ['change'],
    example: {tag: 'u2-choice-input', props: {label: 'Series', items: ['a', 'b', 'c'], value: 'b'}},
  },
  {
    tag: 'u2-choice-chip',
    category: 'Inputs',
    create: (props) => new ChoiceChip({label: props.label as string | undefined, name: props.name as string | undefined,
      items: itemList(props.items), value: (props.value as string | undefined) ?? null,
      disabled: props.disabled as Reason | undefined, invalid: props.invalid as Reason | undefined}),
    description: 'A small choice inline, reading "label value ▾": a press opens a menu of the choices, the current ' +
      'one checked.',
    usage: 'For a choice among others in a row or a heading (an R group\'s cap, a sort): two clicks to change it. In ' +
      'code an item takes a `description` (its grey detail in the menu and the chip\'s tooltip), `invalid` turns the ' +
      'chip red with its reason and `actions` add entries under the choices. A field of a form is `u2-choice-input`.',
    props: [
      {name: 'label', type: 'string', structural: true},
      {name: 'name', type: 'string', description: 'The automation id.'},
      {name: 'value', type: 'string'},
      {name: 'items', type: 'string_list', structural: true},
      {name: 'disabled', type: 'string', structural: true,
        description: 'Why the chip is off — greyed out, its menu closed; empty while on.'},
      {name: 'invalid', type: 'string', structural: true,
        description: 'Why the current choice cannot stand — the chip turns red; empty while it is fine.'},
    ],
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    example: {tag: 'u2-choice-chip', props: {label: 'Sort', items: ['Name', 'Created'], value: 'Name'}},
  },
  {
    tag: 'u2-multi-choice-input',
    category: 'Inputs',
    create: (props) => {
      const input = new MultiChoiceInput({
        ...inputOptions<string[]>(props),
        items: itemList(props.items),
        emptyText: props.emptyText as string | undefined,
        showSummaryCheckbox: props.showSummaryCheckbox as boolean | undefined,
      });
      boundItems(input, props.items);
      return input;
    },
    description: 'Checkbox list; the value holds the checked items in item order.',
    props: inputProps('string_list', {name: 'items', type: 'string_list', bindable: true},
      {name: 'emptyText', type: 'string', description: 'Stands in for an empty list; \'No items\' by default.'},
      {name: 'showSummaryCheckbox', type: 'bool', structural: true,
        description: 'Tri-state "<N> of <M>" summary row that toggles all/none; the list collapses on click.'}),
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    events: ['change'],
    example: {tag: 'u2-multi-choice-input', props: {label: 'Tags', items: ['acid', 'base'], value: ['acid']}},
  },
  {
    tag: 'u2-button',
    category: 'Actions',
    create: (props) => new SpecButton(props),
    description: 'Push button; `on: {"click": "cmd:<name>"}` runs a context command. Inside a u2-form it goes ' +
      'to the form\'s button row.',
    usage: 'Off with the reason: `disabled` takes why (bind it to a context value that is null while the button ' +
      'can run) and shows it as the tooltip — say what is missing rather than leaving the user to guess. In code, ' +
      '`button(text, onClick, {disabled: reason})`; a dialog\'s OK takes the same through `okDisabled(reason)`.',
    props: [
      {name: 'text', type: 'string'},
      {name: 'primary', type: 'bool', description: 'Accent styling for the main action.'},
      {name: 'icon', type: 'string', inputType: 'Icon',
        description: 'Font Awesome icon name shown before the text.'},
      {name: 'disabled', type: 'string', description: 'Why the button is off: its tooltip, said on keyboard ' +
        'focus too; a press does nothing. Empty or null: on.'},
    ],
    defaults: {text: 'Button'},
    events: ['click'],
    example: {tag: 'u2-button', props: {text: 'Save', primary: true}, on: {click: 'cmd:save'}},
  },
  container('u2-panel', () => {
    const el = panel();
    el.dataset.u2 = 'panel';
    return el;
  }, 'Padded vertical container.'),
  container('u2-div-v', () => {
    const el = divV();
    el.dataset.u2 = 'div-v';
    return el;
  }, 'Vertical stack.'),
  container('u2-div-h', () => {
    const el = divH();
    el.dataset.u2 = 'div-h';
    return el;
  }, 'Horizontal row.'),
  {
    tag: 'u2-form-buttons',
    category: 'Containers',
    create: (props) => new SpecFormButtons(props.align as 'end' | 'editors' | 'start' | undefined),
    description: 'A row of buttons; as a child of a u2-form it is the form\'s button row, under the rows.',
    usage: 'A form\'s Save and Cancel: `{"tag": "u2-form-buttons", "children": [{"tag": "u2-button", "props": ' +
      '{"text": "Save", "primary": true}, "bind": {"disabled": "$.saveProblem"}}]}` — the button off with the ' +
      'reason it is off.',
    props: [{name: 'align', type: 'string', choices: ['end', 'editors', 'start'], structural: true,
      description: 'In a form: the buttons at the end (default), starting where the editors do, or at the start.'}],
    acceptsChildren: true,
    example: {tag: 'u2-form-buttons', children: [{tag: 'u2-button', props: {text: 'Cancel'}},
      {tag: 'u2-button', props: {text: 'Save', primary: true}}]},
  },
  {
    tag: 'u2-splitter',
    category: 'Containers',
    create: (props) => splitter(props, []),
    createWithChildren: (props, children) => splitter(props, children),
    description: 'Resizable panels separated by draggable sashes; every spec child becomes a panel.',
    usage: 'Master-detail layouts (list left, details right): see the entity-browser recipe. In a ' +
      'platform view, consider handing details to the context panel (`grok.shell.o`) instead of a pane.',
    props: [
      {name: 'direction', type: 'string', choices: ['horizontal', 'vertical'], structural: true,
        description: '"horizontal" (default) or "vertical".'},
      {name: 'sizes', type: 'object', description: 'Fraction per panel, normalized; equal shares by default.'},
      {name: 'minSize', type: 'double', structural: true, description: 'Smallest panel size in pixels; 60 by default.'},
    ],
    acceptsChildren: true,
    defaults: {direction: 'horizontal'},
    // two, not one: a single-panel splitter has no sash to drag and is a panel with extra steps
    defaultChildren: [{tag: 'u2-panel'}, {tag: 'u2-panel'}],
    designerActions: [{text: 'Add panel',
      produce: () => ({op: 'add-child', node: {tag: 'u2-panel'}})}],
    example: {tag: 'u2-splitter', props: {direction: 'horizontal', sizes: [0.3, 0.7]}, children: [
      {tag: 'u2-panel', children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
      {tag: 'u2-panel', children: [{tag: 'u2-text-area', props: {label: 'Notes'}}]},
    ]},
  },
  {
    tag: 'u2-wizard',
    category: 'Containers',
    create: () => wizard([], []),
    createWithChildren: (_props, children, nodes) => wizard(children, nodes),
    // a designer rebuild keeps the user's place: the open step and the visited ones, by id
    carry: (_node, prev, next) => {
      if (prev instanceof Wizard && next instanceof Wizard)
        next.resumeFrom(prev);
    },
    setChildProp: (parent, index, name, value) => {
      if (name === 'title')
        (parent as Wizard).setStepTitle(index, typeof value === 'string' ? value : `Step ${index + 1}`);
    },
    description: 'Steps on a rail with BACK / NEXT / FINISH, one per spec child; the child node carries the ' +
      'step title. Panels are hidden, never rebuilt, so a step\'s inputs keep their state.',
    usage: 'For a flow whose steps are known up front. A spec expresses the steps, reacts to `finish`, ' +
      '`cancel` and `step` through `on`, and reads `currentStep`; gates (`canProceed`), a step\'s `commit`, ' +
      'its `actions` and dialog mode are code options.',
    props: [
      {name: 'currentStep', type: 'string', bindable: true,
        description: 'Id of the open step (read-only): bind a readout or a source param to it.'},
    ],
    childProps: [
      {name: 'title', type: 'string', description: 'Step title; numbered when absent.'},
      {name: 'id', type: 'string', structural: true, description: 'Step id for `goTo`; `step-<index>` when absent.'},
      {name: 'done', type: 'bool', structural: true,
        description: 'A terminal report step: no BACK, the closing button reads CLOSE.'},
    ],
    events: ['finish', 'cancel', 'step'],
    acceptsChildren: true,
    defaultChildren: [{tag: 'u2-panel', props: {title: 'Step 1'}}, {tag: 'u2-panel', props: {title: 'Step 2'}}],
    example: {tag: 'u2-wizard', children: [
      {tag: 'u2-panel', props: {title: 'Source'}, children: [{tag: 'u2-text-input', props: {label: 'File'}}]},
      {tag: 'u2-panel', props: {title: 'Options'}, children: [{tag: 'u2-bool-input', props: {label: 'Header row'}}]},
    ]},
  },
  {
    tag: 'u2-accordion',
    category: 'Containers',
    create: () => new Accordion(),
    createWithChildren: (_props, children, nodes) => {
      const accordion = new Accordion();
      for (let i = 0; i < children.length; i++)
        accordion.addPane(childTitle(nodes[i], `Pane ${i + 1}`), element(children[i]), false, childIcon(nodes[i]));
      return accordion;
    },
    setChildProp: (parent, index, name, value) => {
      if (name === 'title')
        (parent as Accordion).setPaneTitle(index, typeof value === 'string' ? value : `Pane ${index + 1}`);
    },
    description: 'Independently expanding panes, one per spec child; the child node carries the pane title.',
    usage: 'In code, `addPane(title, content, {hint, validity, persistKey})`: a folded pane says what it ' +
      'holds (`hint`, a signal or an element for a live summary), shows the problem inside it ' +
      '(`validity`) and remembers the user\'s open or close (`persistKey`; `pane.userSet` until then). ' +
      '`new Accordion({flush: true})` lines the titles and the content up with the captions of a form ' +
      'around it, the chevron at the header\'s end — the panes of a docked editor.',
    props: [],
    childProps: [
      {name: 'title', type: 'string', description: 'Pane title; numbered when absent.'},
      {name: 'icon', type: 'string', inputType: 'Icon', structural: true,
        description: 'Font Awesome icon name shown before the title.'},
    ],
    acceptsChildren: true,
    defaultChildren: [{tag: 'u2-panel', props: {title: 'Pane 1'}}],
    designerActions: [{text: 'Add pane', produce: (node) => ({op: 'add-child',
      node: {tag: 'u2-panel', props: {title: `Pane ${(node.children?.length ?? 0) + 1}`}}})}],
    example: {tag: 'u2-accordion', children: [
      {tag: 'u2-panel', props: {title: 'General'}, children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
      {tag: 'u2-panel', props: {title: 'Advanced'}, children: [{tag: 'u2-bool-input', props: {label: 'Active'}}]},
    ]},
  },
  {
    tag: 'u2-section',
    category: 'Containers',
    create: (props) => section(props, []),
    createWithChildren: (props, children) => section(props, children),
    description: 'One collapsible titled section: a form-category-styled header whose hover ' +
      'chevron sits left of the caption; clicking the header folds the body.',
    usage: 'For grouping form fields or page content under a heading. Inputs inside a ' +
      'u2-section child of a u2-form are adopted by the section, not the form — they do not ' +
      'join the form\'s label column. A managed pane set belongs in u2-accordion. In code, a folded ' +
      'section says what it holds with `hint`, shows the problem inside it with `validity`, and remembers ' +
      'the user\'s open or close with `persistKey`; `form.add(input, section.body)` keeps an input in the ' +
      'form\'s label column.',
    props: [
      {name: 'title', type: 'string', bindable: true},
      {name: 'expanded', type: 'bool', bindable: true, twoWay: true},
      {name: 'collapsible', type: 'bool', structural: true,
        description: 'false renders a plain heading and body — no chevron, no click.'},
    ],
    acceptsChildren: true,
    defaults: {title: 'Section'},
    example: {tag: 'u2-section', props: {title: 'Study'}, children: [
      {tag: 'u2-choice-input', props: {label: 'Country', items: ['USA', 'France']}},
    ]},
  },
  {
    tag: 'u2-page',
    category: 'Containers',
    create: (props) => page(props, []),
    createWithChildren: (props, children) => page(props, children),
    description: 'A view\'s content: a body that alone scrolls, every spec child in it one under the other.',
    usage: 'For the content of a view or a dialog-sized editor. In code, `new Page({header, sections, footer, ' +
      'maxWidth})` adds a header that stays put, titled sections with their actions, and an action row pinned ' +
      'under the body (buttons on the right, a status line beside them); the forms on a page share one label ' +
      'column. Spacing is the `--u2-page-*` properties, set on an element of your own.',
    props: [
      {name: 'maxWidth', type: 'double', structural: true,
        description: 'The body\'s widest, in px; its forms go as wide.'},
    ],
    acceptsChildren: true,
    example: {tag: 'u2-page', props: {maxWidth: 640}, children: [
      {tag: 'u2-text-input', props: {label: 'Name'}},
    ]},
  },
  {
    tag: 'u2-tabs',
    category: 'Containers',
    create: (props) => new TabStrip(tabStrip(props)),
    createWithChildren: (props, children, nodes) => {
      const tabs = new TabStrip(tabStrip(props));
      for (let i = 0; i < children.length; i++) {
        tabs.addTab({id: `tab-${i}`, title: childTitle(nodes[i], `Tab ${i + 1}`), icon: childIcon(nodes[i]),
          content: element(children[i])});
      }
      if (typeof props.activeTab === 'string')
        tabs.activeTab.value = props.activeTab;
      return tabs;
    },
    // a designer rebuild keeps the user's tab; an explicit activeTab — literal or bound — stays the
    // document's to decide
    carry: (node, prev, next) => {
      if (node.props?.activeTab !== undefined || node.bind?.activeTab !== undefined ||
          !(prev instanceof TabStrip) || !(next instanceof TabStrip))
        return;
      const active = prev.activeTab.peek();
      const index = active !== null && active.startsWith('tab-') ? Number(active.slice(4)) : NaN;
      if (index < (node.children?.length ?? 0))
        next.activeTab.value = active;
    },
    setChildProp: (parent, index, name, value) => {
      if (name === 'title')
        (parent as TabStrip).setTabTitle(`tab-${index}`, typeof value === 'string' ? value : `Tab ${index + 1}`);
    },
    description: 'Tab strip with one tab per spec child; the child node carries the tab label.',
    usage: 'In code, give a tab a `validity` signal and it is marked while that says something; a tab ' +
      'whose content is a Form answers with the form\'s shown error. A tab\'s content is placed on its first ' +
      'activation; `eager: true` places element contents at once, hidden, so a form on a tab not yet visited ' +
      'measures its captions (a `LabelColumn` across the tabs).',
    props: [
      {name: 'orientation', type: 'string', choices: ['horizontal', 'vertical'], structural: true,
        description: '"horizontal" (default) or "vertical" — a header column on the left.'},
      {name: 'variant', type: 'string', choices: ['platform', 'document'], structural: true,
        description: '"platform" (default) is the ui.tabControl look; "document" the IDE document-tab skin.'},
      {name: 'activeTab', type: 'string', bindable: true, twoWay: true,
        description: 'Id of the selected tab (spec children get tab-0, tab-1, …).'},
    ],
    childProps: [
      {name: 'title', type: 'string', description: 'Tab label; numbered when absent.'},
      {name: 'icon', type: 'string', inputType: 'Icon', structural: true,
        description: 'Font Awesome icon name shown before the title.'},
    ],
    acceptsChildren: true,
    defaults: {orientation: 'horizontal', variant: 'platform'},
    defaultChildren: [{tag: 'u2-panel', props: {title: 'Tab 1'}}],
    designerActions: [{text: 'Add tab', produce: (node) => ({op: 'add-child',
      node: {tag: 'u2-panel', props: {title: `Tab ${(node.children?.length ?? 0) + 1}`}}})}],
    example: {tag: 'u2-tabs', children: [
      {tag: 'u2-panel', props: {title: 'Data', icon: 'table'},
        children: [{tag: 'u2-text-input', props: {label: 'Table'}}]},
      {tag: 'u2-panel', props: {title: 'Style'}, children: [{tag: 'u2-bool-input', props: {label: 'Legend'}}]},
    ]},
  },
  {
    tag: 'u2-form',
    category: 'Containers',
    create: (props) => new Form({
      layout: props.layout as FormLayout | undefined,
      captionAlign: props.captionAlign as LiveOption<'right' | 'left'> | undefined,
    }),
    adopt: (parent, child) => {
      const form = parent as Form;
      if (child instanceof Input)
        form.add(child);
      else if (child instanceof SpecFormButtons)
        form.addButtons((row) => row.append(child.root), {align: child.align});
      else
        form.addElement(element(child));
    },
    description: 'Vertical input layout: spec children that are inputs join the form and share its ' +
      'label column, validity and Enter navigation; a u2-form-buttons child is its button row; anything else ' +
      'renders as a plain row in place.',
    usage: 'Labels share a left column and values stay left-aligned — never center form content. ' +
      'For forms over Property metadata use `propertyForm`/`objectForm` (u2/dg) instead of hand-building. ' +
      'Forms that must line up — the tabs of one editor — share a `LabelColumn` (`labelColumn` option); ' +
      '`labelWidth` fixes the column. The labels align before the first paint; `layoutNow()` for a test ' +
      'that reads the layout in the same task. `addColumns(inputs, {minWidth})` puts short fields two to ' +
      'a row, each editor filling its column, one to a row where two would go under `minWidth`. ' +
      '`quiet: true` drops the underlines at rest — for an inspector, not a form people fill in. Beside a ' +
      'Save, `statusLine(problem, ' +
      'otherwise)` says what Save will do or why it cannot yet.',
    props: [
      {name: 'layout', type: 'string', choices: ['auto', 'normal', 'wide', 'tall'], structural: true,
        description: '"auto" (default) picks wide or tall by fit — editors fill the row while ' +
          'the label column plus the minimum editor width fit, captions move above when they ' +
          'do not; "normal" never switches; "wide" stretches editors across the row; "tall" ' +
          'puts captions above their editors.'},
      {name: 'captionAlign', type: 'string', choices: ['right', 'left'], bindable: true,
        description: 'Caption alignment in caption-left layouts; right is the platform default.'},
    ],
    acceptsChildren: true,
    defaults: {layout: 'auto'},
    example: {tag: 'u2-form', children: [
      {tag: 'u2-text-input', props: {label: 'Name'}},
      {tag: 'u2-number-input', props: {label: 'Amount', value: 1}},
    ]},
  },
  {
    tag: 'u2-property-grid',
    category: 'Display',
    create: (props) => {
      const grid = new PropertyGrid();
      grid.setProperties((props.properties as PropDescriptor[]) ?? [],
        (props.values as Record<string, unknown>) ?? {});
      return grid;
    },
    description: 'Two-column property editor over a JSON descriptor list.',
    usage: 'For editing a live object that carries Property metadata, prefer `propertyForm` — this ' +
      'grid is for descriptor-driven settings panels.',
    props: [
      {name: 'properties', type: 'object',
        description: 'Row descriptors: name, type (string/int/double/bool/choice/string_list), ' +
          'category, choices, min, max, description, readonly.'},
      {name: 'values', type: 'object', description: 'Initial value per property name.'},
    ],
    example: {tag: 'u2-property-grid', props: {
      properties: [{name: 'Title', type: 'string'}, {name: 'Rows', type: 'int', min: 0}],
      values: {Title: 'Demo', Rows: 20},
    }},
  },
  {
    tag: 'u2-icon',
    category: 'Display',
    create: (props) => new SpecIcon(props),
    description: 'A Font Awesome icon.',
    props: [
      {name: 'name', type: 'string', inputType: 'Icon', structural: true,
        description: 'Font Awesome icon name; read once — `name` is also the node\'s identity accessor.'},
      {name: 'variant', type: 'string', choices: ['light', 'regular', 'solid', 'brands'],
        description: 'Face; light by default, brands picked on its own for a brand name.'},
      {name: 'size', type: 'string', choices: ['small', 'medium', 'large']},
    ],
    defaults: {name: 'star'},
    example: {tag: 'u2-icon', props: {name: 'flask', size: 'large'}},
  },
  {
    tag: 'u2-breadcrumbs',
    category: 'Display',
    create: (props) => new Breadcrumbs({items: (props.items as string[]) ?? []}),
    description: 'Path bar that collapses its middle segments when it overflows.',
    usage: 'For navigation within your own view. The shell already renders the view\'s own ' +
      'breadcrumb from its name/path — don\'t duplicate it.',
    props: [{name: 'items', type: 'string_list'}],
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    events: ['click'],
    designerActions: [{text: 'Add item', produce: (node) => {
      const items = (node.props?.items ?? []) as string[];
      return {op: 'set-prop', name: 'items', value: [...items, `Item ${items.length + 1}`]};
    }}],
    example: {tag: 'u2-breadcrumbs', props: {items: ['Home', 'Projects', 'Demo']}},
  },
  {
    tag: 'u2-slider-input',
    category: 'Inputs',
    create: (props) => new SliderInput({
      ...inputOptions<number | null>(props),
      min: props.min as number | undefined,
      max: props.max as number | undefined,
      step: props.step as number | undefined,
    }),
    description: 'Bare range track; the value shows in a tooltip on hover and while dragging.',
    usage: 'For a value that is only meaningful within known bounds. When the exact number matters ' +
      'too, use `u2-number-input` — it carries the slider as an option.',
    props: inputProps('double',
      {name: 'min', type: 'double'},
      {name: 'max', type: 'double'},
      {name: 'step', type: 'double', description: '(max - min) / 100 by default.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-slider-input', props: {label: 'Opacity', value: 50, min: 0, max: 100}},
  },
  {
    tag: 'u2-range-slider',
    category: 'Inputs',
    create: (props) => new RangeSlider({
      min: props.min as number | undefined,
      max: props.max as number | undefined,
      step: props.step as number | undefined,
      minRange: props.minRange as number | undefined,
      vertical: props.vertical as boolean | undefined,
      lo: props.lo as number | Signal<number> | undefined,
      hi: props.hi as number | Signal<number> | undefined,
      label: props.label as string | undefined,
    }),
    description: 'Two-handle range selector: drag a handle to move one end, the band between ' +
      'them to move the whole window.',
    usage: 'For picking an interval within known bounds — filters, zoom windows, thresholds ' +
      'with a low and a high end. A single number belongs in `u2-slider-input`.',
    props: [
      {name: 'lo', type: 'double', bindable: true, twoWay: true,
        description: 'Lower end of the selected range.'},
      {name: 'hi', type: 'double', bindable: true, twoWay: true,
        description: 'Upper end of the selected range.'},
      {name: 'min', type: 'double'},
      {name: 'max', type: 'double'},
      {name: 'step', type: 'double', description: 'Snap increment for dragging; continuous when absent.'},
      {name: 'minRange', type: 'double', description: 'Smallest allowed hi − lo; 0 by default.'},
      {name: 'vertical', type: 'bool', structural: true, description: 'Bottom-to-top track.'},
      {name: 'label', type: 'string', description: 'Accessible name of the pair (aria-label).'},
    ],
    defaults: {min: 0, max: 100, lo: 20, hi: 60},
    events: ['input', 'change'],
    example: {tag: 'u2-range-slider', props: {label: 'Price range', min: 0, max: 100, lo: 20, hi: 60}},
  },
  {
    tag: 'u2-radio-input',
    category: 'Inputs',
    create: (props) => new RadioInput({
      ...inputOptions<string | null>(props),
      items: (props.items as string[]) ?? [],
      buttons: props.buttons as boolean | undefined,
      itemTooltips: props.itemTooltips as Record<string, string> | undefined,
    }),
    description: 'Radio group over native radios; one item at a time.',
    usage: 'For three to five always-visible options. Above that use `u2-choice-input`.',
    props: inputProps('string',
      {name: 'items', type: 'string_list'},
      {name: 'buttons', type: 'bool', structural: true, description: 'Raised-button variant instead of radio dots.'},
      {name: 'itemTooltips', type: 'object', description: 'Hover text per item.'}),
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    events: ['change'],
    example: {tag: 'u2-radio-input', props: {label: 'Stage', items: ['Discovery', 'Preclinical'],
      value: 'Discovery'}},
  },
  {
    tag: 'u2-color-input',
    category: 'Inputs',
    create: (props) => new ColorInput({
      ...inputOptions<string>(props),
      swatchOnly: props.swatchOnly as boolean | undefined,
    }),
    description: 'Hex field plus a swatch that opens the browser color picker.',
    props: inputProps('string',
      {name: 'swatchOnly', type: 'bool', structural: true, description: 'Swatch alone, no hex field.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-color-input', props: {label: 'Color', value: '#1f77b4'}},
  },
  {
    tag: 'u2-icon-input',
    category: 'Inputs',
    create: (props) => new IconInput(inputOptions<string>(props)),
    description: 'Font Awesome icon name, picked from a searchable grid of every icon the platform ships.',
    props: inputProps('string').map((p) => p.name === 'value' ? {...p, inputType: 'Icon'} : p),
    events: ['input', 'change'],
    example: {tag: 'u2-icon-input', props: {label: 'Icon', value: 'star'}},
  },
  {
    tag: 'u2-list-input',
    category: 'Inputs',
    create: (props) => new ListInput({
      ...inputOptions<string[]>(props),
      placeholder: props.placeholder as LiveOption<string> | undefined,
    }),
    description: 'Comma-separated list; items carrying a comma are quoted. Expands to a textarea.',
    props: inputProps('string_list', {name: 'placeholder', type: 'string', bindable: true}),
    events: ['input', 'change'],
    example: {tag: 'u2-list-input', props: {label: 'Synonyms', value: ['ASA', '2-acetoxybenzoic acid']}},
  },
  {
    tag: 'u2-map-input',
    category: 'Inputs',
    create: (props) => new MapInput(inputOptions<Record<string, string>>(props)),
    description: 'Key/value rows with per-row add and remove; duplicate keys are invalid, and a ' +
      'row with an empty key stays out of the value.',
    props: inputProps('object'),
    events: ['input', 'change'],
    example: {tag: 'u2-map-input', props: {label: 'Params', value: {solvent: 'DMSO'}}},
  },
  {
    tag: 'u2-records-input',
    category: 'Inputs',
    create: (props) => new RecordsInput<Record<string, unknown>>({
      ...inputOptions<Record<string, unknown>[]>(props),
      columns: itemList(props.fields).map((name) => ({name,
        input: (value) => new TextInput({inline: true, value: (value as string | undefined) ?? ''})})),
      addText: props.addText as string | undefined,
      orderable: props.orderable as boolean | undefined,
    }),
    description: 'An editable list of records: a row each, an inline input per field under a header, a ' +
      'remove button per row and an add link; orderable rows move by a handle or Alt+↑/↓. The value is the ' +
      'records; a failing cell names its row in the verdict ("Row 2: Name: Value can\'t be empty").',
    usage: 'For a short list of typed rows edited in place — a template\'s fields, its stages. In code each ' +
      'column gives the input its cells get (`{name, header, input: (value) => new ChoiceInput({inline: ' +
      'true, …}), width}`); a spec gives text fields by name. Validate the list as a whole with ' +
      '`addValidator` (duplicate names).',
    props: inputProps('object',
      {name: 'fields', type: 'string_list', structural: true, description: 'The columns: text fields, in order.'},
      {name: 'addText', type: 'string', structural: true, description: 'The add link\'s text; Add by default.'},
      {name: 'orderable', type: 'bool', structural: true, description: 'Rows move by a handle or Alt+↑/↓.'}),
    defaults: {fields: ['name']},
    events: ['change'],
    example: {tag: 'u2-records-input', props: {label: 'Stages', fields: ['name'], addText: 'Add stage',
      orderable: true, value: [{name: 'Design'}, {name: 'Synthesis'}, {name: 'Assay'}]}},
  },
  {
    tag: 'u2-checklist-input',
    category: 'Inputs',
    create: (props) => {
      const input = new ChecklistInput({
        ...inputOptions<string[]>(props),
        items: itemList(props.items),
        orderable: props.orderable as boolean | undefined,
      });
      boundItems(input, props.items);
      return input;
    },
    description: 'Checkboxes the user can also put in order: the value is the checked items in the order ' +
      'shown, `order` every item\'s.',
    usage: 'For picking a subset where the order matters too — the functions to run and their order, the ' +
      'columns to show. Without `orderable` it is a plain check list (`u2-multi-choice-input` adds a ' +
      'summary checkbox).',
    props: inputProps('string_list', {name: 'items', type: 'string_list', bindable: true},
      {name: 'orderable', type: 'bool', structural: true, description: 'Rows move by a handle or Alt+↑/↓.'}),
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    events: ['change'],
    example: {tag: 'u2-checklist-input', props: {label: 'Run', items: ['Descriptors', 'Toxicity', 'Solubility'],
      value: ['Toxicity', 'Descriptors'], orderable: true}},
  },
  {
    tag: 'u2-chips-input',
    category: 'Inputs',
    create: (props) => {
      const input = new ChipsInput({
        ...inputOptions<string[]>(props),
        items: itemList(props.items),
        emptyText: props.emptyText as string | undefined,
      });
      boundItems(input, props.items);
      return input;
    },
    description: 'Toggle chips, one per item; the value holds the pressed ones in item order.',
    usage: 'For a handful of names picked inline (the groups a column is visible to). A long list ' +
      'belongs in `u2-multi-select`, a tall one in `u2-multi-choice-input`.',
    props: inputProps('string_list', {name: 'items', type: 'string_list', bindable: true},
      {name: 'emptyText', type: 'string', description: 'Stands in for an empty list; \'No items\' by default.'}),
    defaults: {items: ['Item 1', 'Item 2', 'Item 3']},
    events: ['change'],
    example: {tag: 'u2-chips-input', props: {label: 'Visible to', items: ['Sales', 'Chemists'], value: ['Sales']}},
  },
  {
    tag: 'u2-access-grid',
    category: 'Inputs',
    create: (props) => new AccessGrid({
      ...inputOptions<AccessRow[]>(props),
      capabilities: itemList(props.capabilities),
      principals: props.principals instanceof Signal ? props.principals as Signal<string[]> :
        itemList(props.principals),
      principalLabel: props.principalLabel as string | undefined,
      addText: props.addText as string | undefined,
    }),
    description: 'Who may do what: a row per principal, a checkbox per capability, a remove button, and ' +
      'a picker that adds a principal. The value is the rows — `[{principal, can: {view: true, …}}]`.',
    usage: 'For the access a schema or a table grants (groups × View / Edit / Delete). Inherited rows ' +
      'and locked capabilities are code options; a spec passes the columns and the picker\'s list.',
    props: inputProps('object',
      {name: 'capabilities', type: 'string_list', structural: true, description: 'The columns, in order.'},
      {name: 'principals', type: 'string_list', bindable: true, description: 'What the add picker offers.'},
      {name: 'principalLabel', type: 'string', description: 'The first column\'s header; Group by default.'},
      {name: 'addText', type: 'string', description: 'The picker\'s placeholder.'}),
    defaults: {capabilities: ['view', 'edit', 'delete'], principals: ['All users', 'Developers']},
    events: ['change'],
    example: {tag: 'u2-access-grid', props: {label: 'Access', capabilities: ['view', 'edit', 'delete'],
      principals: ['Sales', 'Developers'], value: [{principal: 'Sales', can: {view: true}}]}},
  },
  {
    tag: 'u2-qnum-input',
    category: 'Inputs',
    create: (props) => new QNumInput({
      ...inputOptions<number | null>(props),
      placeholder: props.placeholder as LiveOption<string> | undefined,
    }),
    description: 'Qualified number editor — "<5.2", "10", ">1e3". The value is the packed double ' +
      'the platform stores for a qnum: the qualifier rides in the two lowest mantissa bits.',
    usage: 'For measurements reported against a detection limit. A plain number belongs in ' +
      '`u2-number-input` — this editor spends two bits of precision on the qualifier.',
    props: inputProps('double', {name: 'placeholder', type: 'string', bindable: true}),
    events: ['input', 'change'],
    example: {tag: 'u2-qnum-input', props: {label: 'IC50', value: QNum.less(5.2)}},
  },
  {
    tag: 'u2-font-input',
    category: 'Inputs',
    create: (props) => new FontInput({
      ...inputOptions<string>(props),
      families: props.families as string[] | undefined,
      sizes: props.sizes as number[] | undefined,
    }),
    description: 'Font editor over the platform font string — size box with presets, bold and ' +
      'italic toggles, family select.',
    usage: 'For the "<weight> <style> <size>px <family>" strings viewer looks carry; anything ' +
      'else that needs a font belongs in CSS.',
    props: inputProps('string',
      {name: 'families', type: 'string_list', structural: true,
        description: 'Families on offer; four defaults otherwise.'},
      {name: 'sizes', type: 'object', structural: true, description: 'Sizes in the popup under the size box.'}),
    events: ['change'],
    example: {tag: 'u2-font-input', props: {label: 'Font', value: 'bold normal 14px "Roboto"'}},
  },
  {
    tag: 'u2-image-input',
    category: 'Inputs',
    create: (props) => new ImageInput({
      ...inputOptions<string>(props),
      alt: props.alt as string | undefined,
    }),
    description: 'Image preview by URL or data URL; the REMOVE link appears on hover and clears it.',
    props: inputProps('string', {name: 'alt', type: 'string', description: 'Alternative text.'}),
    events: ['change'],
    example: {tag: 'u2-image-input', props: {label: 'Logo', value: 'https://datagrok.ai/img/logo.svg'}},
  },
  {
    tag: 'u2-card',
    category: 'Containers',
    create: (props) => card(props, []),
    createWithChildren: (props, children) => card(props, children),
    description: 'Surface with optional header (title, subtitle, icon, actions), media image, body and footer.',
    usage: 'For dashboard tiles and item previews; a card grid is a collection of these. ' +
      'A KPI number belongs in `u2-stat-card`; a plain padded area in a panel.',
    props: [
      {name: 'title', type: 'string', bindable: true},
      {name: 'subtitle', type: 'string', bindable: true},
      {name: 'icon', type: 'string', inputType: 'Icon', structural: true},
      {name: 'media', type: 'string', structural: true, description: 'Image URL shown above the body.'},
      {name: 'clickable', type: 'bool', structural: true, description: 'The whole card acts as a button.'},
      {name: 'selectable', type: 'bool', structural: true, description: 'Click toggles the selection ring.'},
      {name: 'selected', type: 'bool', bindable: true, twoWay: true},
    ],
    acceptsChildren: true,
    events: ['click'],
    defaults: {title: 'Card'},
    example: {tag: 'u2-card', props: {title: 'Aspirin', subtitle: 'NSAID'},
      children: [{tag: 'p', props: {text: 'Acetylsalicylic acid, 180.16 g/mol.'}}]},
  },
  {
    tag: 'u2-stat-card',
    category: 'Display',
    create: (props) => new StatCard({
      label: (props.label as string) ?? '',
      value: props.value as string | Signal<string | number | undefined> | undefined,
      delta: props.delta as number | Signal<number | undefined> | undefined,
      deltaInverted: props.deltaInverted as boolean | undefined,
      icon: props.icon as string | undefined,
    }),
    description: 'KPI card: large formatted value, label, optional up/down delta and icon.',
    usage: 'For dashboard metrics. Bind `value` to a source output for live numbers; ' +
      'loading states come from the source, not from this card.',
    props: [
      {name: 'label', type: 'string', bindable: true},
      {name: 'value', type: 'string', bindable: true,
        description: 'Literal values are strings ("1.2M"); a bound signal may hold a number, formatted by the card.'},
      {name: 'delta', type: 'double', bindable: true,
        description: 'Fractional change; 0.12 renders "+12%" in the success color.'},
      {name: 'deltaInverted', type: 'bool', structural: true, description: 'Negative is good (error rates).'},
      {name: 'icon', type: 'string', inputType: 'Icon', structural: true},
    ],
    defaults: {label: 'Metric', value: '0'},
    example: {tag: 'u2-stat-card', props: {label: 'Revenue', value: '1.2M', delta: 0.12, icon: 'chart-line'}},
  },
  {
    tag: 'u2-app-header',
    category: 'Display',
    create: (props) => new AppHeader({
      icon: (props.icon as string | undefined) ?? 'rocket',
      iconSize: props.iconSize as number | undefined,
      title: props.title as Text | undefined,
      subtitle: props.subtitle as Text | undefined,
      items: itemList(props.items),
      learnMore: props.learnMore as string | undefined,
      bottomLine: props.bottomLine as boolean | undefined,
    }),
    description: 'An app home\'s hero: a big icon beside the bullets of what the app does, and a Learn more link.',
    usage: 'At the top of an app\'s home view, above what the user starts from. A titled header puts the ' +
      'app name and a one-line subtitle beside a smaller icon; an item may be an element or a link ' +
      '(`AppHeader` in code).',
    props: [
      {name: 'icon', type: 'string', inputType: 'Icon', structural: true,
        description: 'An image URL or a Font Awesome name.'},
      {name: 'iconSize', type: 'int', structural: true, description: 'In pixels: 150, or 90 beside a title.'},
      {name: 'title', type: 'string', bindable: true},
      {name: 'subtitle', type: 'string', bindable: true},
      {name: 'items', type: 'string_list', structural: true, description: 'The bullets, in order.'},
      {name: 'learnMore', type: 'string', structural: true, description: 'The URL Learn more opens in a new tab.'},
      {name: 'bottomLine', type: 'bool', structural: true, description: 'A rule under the header.'},
    ],
    defaults: {icon: 'rocket', items: ['What the app does']},
    example: {tag: 'u2-app-header', props: {icon: 'flask',
      items: ['Sketch molecules in the spreadsheet', 'Calculate molecular properties'],
      learnMore: 'https://datagrok.ai/help'}},
  },
  {
    tag: 'u2-data-table',
    category: 'Display',
    create: (props) => {
      const table = new DataTable<Record<string, unknown>>({
        columns: dataTableColumns(props.columns),
        rowHeight: props.rowHeight as number | undefined,
      });
      const items = props.items;
      if (items instanceof Signal || Array.isArray(items))
        table.setItems(items as Signal<Record<string, unknown>[]> | Record<string, unknown>[]);
      return table;
    },
    description: 'Virtualized table: pooled rows under a sticky header, one column per named ' +
      'property, only the visible window in the DOM.',
    usage: 'For a list of records with several columns — bind `items` to a source\'s rows ' +
      '(`$.orders.rows`). Prefer `u2-domain-grid` for editing a domain table, and `u2-domain-list` ' +
      'for one line or a card per row; a handful of rows with rich cells belongs in a plain table.',
    props: [
      {name: 'columns', type: 'string_list',
        description: 'The item properties to show, in order; `name:Header` gives a column its own header.'},
      {name: 'items', type: 'object', bindable: true,
        description: 'The records to show — bind a source\'s rows (`$.substances.rows`).'},
      {name: 'rowHeight', type: 'int', structural: true, description: 'Row and header height in pixels (default 24).'},
    ],
    defaults: {columns: ['name', 'value']},
    example: {tag: 'u2-data-table', props: {columns: ['name', 'cas:CAS number'],
      items: [{name: 'Aspirin', cas: '50-78-2'}, {name: 'Caffeine', cas: '58-08-2'}]}},
  },
  {
    tag: 'u2-catalog-table',
    category: 'Display',
    create: (props) => {
      const columns = dataTableColumns(props.columns).map((c, i): CatalogColumn<Record<string, unknown>> =>
        ({name: c.name, header: c.header, value: (r) => r[c.name], pinned: i === 0}));
      const key = (props.key as string | undefined) ?? columns[0]?.name ?? 'name';
      return new CatalogTable<Record<string, unknown>>({
        items: (props.items instanceof Signal || Array.isArray(props.items) ? props.items : []) as
          Record<string, unknown>[] | Signal<Record<string, unknown>[]>,
        keyOf: (r) => String(r[key] ?? ''),
        columns,
        groupings: itemList(props.groupings).map((name) =>
          ({name, label: columns.find((c) => c.name === name)?.header ?? name, key: (r) => r[name]})),
        persist: props.persistKey as string | undefined,
      });
    },
    description: 'A table of records people come back to: double-click a header to sort, group the rows into ' +
      'foldable sections with counts, choose the columns; the choices are kept.',
    usage: 'For a catalog of records with a name and a few attributes (campaigns, profiles, templates) where the ' +
      'whole list fits the page. In code, `CatalogTable` takes columns with `value`/`render`, groupings, row ' +
      '`actions`, an `AsyncValue` of items and a `persist` adapter; `chrome()` gives the columns, group and refresh ' +
      'icons for the heading. Thousands of rows belong in `u2-data-table`. Recipe: docs/recipes/catalog-table.md.',
    props: [
      {name: 'columns', type: 'string_list', structural: true,
        description: 'The record properties to show, in order; `name:Header` gives one its own header. The first is ' +
          'the record\'s name: always shown, first.'},
      {name: 'items', type: 'object', bindable: true,
        description: 'The records — bind a source\'s rows (`$.campaigns.rows`).'},
      {name: 'key', type: 'string', structural: true,
        description: 'The property that identifies a record; the first column by default.'},
      {name: 'groupings', type: 'string_list', structural: true,
        description: 'The columns the group menu offers.'},
      {name: 'persistKey', type: 'string', structural: true,
        description: 'Keeps the sort, grouping, columns and folds in the browser under this key.'},
    ],
    defaults: {columns: ['name:Name', 'status:Status']},
    example: {tag: 'u2-catalog-table', props: {columns: ['name:Name', 'status:Status', 'molecules:Molecules'],
      groupings: ['status'], items: [{name: 'Kinase hits', status: 'In progress', molecules: 120},
        {name: 'GPCR series', status: 'Submitted', molecules: 48},
        {name: 'Fragments', status: 'In progress', molecules: 9}]}},
  },
  {
    tag: 'u2-filter-builder',
    category: 'Inputs',
    create: (props) => {
      const options = inputOptions<FilterGroup>(props);
      const template = filterTemplate(props.template);
      if (template && options.value === undefined && options.bind === undefined)
        options.value = Filters.applyTemplate(template);
      return new FilterBuilder({
        ...options,
        schema: filterSchema(props.schema),
        mode: props.mode as FilterMode | Signal<FilterMode> | undefined,
        orientation: props.orientation === 'horizontal' ? 'horizontal' : 'vertical',
        template,
        showQuery: props.showQuery as boolean | undefined,
        target: props.target as FilterTarget | undefined,
      });
    },
    description: 'Schema-driven query builder: condition rows (property, operator, value) joined by ' +
      'one and/or toggle, "+" adds and "−" removes; the value is a filter tree, `query` its string.',
    usage: 'Use for a user-editable condition over a schema (columns, entity properties, a domain ' +
      'table). Emit the filter as `value` (`{op: "and", nodes: [{property, operator, value}]}`) or ' +
      'read `query`, the canonical string. Grammar: `prop op value` joined by `and`/`or`, parentheses ' +
      'and `not (…)` for groups; ops `= != > >= < <= like !like starts ends matches in not in between`, ' +
      '`= null` / `!= null` for empty; values: numbers, "strings", true/false/null, ISO dates, relative ' +
      'spans `-1w 2d now`. Examples: `age > 30 and sex = "F"`; `created > -1w and (status in ("Open", ' +
      '"Blocked") or owner = @current)`; `mw between 200 and 500`. Prefer `u2-filter-query-input` for a ' +
      'one-line search box; this control for visible criteria rows.',
    props: inputProps('object',
      {name: 'schema', type: 'object', structural: true,
        description: '`{properties: IProperty[], values?: {<name>: scalar[]}}` — what can be filtered on.'},
      {name: 'mode', type: 'string', choices: ['simple', 'advanced'], bindable: true, twoWay: true,
        description: 'simple: one flat and/or group; advanced: nested groups, not, drag handles.'},
      {name: 'orientation', type: 'string', choices: ['vertical', 'horizontal'], structural: true,
        description: 'horizontal lays the rows out inline with a connector chip between them.'},
      {name: 'template', type: 'object', structural: true,
        description: '`{root, allowedProperties?, allowedOperators?, allowAdvanced?, allowAdd?}` — locks ' +
          'and limits: `allowAdvanced: false` keeps the builder simple, `allowAdd: false` fixes the rows.'},
      {name: 'showQuery', type: 'bool', structural: true,
        description: 'Footer line with the canonical query string; off by default.'},
      {name: 'target', type: 'string', choices: ['domain', 'dataframe'], structural: true,
        description: 'What the tree must be expressible for; validation only.'},
      {name: 'query', type: 'string', bindable: true, description: 'The canonical query string (read-only).'}),
    events: ['change'],
    example: {tag: 'u2-filter-builder', props: {
      label: 'Criteria',
      schema: {properties: [{name: 'name', type: 'string'}, {name: 'age', type: 'int', min: 0, max: 120},
        {name: 'sex', type: 'string', choices: ['F', 'M']}]},
      value: {op: 'and', nodes: [{property: 'age', operator: '>', value: 30}]},
      showQuery: true,
    }},
  },
  {
    tag: 'u2-filter-query-input',
    category: 'Inputs',
    create: (props) => new FilterQueryInput({
      ...inputOptions<FilterGroup>(props),
      schema: filterSchema(props.schema),
      placeholder: props.placeholder as string | undefined,
      target: props.target as FilterTarget | undefined,
    }),
    description: 'One-line query box over a schema with completion at the caret (properties, operators, ' +
      'values, and/or); Enter or blur parses the text into the filter tree, `query` is the text.',
    usage: 'Use for a compact search box over a schema (columns, entity properties, a domain table) ' +
      'where the user types the condition. Emit the filter as `value` (`{op: "and", nodes: [{property, ' +
      'operator, value}]}`) or read `query`, the text. Grammar: `prop op value` joined by `and`/`or`, ' +
      'parentheses and `not (…)` for groups; ops `= != > >= < <= like !like starts ends matches in not in ' +
      'between`, `= null` / `!= null` for empty; values: numbers, "strings", true/false/null, ISO dates, ' +
      'relative spans `-1w 2d now`. Examples: `age > 30 and sex = "F"`; `created > -1w and (status in ' +
      '("Open", "Blocked") or owner = @current)`; `mw between 200 and 500`. Prefer `u2-filter-builder` ' +
      'for visible criteria rows; this control for a one-line search box.',
    props: inputProps('object',
      {name: 'schema', type: 'object', structural: true,
        description: '`{properties: IProperty[], values?: {<name>: scalar[]}}` — what can be filtered on.'},
      {name: 'placeholder', type: 'string'},
      {name: 'target', type: 'string', choices: ['domain', 'dataframe'], structural: true,
        description: 'What the tree must be expressible for; validation only.'},
      {name: 'query', type: 'string', bindable: true,
        description: 'The draft text (per keystroke, may be unparsable; read-only).'}),
    events: ['change'],
    example: {tag: 'u2-filter-query-input', props: {
      label: 'Search',
      schema: {properties: [{name: 'name', type: 'string'}, {name: 'age', type: 'int', min: 0, max: 120},
        {name: 'sex', type: 'string', choices: ['F', 'M']}]},
      value: {op: 'and', nodes: [{property: 'age', operator: '>', value: 30}]},
      placeholder: 'age > 30 and sex = "F"',
    }},
  },
  {
    tag: 'u2-state',
    category: 'Data',
    visual: false,
    createComponent: (props) => new StateSource({type: props.type as string, initial: props.initial}),
    description: 'A value the form owns: one writable signal several controls share.',
    usage: 'For a draft, a flag or a selection that outlives no single control — anything two ' +
      'nodes must agree on. Data that comes from somewhere belongs in a source, not here.',
    props: [
      {name: 'type', type: 'string', structural: true,
        choices: ['string', 'int', 'double', 'bool', 'string_list', 'object'],
        description: 'What the value holds; `initial` is checked against it.'},
      {name: 'initial', type: 'object', structural: true, description: 'The value it starts at, and what ' +
        '`clear` returns it to; the type\'s empty value otherwise.'},
    ],
    defaults: {type: 'string'},
    example: {tag: 'u2-state', name: 'draft', props: {type: 'string', initial: ''}},
  },
];

/** Idempotent per registry — the gallery, the manifest script and tests all call it freely. */
export function registerAll(reg: Registry = globalRegistry): void {
  if (reg.builtins)
    return;
  reg.builtins = true;
  for (const meta of METAS)
    reg.register(meta);
  registerDataSources(reg);
}
