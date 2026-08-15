/* Every component a dg-ui/1 spec can build, with the prop metadata that validates specs and feeds
   manifest.json. Left out: components whose options are functions (list, tree, combobox, dialog,
   toolbar, async view). */
import {Signal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {button, divH, divV, panel} from '../core/elements.js';
import {Input, InputOptions} from '../core/input-base.js';
import {TextInput, TextArea} from '../components/text-input.js';
import {BoolInput} from '../components/bool-input.js';
import {NumberInput} from '../components/number-input.js';
import {ChoiceInput, MultiChoiceInput} from '../components/choice-input.js';
import {Form} from '../components/form.js';
import {Splitter} from '../components/splitter.js';
import {Accordion} from '../components/accordion.js';
import {TabStrip} from '../components/tabs.js';
import {PropertyGrid, PropDescriptor} from '../components/property-grid.js';
import {Breadcrumbs} from '../components/breadcrumbs.js';
import {ComponentMeta, PropMeta, Registry, registry as globalRegistry} from './registry.js';
import type {SpecNode} from './spec.js';

type Props = Record<string, unknown>;
type Child = Component | HTMLElement;

/** A bound `value` arrives as the context signal itself, so the input adopts it instead of copying:
 * no bridge, no echo, and the spec's two-way contract holds through every edit. */
function inputOptions<T>(props: Props): InputOptions<T> {
  const value = props.value;
  const bound = value instanceof Signal;
  return {
    label: props.label as string | undefined,
    name: props.name as string | undefined,
    inline: props.inline as boolean | undefined,
    tooltipText: props.tooltipText as string | undefined,
    value: bound ? undefined : value as T,
    bind: bound ? value as Signal<T> : undefined,
  };
}

function inputProps(value: PropMeta['type'], ...extra: PropMeta[]): PropMeta[] {
  return [
    {name: 'label', type: 'string'},
    {name: 'name', type: 'string', description: 'Stable key for forms and dumps; defaults to the label.'},
    {name: 'value', type: value, bindable: true, twoWay: true},
    {name: 'inline', type: 'bool', description: 'Compact one-row variant without a label.'},
    {name: 'tooltipText', type: 'string'},
    ...extra,
  ];
}

function element(child: Child): HTMLElement {
  return child instanceof Component ? child.root : child;
}

/** The title a container reads off its child's spec node — the one prop a parent owns. */
function childTitle(node: SpecNode | undefined, fallback: string): string {
  const title = node?.props?.title;
  return typeof title === 'string' ? title : fallback;
}

function splitter(props: Props, children: Child[]): Splitter {
  return new Splitter(children.map(element), {
    direction: props.direction === 'vertical' ? 'vertical' : 'horizontal',
    sizes: props.sizes as number[] | undefined,
    minSize: props.minSize as number | undefined,
  });
}

function container(tag: string, create: () => HTMLElement, description: string): ComponentMeta {
  return {
    tag,
    create: (props) => {
      const el = create();
      const cls = props.cls as string | undefined;
      if (cls)
        el.classList.add(...cls.split(' ').filter((c) => c));
      return new Component(el);
    },
    description,
    props: [{name: 'cls', type: 'string', description: 'Extra layout class.'}],
    acceptsChildren: true,
    example: {tag, children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
  };
}

const METAS: ComponentMeta[] = [
  {
    tag: 'u2-text-input',
    create: (props) => new TextInput({
      ...inputOptions<string>(props),
      placeholder: props.placeholder as string | undefined,
      search: props.search as boolean | undefined,
    }),
    description: 'Single-line text editor with validation.',
    usage: 'With `search: true` this is the filter box; in a platform view the main filter belongs ' +
      'in the view ribbon (`appView({ribbon})`), not inside the content area.',
    props: inputProps('string',
      {name: 'placeholder', type: 'string'},
      {name: 'search', type: 'bool', description: 'Magnifier affordance plus a clear button.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-text-input', props: {label: 'Name', value: 'Aspirin'}},
  },
  {
    tag: 'u2-text-area',
    create: (props) => new TextArea({
      ...inputOptions<string>(props),
      placeholder: props.placeholder as string | undefined,
      autoGrow: props.autoGrow as boolean | undefined,
    }),
    description: 'Multi-line text editor.',
    props: inputProps('string',
      {name: 'placeholder', type: 'string'},
      {name: 'autoGrow', type: 'bool', description: 'Grows with the text instead of scrolling.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-text-area', props: {label: 'Description', placeholder: 'Up to 40 characters'}},
  },
  {
    tag: 'u2-bool-input',
    create: (props) => new BoolInput({
      ...inputOptions<boolean>(props),
      switch: props.switch as boolean | undefined,
    }),
    description: 'Checkbox or switch.',
    props: inputProps('bool', {name: 'switch', type: 'bool', description: 'Compact toggle instead of a checkbox.'}),
    events: ['change'],
    example: {tag: 'u2-bool-input', props: {label: 'Active', value: true, switch: true}},
  },
  {
    tag: 'u2-number-input',
    create: (props) => new NumberInput({
      ...inputOptions<number | null>(props),
      mode: props.mode as 'int' | 'float' | undefined,
      min: props.min as number | undefined,
      max: props.max as number | undefined,
      step: props.step as number | undefined,
    }),
    description: 'Numeric editor with a spinner; text that does not parse only marks the input invalid.',
    props: inputProps('float',
      {name: 'mode', type: 'string', description: '"int" or "float" (default).'},
      {name: 'min', type: 'float'},
      {name: 'max', type: 'float'},
      {name: 'step', type: 'float', description: 'Spinner and arrow-key increment; 1 by default.'}),
    events: ['input', 'change'],
    example: {tag: 'u2-number-input', props: {label: 'Amount', value: 10, min: 0, max: 100}},
  },
  {
    tag: 'u2-choice-input',
    create: (props) => new ChoiceInput({
      ...inputOptions<string | null>(props),
      items: (props.items as string[]) ?? [],
      nullable: props.nullable as boolean | undefined,
    }),
    description: 'Single choice over a native select.',
    props: inputProps('string',
      {name: 'items', type: 'string[]'},
      {name: 'nullable', type: 'bool', description: 'Offers an empty option; true by default.'}),
    events: ['change'],
    example: {tag: 'u2-choice-input', props: {label: 'Series', items: ['a', 'b', 'c'], value: 'b'}},
  },
  {
    tag: 'u2-multi-choice-input',
    create: (props) => new MultiChoiceInput({
      ...inputOptions<string[]>(props),
      items: (props.items as string[]) ?? [],
    }),
    description: 'Checkbox list; the value holds the checked items in item order.',
    props: inputProps('string[]', {name: 'items', type: 'string[]'}),
    events: ['change'],
    example: {tag: 'u2-multi-choice-input', props: {label: 'Tags', items: ['acid', 'base'], value: ['acid']}},
  },
  {
    tag: 'u2-button',
    create: (props) => new Component(button((props.text as string | undefined) ?? '', () => {},
      {primary: props.primary === true})),
    description: 'Push button; `on: {"click": "cmd:<name>"}` runs a context command.',
    props: [
      {name: 'text', type: 'string'},
      {name: 'primary', type: 'bool', description: 'Accent styling for the main action.'},
    ],
    events: ['click'],
    example: {tag: 'u2-button', props: {text: 'Save', primary: true}, on: {click: 'cmd:save'}},
  },
  container('u2-panel', panel, 'Padded vertical container.'),
  container('u2-div-v', divV, 'Vertical stack.'),
  container('u2-div-h', divH, 'Horizontal row.'),
  {
    tag: 'u2-splitter',
    create: (props) => splitter(props, []),
    createWithChildren: (props, children) => splitter(props, children),
    description: 'Resizable panels separated by draggable sashes; every spec child becomes a panel.',
    usage: 'Master-detail layouts (list left, details right): see the entity-browser recipe. In a ' +
      'platform view, consider handing details to the context panel (`grok.shell.o`) instead of a pane.',
    props: [
      {name: 'direction', type: 'string', description: '"horizontal" (default) or "vertical".'},
      {name: 'sizes', type: 'json', description: 'Fraction per panel, normalized; equal shares by default.'},
      {name: 'minSize', type: 'float', description: 'Smallest panel size in pixels; 60 by default.'},
    ],
    acceptsChildren: true,
    example: {tag: 'u2-splitter', props: {direction: 'horizontal', sizes: [0.3, 0.7]}, children: [
      {tag: 'u2-panel', children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
      {tag: 'u2-panel', children: [{tag: 'u2-text-area', props: {label: 'Notes'}}]},
    ]},
  },
  {
    tag: 'u2-accordion',
    create: () => new Accordion(),
    createWithChildren: (_props, children, nodes) => {
      const accordion = new Accordion();
      for (let i = 0; i < children.length; i++)
        accordion.addPane(childTitle(nodes[i], `Pane ${i + 1}`), element(children[i]));
      return accordion;
    },
    description: 'Independently expanding panes, one per spec child; the child node carries the pane title.',
    props: [],
    childProps: [{name: 'title', type: 'string', description: 'Pane title; numbered when absent.'}],
    acceptsChildren: true,
    example: {tag: 'u2-accordion', children: [
      {tag: 'u2-panel', props: {title: 'General'}, children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
      {tag: 'u2-panel', props: {title: 'Advanced'}, children: [{tag: 'u2-bool-input', props: {label: 'Active'}}]},
    ]},
  },
  {
    tag: 'u2-tabs',
    create: () => new TabStrip(),
    createWithChildren: (_props, children, nodes) => {
      const tabs = new TabStrip();
      for (let i = 0; i < children.length; i++)
        tabs.addTab({id: `tab-${i}`, label: childTitle(nodes[i], `Tab ${i + 1}`), content: element(children[i])});
      return tabs;
    },
    description: 'Tab strip with one tab per spec child; the child node carries the tab label.',
    props: [],
    childProps: [{name: 'title', type: 'string', description: 'Tab label; numbered when absent.'}],
    acceptsChildren: true,
    example: {tag: 'u2-tabs', children: [
      {tag: 'u2-panel', props: {title: 'Data'}, children: [{tag: 'u2-text-input', props: {label: 'Table'}}]},
      {tag: 'u2-panel', props: {title: 'Style'}, children: [{tag: 'u2-bool-input', props: {label: 'Legend'}}]},
    ]},
  },
  {
    tag: 'u2-form',
    create: (props) => new Form({
      condensed: props.condensed as boolean | undefined,
      wide: props.wide as boolean | undefined,
    }),
    adopt: (parent, child, index) => {
      const form = parent as Form;
      if (child instanceof Input) {
        form.add(child);
        return;
      }
      console.warn(`u2 spec: u2-form: child ${index} is not an input — appended to the form root`);
      form.root.append(element(child));
    },
    description: 'Vertical input layout: spec children that are inputs join the form and share its ' +
      'label column, validity and Enter navigation; anything else lands in the form root.',
    usage: 'Labels share a left column and values stay left-aligned — never center form content. ' +
      'For forms over Property metadata use `propertyForm`/`objectForm` (u2/dg) instead of hand-building.',
    props: [
      {name: 'condensed', type: 'bool', description: 'Tighter rows.'},
      {name: 'wide', type: 'bool', description: 'Two columns.'},
    ],
    acceptsChildren: true,
    example: {tag: 'u2-form', children: [
      {tag: 'u2-text-input', props: {label: 'Name'}},
      {tag: 'u2-number-input', props: {label: 'Amount', value: 1}},
    ]},
  },
  {
    tag: 'u2-property-grid',
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
      {name: 'properties', type: 'json',
        description: 'Row descriptors: name, type (string/int/float/bool/choice), category, choices, ' +
          'min, max, description, readonly.'},
      {name: 'values', type: 'json', description: 'Initial value per property name.'},
    ],
    example: {tag: 'u2-property-grid', props: {
      properties: [{name: 'Title', type: 'string'}, {name: 'Rows', type: 'int', min: 0}],
      values: {Title: 'Demo', Rows: 20},
    }},
  },
  {
    tag: 'u2-breadcrumbs',
    create: (props) => new Breadcrumbs({items: (props.items as string[]) ?? []}),
    description: 'Path bar that collapses its middle segments when it overflows.',
    usage: 'For navigation within your own view. The shell already renders the view\'s own ' +
      'breadcrumb from its name/path — don\'t duplicate it.',
    props: [{name: 'items', type: 'string[]'}],
    events: ['click'],
    example: {tag: 'u2-breadcrumbs', props: {items: ['Home', 'Projects', 'Demo']}},
  },
];

const registered = new WeakSet<Registry>();

/** Idempotent per registry — the gallery, the manifest script and tests all call it freely. */
export function registerAll(reg: Registry = globalRegistry): void {
  if (registered.has(reg))
    return;
  registered.add(reg);
  for (const meta of METAS)
    reg.register(meta);
}
