/* The Open menu's sample forms (F12): five conceptually distinct specs, data-only and
   headless-importable. Every one parses and renders with zero broken nodes against a bare
   `SpecContext` — no `bind`, no `on` — so a sample never depends on the hosting demo's context
   (pinned by tests/samples.test.js). The first entry is also what 'New form' opens (F3). */
import {SPEC_SCHEMA} from '../../spec/registry.js';

export const SAMPLES: {name: string, spec: object}[] = [
  {name: 'Blank', spec: {$schema: SPEC_SCHEMA, root: {tag: 'u2-form', name: 'form1'}}},
  {name: 'Data entry', spec: {$schema: SPEC_SCHEMA, root: {
    tag: 'u2-panel', name: 'dataEntry', children: [
      {tag: 'u2-form', name: 'entryForm', children: [
        {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}},
        {tag: 'u2-number-input', name: 'amountInput', props: {label: 'Amount', value: 1, min: 0, max: 100}},
        {tag: 'u2-choice-input', name: 'categoryInput',
          props: {label: 'Category', items: ['Reagent', 'Solvent', 'Catalyst']}},
        {tag: 'u2-bool-input', name: 'activeInput', props: {label: 'Active', value: true}},
        {tag: 'u2-text-area', name: 'notesInput', props: {label: 'Notes'}},
      ]},
      {tag: 'u2-div-h', name: 'entryButtons', children: [
        {tag: 'u2-button', name: 'saveButton', props: {text: 'Save', primary: true}},
        {tag: 'u2-button', name: 'cancelButton', props: {text: 'Cancel'}},
      ]},
    ]}}},
  {name: 'Master–detail', spec: {$schema: SPEC_SCHEMA, root: {
    tag: 'u2-splitter', name: 'masterDetail', props: {sizes: [0.3, 0.7]}, children: [
      {tag: 'u2-panel', name: 'masterPane', children: [
        {tag: 'u2-text-input', name: 'searchInput', props: {label: 'Search', search: true, inline: true}},
        {tag: 'u2-multi-choice-input', name: 'itemList',
          props: {label: 'Items', items: ['Aspirin', 'Ibuprofen', 'Paracetamol', 'Naproxen']}},
      ]},
      {tag: 'u2-panel', name: 'detailPane', children: [
        {tag: 'h2', name: 'detailTitle', props: {text: 'Details'}},
        {tag: 'u2-form', name: 'detailForm', children: [
          {tag: 'u2-text-input', name: 'detailName', props: {label: 'Name', value: 'Aspirin'}},
          {tag: 'u2-text-area', name: 'detailNotes', props: {label: 'Notes'}},
        ]},
      ]},
    ]}}},
  {name: 'Settings', spec: {$schema: SPEC_SCHEMA, root: {
    tag: 'u2-tabs', name: 'settingsTabs', children: [
      {tag: 'u2-panel', name: 'generalPane', props: {title: 'General'}, children: [
        {tag: 'u2-form', name: 'generalForm', children: [
          {tag: 'u2-text-input', name: 'titleInput', props: {label: 'Title', value: 'My dashboard'}},
          {tag: 'u2-bool-input', name: 'autoSaveInput', props: {label: 'Auto-save', value: true}},
        ]},
      ]},
      {tag: 'u2-panel', name: 'appearancePane', props: {title: 'Appearance'}, children: [
        {tag: 'u2-form', name: 'appearanceForm', children: [
          {tag: 'u2-color-input', name: 'accentInput', props: {label: 'Accent', value: '#1f77b4'}},
          {tag: 'u2-slider-input', name: 'opacityInput', props: {label: 'Opacity', value: 80, min: 0, max: 100}},
        ]},
      ]},
      {tag: 'u2-panel', name: 'advancedPane', props: {title: 'Advanced'}, children: [
        {tag: 'u2-form', name: 'advancedForm', children: [
          {tag: 'u2-number-input', name: 'cacheInput', props: {label: 'Cache size', value: 100, min: 0, max: 1000}},
          {tag: 'u2-choice-input', name: 'logLevelInput',
            props: {label: 'Log level', items: ['error', 'warning', 'info', 'debug'], value: 'info'}},
        ]},
      ]},
    ]}}},
  {name: 'Dashboard', spec: {$schema: SPEC_SCHEMA, root: {
    tag: 'u2-panel', name: 'dashboard', children: [
      {tag: 'u2-breadcrumbs', name: 'crumbs', props: {items: ['Home', 'Projects', 'Demo']}},
      {tag: 'h2', name: 'dashboardTitle', props: {text: 'Overview'}},
      {tag: 'u2-div-h', name: 'dashboardRow', children: [
        {tag: 'u2-panel', name: 'settingsCard', children: [
          {tag: 'u2-property-grid', name: 'settingsGrid', props: {
            properties: [{name: 'Title', type: 'string'}, {name: 'Rows', type: 'int', min: 0}],
            values: {Title: 'Demo', Rows: 20},
          }},
        ]},
        {tag: 'u2-panel', name: 'logoCard', children: [
          {tag: 'u2-image-input', name: 'logoInput', props: {label: 'Logo'}},
        ]},
        {tag: 'u2-panel', name: 'summaryCard', children: [
          {tag: 'u2-form', name: 'summaryForm', children: [
            {tag: 'u2-text-input', name: 'ownerInput', props: {label: 'Owner', value: 'admin'}},
            {tag: 'u2-number-input', name: 'countInput', props: {label: 'Items', value: 3}},
          ]},
        ]},
      ]},
    ]}}},
];
