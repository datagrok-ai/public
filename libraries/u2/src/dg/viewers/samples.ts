/* The Open menu's viewer sample (VP-15/VP-16): platform-only, so it lives apart from the core
   `SAMPLES` the bare-context gallery test pins. Its `orders` is the `sample` policy over three
   rows — real DataFrame rows in a real grid at design time, without running anything. */
import {SPEC_SCHEMA} from '../../spec/registry.js';
import type {Registry} from '../../spec/registry.js';

export const VIEWER_SAMPLES: {name: string, spec: object}[] = [
  {name: 'Master–detail (grid + filters)', spec: {$schema: SPEC_SCHEMA,
    components: [{tag: 'u2-func-source', name: 'orders',
      props: {func: 'demoOrders', params: {days: 30}, designData: 'sample', debounce: 50,
        sample: [{orderId: 1001, customer: 'Aspirin Labs', city: 'Kyiv', total: 1240},
          {orderId: 1002, customer: 'Bayer', city: 'Lviv', total: 380},
          {orderId: 1003, customer: 'Roche', city: 'Basel', total: 2150}]}}],
    root: {tag: 'u2-splitter', name: 'masterDetail', props: {sizes: [0.25, 0.45, 0.3]}, children: [
      {tag: 'u2-viewer-filters', name: 'filters', bind: {table: '$.orders'},
        props: {columnNames: ['city', 'total']}},
      {tag: 'u2-viewer-grid', name: 'grid', bind: {table: '$.orders'}, props: {allowEdit: false}},
      {tag: 'u2-panel', name: 'detailPane', children: [
        {tag: 'u2-viewer-form', name: 'form', bind: {table: '$.orders'}},
        {tag: 'u2-form', name: 'detailForm', children: [
          {tag: 'u2-text-input', name: 'customerField', props: {label: 'Customer'},
            bind: {value: '$.orders.currentRow.customer'}},
          {tag: 'u2-number-input', name: 'rowIdx', props: {label: 'Row', mode: 'int'},
            bind: {value: '$.orders.currentRowIdx'}}]}]}]}}},
];

/** The viewer samples, where the registry has the viewer tags — none on the core registry. */
export function platformSamples(reg: Registry | undefined): {name: string, spec: object}[] {
  return reg?.get('u2-viewer-grid') !== undefined ? VIEWER_SAMPLES : [];
}
