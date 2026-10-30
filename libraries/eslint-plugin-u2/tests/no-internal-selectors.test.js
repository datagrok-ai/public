'use strict';

const tester = require('./tester');
const rule = require('../rules/no-internal-selectors');

const internal = (code, count = 1) => ({code, errors: Array(count).fill({messageId: 'internal'})});

tester.run('no-internal-selectors', rule, {
  valid: [
    `root.querySelector('.hit-triage-compute-descriptors');`,
    `el.classList.add('hit-triage-row', 'hit-triage-row-current');`,
    `el.style.setProperty('--u2-row-height', '28px');`,
    `if (el.matches(':hover')) open();`,
    `const cls = 'u2-tree'; log(cls);`,
    `document.querySelectorAll(selector);`,
    `list.root.classList.add('u2-adaptive'); meta.classList.add('u2-p2', 'u2-p1');`,
    `root.querySelectorAll('.u2-adaptive .u2-p2');`,
  ],
  invalid: [
    internal(`await waitUntil(() =>
      root.querySelector('.hit-triage-compute-descriptors .u2-tree:not([hidden])') !== null);`),
    internal(`const label = (r) => r.querySelector('.u2-tree-label')?.textContent ?? '';`),
    internal(`root.querySelector('[data-u2-name="template-picker"] select, [data-u2-name="template-picker"] ' +
      'input');`),
    internal(`if (e.target.closest('[data-u2="card"]') === null) select(null);`),
    internal('const rows = host.querySelectorAll(`.u2-catalog-row[data-key="${key}"]`);'),
    internal(`cell.classList.add('u2-input-fill');`),
    internal(`row.classList.toggle('hit-triage-row', on); row.classList.toggle('u2-list-row-current', on);`),
    internal(`const opts = form.getElementsByClassName('u2-input-options');`),
    internal(`root.querySelector('.u2-adaptive .u2-tree-label');`),
    internal(`el.classList.add('u2-p1-label');`),
  ],
});
