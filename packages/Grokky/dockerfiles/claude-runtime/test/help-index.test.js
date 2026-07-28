// buildHelpIndex — generates workspace/help/INDEX.md, the single Read that replaced an
// unbounded grep chain for help grounding (see docs/LATENCY.md tier 1a). What it indexes,
// and what it deliberately skips, directly determines help-turn round-trip count.

const test = require('node:test');
const assert = require('node:assert');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {buildHelpIndex} = require('../dist/help-index.js');

// Builds a throwaway workspace/help tree from {relPath: contents} and returns the INDEX.md text.
function indexOf(files) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'helpidx-'));
  for (const [rel, body] of Object.entries(files)) {
    const full = path.join(root, 'help', rel);
    fs.mkdirSync(path.dirname(full), {recursive: true});
    fs.writeFileSync(full, body);
  }
  buildHelpIndex(root);
  const out = fs.readFileSync(path.join(root, 'help', 'INDEX.md'), 'utf8');
  fs.rmSync(root, {recursive: true, force: true});
  return out;
}

test('frontmatter title, description and keywords all reach the index line', () => {
  const idx = indexOf({
    'visualize/scatter-plot.md': `---
title: Scatter plot
description: Plots two numeric columns against each other
keywords:
  - scatterplot
  - xy chart
---
# Scatter plot
body`,
  });
  assert.match(idx, /visualize\/scatter-plot\.md — Scatter plot/);
  assert.match(idx, /Plots two numeric columns against each other/);
  assert.match(idx, /\[scatterplot, xy chart\]/);
});

test('a page with no frontmatter falls back to its first H1', () => {
  const idx = indexOf({'basics/intro.md': '# Getting started\n\nsome text'});
  assert.match(idx, /basics\/intro\.md — Getting started/);
});

test('a page with neither frontmatter nor H1 falls back to its filename', () => {
  const idx = indexOf({'basics/orphan.md': 'no heading at all'});
  assert.match(idx, /basics\/orphan\.md — orphan/);
});

test('quoted frontmatter values are unquoted', () => {
  const idx = indexOf({'a/b.md': '---\ntitle: "Quoted Title"\ndescription: \'single\'\n---\nx'});
  assert.match(idx, /— Quoted Title — single/);
  assert.ok(!idx.includes('"Quoted Title"'));
});

test('QA scenario pages (-test.md) are excluded so answers cite real docs', () => {
  const idx = indexOf({
    'viewers/scatter.md': '# Scatter',
    'viewers/scatter-test.md': '# Scatter test checklist',
  });
  assert.match(idx, /viewers\/scatter\.md/);
  assert.ok(!idx.includes('scatter-test.md'), '-test.md pages must not be indexed');
});

test('underscore-prefixed files, INDEX.md and CLAUDE.md are excluded', () => {
  const idx = indexOf({
    'real.md': '# Real',
    '_draft.md': '# Draft',
    'CLAUDE.md': '# Agent instructions',
  });
  assert.match(idx, /real\.md/);
  assert.ok(!idx.includes('_draft.md'));
  assert.ok(!idx.includes('CLAUDE.md'));
});

test('internal and asset directories are skipped', () => {
  const idx = indexOf({
    'keep/page.md': '# Keep',
    '_internal/secret.md': '# Secret',
    'uploads/blob.md': '# Blob',
    'img/pic.md': '# Pic',
  });
  assert.match(idx, /keep\/page\.md/);
  for (const skipped of ['secret.md', 'blob.md', 'pic.md'])
    assert.ok(!idx.includes(skipped), `${skipped} should be skipped`);
});

test('pages are grouped under a section header per top-level directory', () => {
  const idx = indexOf({
    'access/login.md': '# Login',
    'visualize/chart.md': '# Chart',
  });
  assert.match(idx, /^## access$/m);
  assert.match(idx, /^## visualize$/m);
});

test('the header tells the model to stop searching when nothing matches', () => {
  const idx = indexOf({'a.md': '# A'});
  assert.match(idx, /If no page covers the topic, the docs do not cover it/);
});

test('mdx pages are indexed alongside md', () => {
  const idx = indexOf({'guide/tour.mdx': '# Tour'});
  assert.match(idx, /guide\/tour\.mdx — Tour/);
});

test('a missing help directory is a no-op rather than a crash', () => {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'helpidx-empty-'));
  assert.doesNotThrow(() => buildHelpIndex(root));
  assert.ok(!fs.existsSync(path.join(root, 'help', 'INDEX.md')));
  fs.rmSync(root, {recursive: true, force: true});
});
