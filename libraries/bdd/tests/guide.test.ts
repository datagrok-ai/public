import assert from 'node:assert/strict';
import {spawnSync} from 'node:child_process';
import {existsSync, mkdtempSync, readFileSync, rmSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {test} from 'node:test';
import {captionOf, parseTitle, slugOf} from '../src/runtime/guide.js';

test('a step title splits into its keyword and text', () => {
  assert.deepEqual(parseTitle('When user clicks on browse tab'), {keyword: 'When', text: 'user clicks on browse tab'});
  assert.deepEqual(parseTitle('And the table should have 5 rows'),
    {keyword: 'And', text: 'the table should have 5 rows'});
  assert.deepEqual(parseTitle('something without a keyword'), {keyword: '', text: 'something without a keyword'});
});

test('a data table is told inline after the step, each two-cell row as name = value', () => {
  assert.equal(captionOf('user adds a scatter plot viewer with:', 'action', [['X', 'AGE'], ['Y', 'HEIGHT']]),
    'Add a scatter plot viewer with X = AGE, Y = HEIGHT');
  assert.equal(captionOf('user fills in:', 'action', [['Name', 'x'], ['Type', 'y']]), 'Fill in Name = x, Type = y');
});

test('an action reads as an instruction: the verb imperative, "inside" as "in", a tree path with separators', () => {
  assert.equal(captionOf('user clicks on "Open local file" icon inside browse toolbar', 'action'),
    'Click on "Open local file" icon in browse toolbar');
  assert.equal(captionOf('user uploads "x.csv" through "Open local file" icon inside browse toolbar', 'action'),
    'Upload "x.csv" through "Open local file" icon in browse toolbar');
  assert.equal(captionOf('user presses Escape', 'action'), 'Press Escape');
  assert.equal(captionOf('user double-clicks on Files---Demo tree node', 'action'),
    'Double-click on Files › Demo tree node');
  assert.equal(captionOf('user hovers over the "stats" area of box plot viewer', 'action'),
    'Hover over the "stats" area of box plot viewer');
  assert.equal(captionOf('user applies the layout', 'action'), 'Apply the layout');
  assert.equal(captionOf('user pushes the button', 'action'), 'Push the button');
  assert.equal(captionOf('user copies the cell', 'action'), 'Copy the cell');
  assert.equal(captionOf('user is logged in', 'setup'), 'User is logged in');
  assert.equal(captionOf('the browse panel is open', 'action'), 'The browse panel is open');
});

test('a check reads as the fact it verifies', () => {
  assert.equal(captionOf('the table should have 5 rows', 'check'), 'The table has 5 rows');
  assert.equal(captionOf('the "browse-import" view should be current', 'check'), 'The "browse-import" view is current');
  assert.equal(captionOf('the browse tree should not be visible', 'check'), 'The browse tree is not visible');
  assert.equal(captionOf('the top menu should list:', 'check'), 'The top menu lists:');
  assert.equal(captionOf('5 rows should be selected', 'check'), '5 rows are selected');
  assert.equal(captionOf('9 rows of table "SPGI-linked1" should pass the filter', 'check'),
    '9 rows of table "SPGI-linked1" pass the filter');
  assert.equal(captionOf('grid should show 9 rows', 'check'), 'Grid shows 9 rows');
  assert.equal(captionOf('the filter should pass exactly the rows where "SEX" is "F"', 'check'),
    'The filter passes exactly the rows where "SEX" is "F"');
  assert.equal(captionOf('the tooltip should not show columns "AGE"', 'check'), 'The tooltip does not show columns "AGE"');
  assert.equal(captionOf('no rows should be selected', 'check'), 'No rows are selected');
  assert.equal(captionOf('the status should be "ready"', 'check'), 'The status is "ready"');
  assert.equal(captionOf('rows 1 to 5 should be selected', 'check'), 'Rows 1 to 5 are selected');
  assert.equal(captionOf('the value of "AGE" column in row 3 should be "5"', 'check'),
    'The value of "AGE" column in row 3 is "5"');
});

test('a scenario name becomes a directory name', () => {
  assert.equal(slugOf('Open local file imports a CSV into a table view'),
    'open-local-file-imports-a-csv-into-a-table-view');
  assert.equal(slugOf('  ---  '), 'untitled');
});

test('the renderer turns a synthetic manifest into a video, a GIF, a thumb and steps.md', {timeout: 120000}, (t) => {
  const python = process.env.BDD_PYTHON ?? (process.platform === 'win32' ? 'py' : 'python3');
  const dir = mkdtempSync(join(tmpdir(), 'bdd-guide-'));
  try {
    const result = spawnSync(python, [join('tool', 'guide-render.py'), '--selftest', dir, '--quiet'],
      {encoding: 'utf8'});
    if (result.error || result.status === 3) {
      t.skip(`no python/Pillow/ffmpeg here: ${result.error?.message ?? result.stderr}`);
      return;
    }
    assert.equal(result.status, 0, result.stderr);
    assert.ok(existsSync(join(dir, 'guide.mp4')) || existsSync(join(dir, 'guide.webm')), 'a video');
    assert.ok(existsSync(join(dir, 'guide.gif')), 'a gif');
    assert.ok(existsSync(join(dir, 'guide-thumb.png')), 'a thumb');
    assert.ok(existsSync(join(dir, 'step-01.png')) && existsSync(join(dir, 'step-02.png')), 'a picture per step');
    const md = readFileSync(join(dir, 'steps.md'), 'utf8');
    assert.match(md, /^# Two synthetic steps/);
    assert.match(md, /1\. Click on the blue button/);
    assert.match(md, /2\. ✓ The page is blue/);
  }
  finally {
    rmSync(dir, {recursive: true, force: true});
  }
});
