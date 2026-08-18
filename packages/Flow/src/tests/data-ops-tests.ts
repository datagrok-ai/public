/** Flow's own data operations: the functions, the nodes they produce,
 *  and the editors that make their parameters fillable. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {effectiveFuncInputs, customEditorFor, funcValidatorOf} from '../utils/func-input-overrides';
import {nodeMissingRequirements} from '../rete/scheme';
import {categorizeFunc} from '../panel/function-browser';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';
import {
  AGGREGATION_TYPES, aggregationProblems, parseAggregations, randomIndices,
  aggregate, deleteColumns, deleteRows, expressionToColumn, extractRows, filterRandomRows,
  filterRows, selectRows, tagColumns, unpivot,
} from '../ops/data-ops';
import {
  aggregationEditor, columnsForAggregation, describeAggregation, pruneAggregations,
  serializeAggregations,
} from '../panel/editors/aggregation-editor';
import {rowConditionEditor} from '../panel/editors/expression-editor';

function typeNameOf(nqName: string): string | null {
  return getRegisteredFuncs().find((f) => {
    try {
      return f.func.nqName === nqName;
    } catch {
      return false;
    }
  })?.nodeTypeName ?? null;
}

function sampleTable(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('team', ['a', 'a', 'b', 'b', 'c', 'c']),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'score', [1, 5, 2, 8, 3, 9]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'weight', [10, 20, 30, 40, 50, 60]),
  ]);
}

category('Flow: data ops', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('Filter Rows keeps exactly the matching rows', async () => {
    const table = sampleTable();
    const result = await filterRows(table, '${score} > 3');
    expect(result.rowCount, 3, 'rows over 3');
    expect(result.col('score')!.toList().join(','), '5,8,9');
    expect(table.rowCount, 6, 'the input table is untouched');
  });

  test('Delete Rows is the exact complement of Filter Rows', async () => {
    const table = sampleTable();
    const kept = await deleteRows(table, '${score} > 3');
    expect(kept.rowCount, 3);
    expect(kept.col('score')!.toList().join(','), '1,2,3');
  });

  test('Extract Rows narrows columns as well as rows', async () => {
    const table = sampleTable();
    const result = await extractRows(table, '${team} == "b"', [table.col('team')!, table.col('score')!]);
    expect(result.rowCount, 2);
    expect(result.columns.names().join(','), 'team,score', 'only the chosen columns');
    const all = await extractRows(table, '${team} == "b"', null);
    expect(all.columns.length, 3, 'a blank column list keeps them all');
  });

  test('Select Rows selects on the table and passes it on', async () => {
    const table = sampleTable();
    const result = await selectRows(table, '${score} > 3', true);
    expect(result === table, true, 'the same table flows on — selection is table state');
    expect(table.selection.trueCount, 3);
    await selectRows(table, '${team} == "a"', true);
    expect(table.selection.trueCount, 2, 'the previous selection was dropped');
    await selectRows(table, '${score} > 8', false);
    expect(table.selection.trueCount, 3, 'and adds to it when asked to');
  });

  test('a non-boolean condition is refused with a readable message', async () => {
    let message = '';
    try {
      await filterRows(sampleTable(), '${score} + 1');
    } catch (e) {
      message = String((e as Error)?.message ?? e);
    }
    expect(message.length > 0, true, 'it throws rather than filtering by a number');
    expect(message.includes('true/false'), true, `names the problem: ${message}`);
  });

  test('an empty condition is refused', async () => {
    let threw = false;
    try {
      await filterRows(sampleTable(), '   ');
    } catch {
      threw = true;
    }
    expect(threw, true);
  });

  test('random sampling is reproducible and drawn without replacement', async () => {
    const a = randomIndices(100, 10, 42);
    const b = randomIndices(100, 10, 42);
    expect(a.join(','), b.join(','), 'the same seed draws the same rows');
    expect(new Set(a).size, 10, 'no duplicates');
    expect(a.every((i) => i >= 0 && i < 100), true, 'all in range');
    expect(randomIndices(100, 10, 7).join(',') !== a.join(','), true, 'a different seed differs');
    expect(randomIndices(5, 50, 1).length, 5);

    const table = sampleTable();
    const sample = filterRandomRows(table, 2, 42);
    expect(sample.rowCount, 2);
    expect(table.rowCount, 6, 'the input table is untouched');
  });

  test('Delete Columns drops the chosen ones and refuses to empty the table', async () => {
    const table = sampleTable();
    const result = deleteColumns(table, [table.col('weight')!]);
    expect(result.columns.names().join(','), 'team,score');
    expect(table.columns.length, 3, 'the input table is untouched');

    let threw = false;
    try {
      deleteColumns(table, table.columns.toList());
    } catch {
      threw = true;
    }
    expect(threw, true, 'deleting every column is refused');
  });

  test('Tag Columns writes the tag and passes the table on', async () => {
    const table = sampleTable();
    const result = tagColumns(table, [table.col('score')!], 'units', 'points');
    expect(result === table, true);
    expect(table.col('score')!.getTag('units'), 'points');
    let threw = false;
    try {
      tagColumns(table, [table.col('score')!], '  ', 'x');
    } catch {
      threw = true;
    }
    expect(threw, true, 'a blank tag name is refused');
  });

  test('Expression To Column adds the computed column', async () => {
    const table = sampleTable();
    const column = await expressionToColumn(table, '${score} * 2', 'doubled', 'int');
    expect(column.name, 'doubled');
    expect(table.col('doubled')!.toList().join(','), '2,10,4,16,6,18');
  });

  test('Aggregate groups and aggregates', async () => {
    const table = sampleTable();
    const result = aggregate(table, [table.col('team')!],
      JSON.stringify([{column: 'score', type: 'sum', name: 'total'}]), null);
    expect(result.rowCount, 3, 'one row per team');
    expect(result.col('total')!.toList().join(','), '6,10,12');
  });

  test('Aggregate pivots when given a pivot column', async () => {
    const table = sampleTable();
    const result = aggregate(table, [], JSON.stringify([{column: 'score', type: 'sum'}]),
      [table.col('team')!]);
    expect(result.rowCount, 1, 'no group-by folds everything into one row');
    expect(result.columns.length > 1, true, 'one result column per pivot value');
  });

  test('Aggregate refuses an incomplete list instead of aggregating nothing', async () => {
    const table = sampleTable();
    for (const bad of ['', '[]', JSON.stringify([{column: '', type: 'avg'}])]) {
      let threw = false;
      try {
        aggregate(table, [table.col('team')!], bad, null);
      } catch {
        threw = true;
      }
      expect(threw, true, `refused: ${bad || '(blank)'}`);
    }
  });

  test('every offered aggregation actually runs', async () => {
    // incl. `concat unique`, a STR_AGG the TS signature does not admit
    const table = sampleTable();
    for (const type of AGGREGATION_TYPES) {
      const column = type === 'count' ? '' : (type === 'concat unique' ? 'team' : 'score');
      const result = aggregate(table, [table.col('team')!],
        JSON.stringify([{column, type, name: 'agg'}]), null);
      expect(result.rowCount, 3, `${type} produced one row per group`);
      expect(result.col('agg') !== null, true, `${type} produced its column`);
    }
  });

  test('Unpivot folds the merged columns into category and value rows', async () => {
    const table = sampleTable();
    const result = await unpivot(table, [table.col('team')!],
      [table.col('score')!, table.col('weight')!], 'metric', 'amount');
    expect(result.rowCount, 12, 'six rows times two merged columns');
    expect(result.col('metric') !== null, true, 'the category column is named as asked');
    expect(result.col('amount') !== null, true, 'and so is the value column');
    let threw = false;
    try {
      await unpivot(table, null, [], 'metric', 'amount');
    } catch {
      threw = true;
    }
    expect(threw, true, 'no merge columns is refused before the call');
  });
});

category('Flow: data op nodes', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  // each declares meta.includeInFlow: true — no allowlist entry needed
  const OPS = ['filterRows', 'deleteRows', 'extractRows', 'selectRows', 'filterRandomRows',
    'selectRandomRows', 'deleteColumns', 'tagColumns', 'expressionToColumn', 'aggregate', 'unpivot'];

  test('all of them are registered as nodes', async () => {
    for (const name of OPS)
      expect(typeNameOf(`Flow:${name}`) !== null, true, `Flow:${name} is a node type`);
  });

  test('the FuncCall-typed originals they replace stay out of the catalog', async () => {
    for (const gone of ['core:FilterRows', 'core:DeleteRows', 'core:SelectRows', 'core:ExtractRows',
      'core:DeleteColumns', 'core:TagColumns', 'core:Subset', 'core:Aggregate'])
      expect(typeNameOf(gone), null, `${gone} is not a node`);
  });

  test('every condition node takes a table plus a plain string condition', async () => {
    for (const name of ['filterRows', 'deleteRows', 'extractRows', 'selectRows']) {
      const typeName = typeNameOf(`Flow:${name}`);
      if (!typeName) continue;
      const e = makeEditor();
      try {
        const node = await addNode(e.flow, typeName);
        const inputs = effectiveFuncInputs(node.dgFunc!);
        expect(String(inputs[0].propertyType), 'dataframe', `${name} leads with a table`);
        const condition = inputs.find((p) => p.name === 'condition');
        expect(String(condition?.propertyType), 'string', `${name}: condition is a plain string`);
        expect(node.requiredInputs.includes('table'), true, `${name}: table required`);
        expect(node.requiredInputs.includes('condition'), true, `${name}: condition required`);
        expect(nodeMissingRequirements(node, () => false).length > 0, true,
          `${name}: a fresh node is not runnable`);
      } finally {
        destroyEditor(e);
      }
    }
  });

  test('column parameters are real column_list slots bound to the table', async () => {
    const cases: Array<[string, string[]]> = [
      ['deleteColumns', ['columns']],
      ['tagColumns', ['columns']],
      ['aggregate', ['groupByColumns', 'pivotColumns']],
      ['unpivot', ['copyColumns', 'mergeColumns']],
      ['extractRows', ['columns']],
    ];
    for (const [name, params] of cases) {
      const typeName = typeNameOf(`Flow:${name}`);
      if (!typeName) continue;
      const e = makeEditor();
      try {
        const node = await addNode(e.flow, typeName);
        const inputs = effectiveFuncInputs(node.dgFunc!);
        const tables = node.properties['columnTables'] as Record<string, string> | undefined;
        for (const param of params) {
          const prop = inputs.find((p) => p.name === param);
          expect(String(prop?.propertyType), 'column_list', `${name}.${param} is a column list`);
          expect(tables?.[param], 'table', `${name}.${param} resolves against the table`);
        }
      } finally {
        destroyEditor(e);
      }
    }
  });

  test('the mandatory column lists are required, the optional ones are not', async () => {
    const required: Array<[string, string]> = [
      ['deleteColumns', 'columns'], ['tagColumns', 'columns'], ['unpivot', 'mergeColumns'],
    ];
    const optional: Array<[string, string]> = [
      ['aggregate', 'groupByColumns'], ['aggregate', 'pivotColumns'],
      ['unpivot', 'copyColumns'], ['extractRows', 'columns'],
    ];
    for (const [name, param] of [...required, ...optional]) {
      const typeName = typeNameOf(`Flow:${name}`);
      if (!typeName) continue;
      const e = makeEditor();
      try {
        const node = await addNode(e.flow, typeName);
        const shouldBeRequired = required.some(([n, p]) => n === name && p === param);
        expect(node.requiredInputs.includes(param), shouldBeRequired,
          `${name}.${param} required = ${shouldBeRequired}`);
      } finally {
        destroyEditor(e);
      }
    }
  });

  test('they land in the toolbox categories their signatures imply', async () => {
    const expected: Record<string, string> = {
      filterRows: 'Transform Tables',
      deleteRows: 'Transform Tables',
      extractRows: 'Transform Tables',
      selectRows: 'Transform Tables',
      filterRandomRows: 'Transform Tables',
      deleteColumns: 'Transform Tables',
      tagColumns: 'Transform Tables',
      aggregate: 'Transform Tables',
      unpivot: 'Transform Tables',
      expressionToColumn: 'Column Operations',
    };
    for (const [name, category] of Object.entries(expected)) {
      const info = getRegisteredFuncs().find((f) => {
        try {
          return f.func.nqName === `Flow:${name}`;
        } catch {
          return false;
        }
      });
      if (!info) continue;
      expect(categorizeFunc(info.func, info.role ?? null, info.packageName), category,
        `Flow:${name} is filed under ${category}`);
    }
  });

  test('nodes are titled in words, not in camelCase', async () => {
    // `//friendlyName:` does not survive publishing, so the header humanizes the raw name
    const expected: Record<string, string> = {
      filterRows: 'Filter Rows',
      filterRandomRows: 'Filter Random Rows',
      deleteColumns: 'Delete Columns',
      expressionToColumn: 'Expression To Column',
      aggregate: 'Aggregate',
      unpivot: 'Unpivot',
    };
    const e = makeEditor();
    try {
      for (const [name, title] of Object.entries(expected)) {
        const typeName = typeNameOf(`Flow:${name}`);
        if (!typeName) continue;
        const node = await addNode(e.flow, typeName);
        expect(node.label, title, `Flow:${name} is titled "${title}"`);
      }
    } finally {
      destroyEditor(e);
    }
  });

  test('choice parameters are literal lists, not free text', async () => {
    const typeName = typeNameOf('Flow:expressionToColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const type = effectiveFuncInputs(node.dgFunc!).find((p) => p.name === 'type')!;
      expect((type.choices ?? []).includes('auto'), true, 'the type parameter offers its choices');
      expect((type.choices ?? []).includes('bool'), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('an Aggregate node with a half-filled list is not runnable', async () => {
    const typeName = typeNameOf('Flow:aggregate');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const validator = funcValidatorOf(node.dgFunc!);
      expect(validator !== undefined, true, 'aggregate has a readiness check');

      // blank is caught by the generic check — the validator exists for a non-blank list that says nothing
      node.inputValues['aggregations'] = JSON.stringify([{column: '', type: 'avg'}]);
      expect(validator!(node).length > 0, true, 'an aggregation over no column blocks the node');
      expect(validator!(node)[0].includes('column'), true, 'and names what is missing');

      node.inputValues['aggregations'] = JSON.stringify([{column: 'score', type: 'avg'}]);
      expect(validator!(node).length, 0, 'a complete list is ready');

      node.inputValues['aggregations'] = JSON.stringify([{column: '', type: 'count'}]);
      expect(validator!(node).length, 0, 'count needs no column');
    } finally {
      destroyEditor(e);
    }
  });

  test('a formula the editor rejects makes the node unready', async () => {
    const typeName = typeNameOf('Flow:filterRows');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const validator = funcValidatorOf(node.dgFunc!);
      expect(validator !== undefined, true, 'a condition has a readiness check');

      node.inputValues['condition'] = '${score} + 1';
      expect(validator!(node).length, 0,
        'nothing has checked it yet — an unverified formula must not block a flow');

      const ed = rowConditionEditor(
        DG.Property.fromOptions({name: 'condition', type: 'string'} as never), {
          inputValue: () => undefined, columns: () => null, table: () => null,
          isConnected: () => true, watch: () => {}, node,
        } as never);
      // validity is decided by evaluation — the hosted editor announces verdicts as DOM events
      const validate = (expression: string, error: string): void => {
        ed.element.dispatchEvent(new CustomEvent('expression-validated', {detail: {expression, error}}));
      };

      validate('${score} + 1', 'The formula is int, not a true/false condition. Compare something');
      expect(ed.element.hasAttribute('data-invalid'), true, 'the editor marks itself');
      const problems = validator!(node);
      expect(problems.length, 1, 'and the node refuses to run');
      expect(problems[0].includes('true/false'), true, 'saying why');

      node.inputValues['condition'] = '${score} > 1';
      expect(validator!(node).length, 0, 'a verdict for other text does not apply');

      validate('${score} > 1', '');
      expect(ed.element.hasAttribute('data-invalid'), false, 'a good formula clears the mark');
      expect(validator!(node).length, 0, 'and the node is runnable again');
      ed.detach!();
    } finally {
      destroyEditor(e);
    }
  });

  test('every expression parameter gets the formula editor', async () => {
    const e = makeEditor();
    try {
      for (const [name, param] of [['filterRows', 'condition'], ['deleteRows', 'condition'],
        ['extractRows', 'condition'], ['selectRows', 'condition'],
        ['expressionToColumn', 'expression'], ['aggregate', 'aggregations']] as Array<[string, string]>) {
        const typeName = typeNameOf(`Flow:${name}`);
        if (!typeName) continue;
        const node = await addNode(e.flow, typeName);
        expect(customEditorFor(node.dgFunc!, param) !== null, true,
          `Flow:${name}.${param} has a custom editor`);
      }
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: data op editors', () => {
  /** Stand-in for the panel context, so editor states are testable without a live node. */
  function fakeCtx(columns: DG.Column[] | null, connected: boolean) {
    const watchers: Array<(v: unknown) => void> = [];
    return {
      ctx: {
        inputValue: () => undefined,
        columns: () => columns,
        table: () => null,
        isConnected: () => connected,
        watch: (_name: string, cb: (v: unknown) => void) => {watchers.push(cb);},
      },
      fire: (): void => watchers.forEach((cb) => cb(null)),
    };
  }

  const stringParam = (name: string): DG.Property =>
    DG.Property.fromOptions({name, type: 'string'} as never);

  test('parseAggregations tolerates anything the panel can hold mid-edit', async () => {
    expect(parseAggregations('').length, 0, 'blank');
    expect(parseAggregations(null).length, 0, 'null');
    expect(parseAggregations('{not json').length, 0, 'malformed');
    expect(parseAggregations('{"a": 1}').length, 0, 'an object, not a list');
    expect(parseAggregations(JSON.stringify([{column: 'x', type: 'nonsense'}])).length, 0,
      'an unknown aggregation is dropped, not thrown on');
    const parsed = parseAggregations(JSON.stringify([{column: 'x', type: 'avg', name: 'mean x'}]));
    expect(parsed.length, 1);
    expect(parsed[0].name, 'mean x');
  });

  test('aggregationProblems names exactly what is missing', async () => {
    expect(aggregationProblems([])[0].includes('at least one'), true);
    expect(aggregationProblems([{column: '', type: 'avg'}]).length, 1);
    expect(aggregationProblems([{column: '', type: 'count'}]).length, 0, 'count needs no column');
    expect(aggregationProblems([{column: 'x', type: 'avg'}]).length, 0);
    expect(aggregationProblems([{column: '', type: 'avg'}, {column: '', type: 'avg'}]).length, 1);
  });

  test('an aggregation only offers columns it can consume', async () => {
    const columns = [
      DG.Column.fromStrings('name', ['a']),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'score', [1]),
    ];
    expect(columnsForAggregation('avg', columns).join(','), 'score', 'numeric aggregations, numbers only');
    expect(columnsForAggregation('unique', columns).join(','), 'name,score', 'counting works on anything');
    expect(columnsForAggregation('count', columns).length, 0, 'count takes no column at all');
    expect(columnsForAggregation('avg', null).length, 0, 'no table, no choices');
  });

  test('aggregations over columns that are gone are pruned', async () => {
    const specs = [{column: 'score', type: 'avg'}, {column: 'ghost', type: 'sum'},
      {column: '', type: 'count'}, {column: '', type: 'avg'}];
    const kept = pruneAggregations(specs, ['score']);
    expect(kept.length, 3, 'the aggregation over a missing column is dropped');
    expect(kept.some((a) => a.column === 'ghost'), false, 'and it is the stale one');
    expect(kept.some((a) => a.type === 'count'), true, 'column-less ones always survive');
    expect(kept.filter((a) => a.column === '' && a.type === 'avg').length, 1,
      'an unfinished row survives');
  });

  test('serializeAggregations stores nothing when there is nothing to store', async () => {
    expect(serializeAggregations([]), '', 'empty stays blank, not "[]"');
    expect(serializeAggregations([{column: 'x', type: 'avg', name: '  '}]),
      JSON.stringify([{column: 'x', type: 'avg'}]), 'a blank result name is dropped');
    expect(serializeAggregations([{column: 'x', type: 'count'}]),
      JSON.stringify([{column: '', type: 'count'}]), 'count never carries a column');
    expect(describeAggregation({column: 'x', type: 'avg', name: 'm'}), 'avg(x) → m');
    expect(describeAggregation({column: '', type: 'count'}), 'count');
  });

  test('the aggregation editor blocks rather than showing combos it cannot fill', async () => {
    const {ctx} = fakeCtx(null, false);
    const ed = aggregationEditor(stringParam('aggregations'), ctx as never);
    ed.setValue('');
    expect(ed.element.getAttribute('data-blocked'), 'true', 'no table means a notice');
    expect(ed.element.querySelectorAll('select').length, 0, 'and zero combos');
    expect(ed.element.textContent!.includes('Connect a table'), true, 'saying which of the two states it is');

    const uncomputed = aggregationEditor(stringParam('aggregations'), fakeCtx(null, true).ctx as never);
    uncomputed.setValue('');
    expect(uncomputed.element.textContent!.includes('not been computed'), true);

    const stored = aggregationEditor(stringParam('aggregations'), fakeCtx(null, false).ctx as never);
    stored.setValue(JSON.stringify([{column: 'score', type: 'avg'}]));
    expect(stored.element.textContent!.includes('avg(score)'), true, 'the stored list is shown');
  });

  test('the aggregation editor renders one row per aggregation and adds new ones', async () => {
    const columns = [
      DG.Column.fromStrings('team', ['a']),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'score', [1]),
    ];
    const {ctx} = fakeCtx(columns, true);
    const ed = aggregationEditor(stringParam('aggregations'), ctx as never);
    let stored = '';
    ed.onChanged = (v): void => {stored = String(v ?? '');};

    ed.setValue(JSON.stringify([{column: 'score', type: 'avg'}]));
    expect(ed.element.hasAttribute('data-blocked'), false, 'with columns it is a real form');
    expect(ed.element.querySelectorAll('.ff-aggregation-row').length, 1, 'one row per aggregation');
    // Aggregation + column = two combos on a normal row.
    expect(ed.element.querySelectorAll('.ff-aggregation-row select').length, 2);

    (ed.element.querySelector('.ff-aggregation-add') as HTMLElement).click();
    expect(ed.element.querySelectorAll('.ff-aggregation-row').length, 2, 'Add appends a row');
    expect(parseAggregations(stored).length, 2, 'and the value is reported');

    (ed.element.querySelectorAll('.ff-aggregation-remove')[1] as HTMLElement).click();
    expect(ed.element.querySelectorAll('.ff-aggregation-row').length, 1, 'the remove icon removes it');
    expect(parseAggregations(ed.getValue()).length, 1);
  });

  test('a column-less aggregation drops its column combo', async () => {
    const columns = [DG.Column.fromList(DG.COLUMN_TYPE.INT, 'score', [1])];
    const ed = aggregationEditor(stringParam('aggregations'), fakeCtx(columns, true).ctx as never);
    ed.setValue(JSON.stringify([{column: '', type: 'count'}]));
    expect(ed.element.querySelectorAll('.ff-aggregation-row select').length, 1,
      'count shows the aggregation combo only');
  });

  test('the expression editor falls back to an editable string input with no table', async () => {
    const {ctx} = fakeCtx(null, false);
    const ed = rowConditionEditor(stringParam('condition'), ctx as never);
    let stored = '';
    ed.onChanged = (v): void => {stored = String(v ?? '');};
    ed.setValue('${age} > 30');

    expect(ed.element.getAttribute('data-mode'), 'plain', 'no table means the plain input');
    const input = ed.element.querySelector('input') as HTMLInputElement;
    expect(input !== null, true, 'and it is a real, editable input');
    expect(input.value, '${age} > 30', 'showing the stored condition');
    expect(ed.element.textContent!.includes('Connect a table'), true, 'with the reason spelled out');

    input.value = '${age} > 40';
    input.dispatchEvent(new Event('input'));
    input.dispatchEvent(new Event('change'));
    expect(ed.getValue(), '${age} > 40', 'edits are kept');
    expect(stored, '${age} > 40', 'and reported');
  });

  test('the expression editor distinguishes unconnected from uncomputed', async () => {
    const ed = rowConditionEditor(stringParam('condition'), fakeCtx(null, true).ctx as never);
    ed.setValue('');
    expect(ed.element.textContent!.includes('not been computed'), true,
      'a wired-but-uncomputed table gets a different instruction');
    expect(ed.element.querySelector('.ff-expression-editor-load') !== null, true);
  });

  test('with a table the expression editor mounts the real formula editor', async () => {
    const table = sampleTable();
    const ed = rowConditionEditor(stringParam('condition'), {
      inputValue: () => undefined,
      columns: () => Array.from(table.columns),
      table: () => table,
      isConnected: () => true,
      watch: () => {},
    } as never);
    let reported = '';
    ed.onChanged = (v): void => {reported = String(v ?? '');};
    ed.setValue('${score} > 3');

    const mounted = await until(() => ed.element.getAttribute('data-mode') === 'formula', 20_000);
    // skip ONLY when PowerPack's editor function is absent — if present, the mount must happen
    if (!mounted && DG.Func.find({package: 'PowerPack', name: 'expressionEditorWidget'}).length === 0) {
      expect(ed.element.getAttribute('data-mode'), 'plain',
        'no PowerPack editor function — it must still fall back to a usable input');
      console.warn('Flow: skipping the live formula-editor mount — PowerPack:expressionEditorWidget ' +
        'is not on this stand. Republish PowerPack to cover it.');
      return;
    }
    expect(mounted, true, 'the formula editor mounted');
    // CodeMirror arrives a beat after the host flips to `formula` — the widget builds it asynchronously
    const cmReady = await until(() => ed.element.querySelector('.cm-editor') !== null, 20_000);
    expect(cmReady, true, 'CodeMirror is mounted');
    expect(ed.getValue(), '${score} > 3', 'seeded with the stored condition');
    expect(reported, '', 'and merely opening it is not an edit');
    ed.detach!();
  });

  test('the expression editor releases its subscription on detach', async () => {
    const ed = rowConditionEditor(stringParam('condition'), fakeCtx(null, false).ctx as never);
    ed.setValue('x');
    expect(typeof ed.detach, 'function', 'it declares cleanup');
    ed.detach!();
    expect(ed.element.hasAttribute('data-mode'), false, 'and the mounted state is torn down');
  });
});
