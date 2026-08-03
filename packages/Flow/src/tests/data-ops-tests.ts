/** Flow's own data operations — the row/column/aggregation verbs written to
 *  replace the platform originals whose predicates are `FuncCall`-typed.
 *
 *  Three things are worth locking down, and they map onto the three categories
 *  below: the functions compute what they claim (against a real backend — the
 *  condition engine IS `core:AddNewColumn`), the nodes they produce cannot be
 *  run half-configured, and the editors that make those parameters fillable
 *  degrade to something usable when the upstream table isn't there yet. */

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

/** A small, deterministic table: three groups, a numeric and a string column. */
function sampleTable(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('team', ['a', 'a', 'b', 'b', 'c', 'c']),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'score', [1, 5, 2, 8, 3, 9]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'weight', [10, 20, 30, 40, 50, 60]),
  ]);
}

// ---------- the functions ----------

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
    // The source is never touched — a pipeline step that mutated its input
    // would make an upstream preview show downstream data.
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
    // No columns chosen means every column — a blank list must not produce an
    // empty table.
    const all = await extractRows(table, '${team} == "b"', null);
    expect(all.columns.length, 3, 'a blank column list keeps them all');
  });

  test('Select Rows selects on the table and passes it on', async () => {
    const table = sampleTable();
    const result = await selectRows(table, '${score} > 3', true);
    expect(result === table, true, 'the same table flows on — selection is table state');
    expect(table.selection.trueCount, 3);
    // Selecting again with `clearSelection` must replace, not accumulate.
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
    // The point of the check: a numeric expression silently matching nothing
    // (or everything) is far worse than a refusal that says what to write.
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
    // A pipeline that returns different rows on every re-run breaks
    // invalidation: a node re-runs and its downstream no longer agrees with it.
    const a = randomIndices(100, 10, 42);
    const b = randomIndices(100, 10, 42);
    expect(a.join(','), b.join(','), 'the same seed draws the same rows');
    expect(new Set(a).size, 10, 'no duplicates');
    expect(a.every((i) => i >= 0 && i < 100), true, 'all in range');
    expect(randomIndices(100, 10, 7).join(',') !== a.join(','), true, 'a different seed differs');
    // Asking for more rows than exist yields the whole table, not a crash.
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
    // A pivot column turns the group-by into a pivot table — the whole reason
    // this node also answers to "pivot".
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
    // The list is hand-written and includes `concat unique`, a STR_AGG rather
    // than an AGG — this is what proves the Dart side takes it.
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

// ---------- the nodes ----------

category('Flow: data op nodes', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  /** Every one of these must be in the catalog without an allowlist entry —
   *  they declare `meta.includeInFlow: true` on themselves. */
  const OPS = ['filterRows', 'deleteRows', 'extractRows', 'selectRows', 'filterRandomRows',
    'selectRandomRows', 'deleteColumns', 'tagColumns', 'expressionToColumn', 'aggregate', 'unpivot'];

  test('all of them are registered as nodes', async () => {
    for (const name of OPS)
      expect(typeNameOf(`Flow:${name}`) !== null, true, `Flow:${name} is a node type`);
  });

  test('the FuncCall-typed originals they replace stay out of the catalog', async () => {
    // These are the entries whose `TableRowFilterCall` / `ColFilterCall`
    // parameters have no editor and no socket — the reason the whole family was
    // commented out in included-funcs.ts.
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
        // Blank on a fresh node — the whole point of the requirement.
        expect(nodeMissingRequirements(node, () => false).length > 0, true,
          `${name}: a fresh node is not runnable`);
      } finally {
        destroyEditor(e);
      }
    }
  });

  test('column parameters are real column_list slots bound to the table', async () => {
    // `column_list` is what gives them Flow's column picker for free — the
    // originals took `ColFilterCall` and `core:Unpivot` takes bare string lists.
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
    // No category is declared anywhere: `categorizeFunc` routes by signature, so
    // this is the check that the signatures themselves are shaped right.
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
      // The one that produces a column rather than a table.
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
    // A `//friendlyName:` annotation does not survive publishing (checked on a
    // live stand), so the header has to humanize the raw name itself.
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

      // A blank value is already caught by the generic required-input check;
      // what needs the validator is a list that is non-blank but says nothing.
      node.inputValues['aggregations'] = JSON.stringify([{column: '', type: 'avg'}]);
      expect(validator!(node).length > 0, true, 'an aggregation over no column blocks the node');
      expect(validator!(node)[0].includes('column'), true, 'and names what is missing');

      node.inputValues['aggregations'] = JSON.stringify([{column: 'score', type: 'avg'}]);
      expect(validator!(node).length, 0, 'a complete list is ready');

      // `count` counts rows, so it legitimately has no column.
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
      // The verdict channel: deciding whether a formula is true/false means
      // evaluating it, so the hosted editor announces each check on the element
      // it is mounted in and the node reads the answer back synchronously.
      const validate = (expression: string, error: string): void => {
        ed.element.dispatchEvent(new CustomEvent('expression-validated', {detail: {expression, error}}));
      };

      validate('${score} + 1', 'The formula is int, not a true/false condition. Compare something');
      expect(ed.element.hasAttribute('data-invalid'), true, 'the editor marks itself');
      const problems = validator!(node);
      expect(problems.length, 1, 'and the node refuses to run');
      expect(problems[0].includes('true/false'), true, 'saying why');

      // A verdict is about ONE piece of text; editing past it makes the state
      // unknown again rather than leaving a stale block in place.
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

// ---------- the editors ----------

category('Flow: data op editors', () => {
  /** Stand-in for the panel context, so an editor's states are testable without
   *  a live node — same shape the MPO mapper's tests use. */
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
    // Two broken entries of the same kind read as one problem, not two.
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
    // A blank column is a row just added, not a stale reference — dropping it
    // would make the Add button look broken.
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

    // Connected but uncomputed is a different instruction from unconnected.
    const uncomputed = aggregationEditor(stringParam('aggregations'), fakeCtx(null, true).ctx as never);
    uncomputed.setValue('');
    expect(uncomputed.element.textContent!.includes('not been computed'), true);

    // Whatever is stored stays readable — a flow loaded from disk must not look
    // empty just because it has not been run.
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

    // Editable is the point: a flow loaded from disk must be fixable without
    // running anything first.
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
    // …and an action to fix it, rather than a dead end.
    expect(ed.element.querySelector('.ff-expression-editor-load') !== null, true);
  });

  test('with a table the expression editor mounts the real formula editor', async () => {
    // The whole point of the parameter shape: with a table in hand the
    // condition is edited in PowerPack's formula editor — column autocomplete,
    // function list, inline validation — not a text box.
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
    // Narrow skip: ONLY when the editor function is absent from this stand's
    // catalog. If it IS there the mount has to happen — a silent fallback would
    // otherwise hide a broken hand-off behind a passing test.
    if (!mounted && DG.Func.find({package: 'PowerPack', name: 'expressionEditorWidget'}).length === 0) {
      expect(ed.element.getAttribute('data-mode'), 'plain',
        'no PowerPack editor function — it must still fall back to a usable input');
      console.warn('Flow: skipping the live formula-editor mount — PowerPack:expressionEditorWidget ' +
        'is not on this stand. Republish PowerPack to cover it.');
      return;
    }
    expect(mounted, true, 'the formula editor mounted');
    // CodeMirror is what the widget mounts; its presence is the proof that the
    // real editor came up rather than an empty host. It arrives a beat after
    // the host flips to `formula` — the widget builds it asynchronously.
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
