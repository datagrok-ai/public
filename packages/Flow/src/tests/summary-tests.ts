/** Tests for the heuristic node/flow summaries (U12) — pure, no platform. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {summarizeNode, summarizeFlow, humanize} from '../summary/summary-generator';
import {FlowNode, FlowConnection} from '../rete/scheme';
import {registerBuiltinNodes, createNode} from '../rete/node-factory';

function funcNode(name: string, typeName: string, inputValues: Record<string, unknown> = {}): FlowNode {
  const n = new FlowNode('node');
  n.dgFuncName = name;
  n.dgTypeName = typeName;
  n.dgNodeType = 'func';
  Object.assign(n.inputValues, inputValues);
  return n;
}
const dgType = (name: string): string => `DG Functions/Uncategorized/${name}`;

category('Flow: summaries', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('humanize turns identifiers into phrases', async () => {
    expect(humanize('addChemRisksColumns'), 'Add chem risks columns');
    expect(humanize('FilterRows'), 'Filter rows');
    expect(humanize('OpenFile'), 'Open file');
    expect(humanize('to_atomic_level'), 'To atomic level');
  });

  test('built-in node summaries (input/output)', async () => {
    const inp = createNode('Inputs/Table Input')!;
    const out = createNode('Outputs/Table Output')!;
    expect(summarizeNode(inp), 'Input “df” (dataframe)');
    expect(summarizeNode(out), 'Outputs “result”');
  });

  test('curated function summaries read the node\'s inputs', async () => {
    const open = funcNode('OpenFile', dgType('OpenFile'), {fullPath: 'System:DemoFiles/demog.csv'});
    expect(summarizeNode(open), 'Loads file demog.csv');

    const col = funcNode('AddNewColumn', dgType('AddNewColumn'), {name: 'AgeMonths', expression: '${AGE} * 12'});
    expect(summarizeNode(col).startsWith('Adds column “AgeMonths”'), true, summarizeNode(col));

    const logp = funcNode('CalculateLogP', dgType('CalculateLogP'));
    expect(summarizeNode(logp), 'Calculates logP');

    const chem = funcNode('addChemPropertiesColumns', dgType('addChemPropertiesColumns'),
      {MW: true, logP: true, HBA: true});
    expect(summarizeNode(chem), 'Computes MW, HBA, logP');
  });

  test('fallback humanizes unknown functions', async () => {
    const n = funcNode('fooBarBaz', dgType('fooBarBaz'));
    expect(summarizeNode(n), 'Foo bar baz');
  });

  test('fallback prefers the description, then a distinct friendly name verbatim', async () => {
    const withDesc = funcNode('getInchis', dgType('getInchis'));
    withDesc.dgFunc = {name: 'getInchis', friendlyName: 'InChI',
      description: 'Computes the InChI identifier for each molecule.'} as never;
    expect(summarizeNode(withDesc), 'Computes the InChI identifier for each molecule.');

    const friendlyOnly = funcNode('getInchis', dgType('getInchis'));
    friendlyOnly.dgFunc = {name: 'getInchis', friendlyName: 'InChI', description: ''} as never;
    expect(summarizeNode(friendlyOnly), 'InChI');

    const noMeta = funcNode('fooBarBaz', dgType('fooBarBaz'));
    noMeta.dgFunc = {name: 'fooBarBaz', friendlyName: 'fooBarBaz', description: ''} as never;
    expect(summarizeNode(noMeta), 'Foo bar baz');

    const curated = funcNode('CalculateLogP', dgType('CalculateLogP'));
    curated.dgFunc = {name: 'CalculateLogP', friendlyName: 'CalculateLogP',
      description: 'Some long docs text.'} as never;
    expect(summarizeNode(curated), 'Calculates logP');
  });

  test('summarizeFlow groups disjoint pipelines, ordered top-then-left', async () => {
    const a = funcNode('OpenFile', dgType('OpenFile'), {fullPath: 'x/demog.csv'});
    a.pos = {x: 0, y: 0};
    const b = funcNode('AddNewColumn', dgType('AddNewColumn'), {name: 'c', expression: '1'});
    b.pos = {x: 240, y: 0};
    const c = funcNode('CalculateLogP', dgType('CalculateLogP'));
    c.pos = {x: 0, y: 300};

    const conns = [{source: a.id, target: b.id, targetInput: 'table'}] as unknown as FlowConnection[];
    const s = summarizeFlow([a, b, c], conns);

    expect(s.pipelineCount, 2, s.text);
    expect(s.nodeCount, 3);
    expect(s.pipelines[0].steps.length, 2, JSON.stringify(s.pipelines[0]));
    expect(s.pipelines[0].steps[0].caption.includes('Loads file'), true, s.pipelines[0].steps[0].caption);
    expect(s.pipelines[1].steps.length, 1, JSON.stringify(s.pipelines[1]));
  });

  test('summarizeFlow orders by dependency and records what feeds what', async () => {
    // b is placed LEFT of a so a positional sort would mis-order.
    const a = funcNode('OpenFile', dgType('OpenFile'), {fullPath: 'x/demog.csv'});
    a.pos = {x: 500, y: 0};
    const b = funcNode('AddNewColumn', dgType('AddNewColumn'), {name: 'c', expression: '1'});
    b.pos = {x: 0, y: 0};
    const conns = [{source: a.id, target: b.id, targetInput: 'table'}] as unknown as FlowConnection[];
    const s = summarizeFlow([a, b], conns);

    expect(s.pipelines[0].steps[0].caption.includes('Loads file'), true, 'producer first');
    expect(s.pipelines[0].steps[1].caption.includes('Adds column'), true, 'consumer second');
    expect(s.pipelines[0].steps[1].inputs.length, 1);
    expect(s.pipelines[0].steps[1].inputs[0].key, 'table');
    expect(s.pipelines[0].steps[1].inputs[0].from, 1);
    expect(s.pipelines[0].steps[0].inputs.length, 0, 'source has no inputs');
  });

  test('summarizeFlow does not truncate node captions', async () => {
    const longExpr = '${AGE} * 12 + ${HEIGHT} / 100 - ${WEIGHT} * 2 + somethingQuiteLong';
    const n = funcNode('AddNewColumn', dgType('AddNewColumn'), {name: 'verbose', expression: longExpr});
    n.pos = {x: 0, y: 0};
    const s = summarizeFlow([n], []);
    expect(s.pipelines[0].steps[0].caption.includes(longExpr), true, s.pipelines[0].steps[0].caption);
  });

  test('summarizeFlow on an empty graph is graceful', async () => {
    const s = summarizeFlow([], []);
    expect(s.pipelineCount, 0);
    expect(s.nodeCount, 0);
  });
});
