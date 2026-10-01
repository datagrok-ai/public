import * as DG from 'datagrok-api/dg';
import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {makeFuncCall} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {resolveSources} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

const VALUES = 'LibTests:TestValueAnnotations';
const LOOKUP = 'LibTests:TestLookupAnnotations';
const created = ['calc', 'bare', 'metric', 'speed', 'region', 'city', 'model'];

async function newCallInputs(initialValues?: Record<string, any>) {
  const {adapter} = await makeFuncCall(VALUES, false, initialValues);
  const fc = adapter.getFuncCall();
  return Object.fromEntries(created.map((name) => [name, fc.inputs[name]]));
}

category('ComputeUtils: Driver value annotations', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('New calls get evaluated and literal defaults, choices stay empty', async () => {
    expectDeepEqual(await newCallInputs(), {
      calc: 4, bare: 'high', metric: 'minkowski', speed: null, region: 'EU', city: null, model: null,
    });
  });

  test('Initial values win over defaults', async () => {
    expectDeepEqual(await newCallInputs({calc: 7, bare: null, speed: 'fast'}), {
      calc: 7, bare: null, metric: 'minkowski', speed: 'fast', region: 'EU', city: null, model: null,
    });
  });

  test('Choices source evaluates once per dependency change', async () => {
    const fc = DG.Func.byName(VALUES).prepare({region: 'EU'});
    const choices = (io: string, value: any): any => resolveSources({
      hasCall: (name: string) => name === 'call',
      getFirst: (name: string) => name === 'call' ? fc : value,
      getMatchedPositions: () => [{path: [], position: 0, ioName: io}],
    } as any, {c: {choices: {input: 'value', call: 'call'}}});
    const first = await choices('city', 'EU-1');
    expectDeepEqual(first.c, {
      items: ['EU-1', 'EU-2'], values: {'EU-1': 'EU-1', 'EU-2': 'EU-2'}, inList: true, row: null, rowErrors: [],
    });
    const cached = choices('city', 'Nope');
    expect(cached instanceof Promise, false);
    expect(cached.c.inList, false);
    fc.setParamValue('region', 'US');
    const refreshed = choices('city', 'EU-1');
    expect(refreshed instanceof Promise, true);
    const {c} = await refreshed;
    expectDeepEqual([c.items, c.inList], [['US-1', 'US-2'], false]);
    expectDeepEqual((await choices('model', 'Volvo')).c.row, {mpg: 30, CYL: 4});
    expect((await choices('model', 'Nope')).c.row === null, true);
    expect(choices('model', null).c.inList, true);
  });

  test('Lookup rows are converted to the input types', async () => {
    const fc = DG.Func.byName(LOOKUP).prepare();
    const lookup = async (value: string) => (await resolveSources({
      hasCall: (name: string) => name === 'call',
      getFirst: (name: string) => name === 'call' ? fc : value,
      getMatchedPositions: () => [{path: [], position: 0, ioName: 'model'}],
    } as any, {c: {choices: {input: 'value', call: 'call'}}})).c;
    const mazda = await lookup('Mazda');
    const volvo = await lookup('Volvo');
    expectDeepEqual([mazda.row, mazda.rowErrors], [
      {cyl: 4, name: '1', flag: true, engine: 'E1'}, ['mpg: 21.5 is not a valid int'],
    ]);
    expectDeepEqual([volvo.row, volvo.rowErrors], [
      {mpg: 30, name: '2', flag: false, engine: 'E2'}, ['cyl: "abc" is not a valid int'],
    ]);
  });

  test('Only the first propagateChoice key gets a lookup', async () => {
    const logger = new DriverLogger();
    const config: PipelineConfiguration = {id: 'pipeline1', type: 'static', steps: [{id: 'step', nqName: LOOKUP}]};
    const pconf: any = await getProcessedConfig(config, logger);
    const links = pconf.steps[0].links;
    expectDeepEqual(links.map((link: any) => link.id), [
      '::model:choices::meta', '::model:choices::validator', '::model:lookup::data',
      '::engine:choices::meta', '::engine:choices::validator',
    ]);
    const byId = Object.fromEntries(links.map((link: any) => [link.id, link]));
    expectDeepEqual(byId['::model:choices::validator'].params.effects[1], {
      effect: 'warning', targets: 'model_target', message: {var: 'model_choices.rowErrors'},
      when: {'!!': {var: 'model_choices'}},
    });
    expectDeepEqual(byId['::engine:choices::validator'].params.effects.map((effect: any) => effect.message),
      ['Not in the list of choices']);
    expectDeepEqual(byId['::model:lookup::data'].to.map((item: any) => item.name),
      ['engine', 'cyl', 'mpg', 'name', 'flag']);
    expectDeepEqual(logger.errors.map((item) => [item.severity, item.message]), [['warning',
      `Step ${LOOKUP}: propagateChoice on 'engine' is ignored, 'model' already fills the step's inputs`]]);
  });

  test('Dynamic choices become rules on the step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{id: 'step', nqName: VALUES}],
      links: [{
        id: 'fixMpg', from: 'in:step/calc', to: 'out:step/mpg',
        handler({controller}) {
          controller.setAll('out', controller.getFirst('in'));
        },
      }],
    };
    const pconf: any = await getProcessedConfig(config);
    const io = Object.fromEntries(pconf.steps[0].io.map((item: any) => [item.id, item]));
    expectDeepEqual([io.city.dynamicChoices, io.model.dynamicChoices, io.speed.dynamicChoices, io.speed.checks],
      [{propagate: false}, {propagate: true}, undefined, {choices: ['slow', 'fast']}]);
    const links = pconf.steps[0].links;
    expectDeepEqual(links.map((link: any) => [link.id, link.type]), [
      ['::city:choices::meta', 'meta'], ['::city:choices::validator', 'validator'],
      ['::model:choices::meta', 'meta'], ['::model:choices::validator', 'validator'], ['::model:lookup::data', 'data'],
    ]);
    const byId = Object.fromEntries(links.map((link: any) => [link.id, link]));
    expectDeepEqual(byId['::city:choices::meta'].params.effects, [
      {effect: 'items', targets: 'city_target', items: {var: 'city_choices.items'}, when: {'!!': {var: 'city_choices'}}},
      {effect: 'meta', targets: 'city_target', meta: {emptyChoice: true}},
    ]);
    expectDeepEqual(byId['::city:choices::validator'].debounce, 0);
    expectDeepEqual(byId['::model:choices::validator'].params.effects[0].message, 'Not in the lookup table');
    const lookup = byId['::model:lookup::data'];
    expectDeepEqual([lookup.runOnInit, lookup.from.map((item: any) => item.name)], [true, ['model', 'call']]);
    expectDeepEqual(lookup.to.map((item: any) => item.name), ['calc', 'bare', 'metric', 'speed', 'region', 'city', 'mpg', 'cyl']);
    expectDeepEqual(lookup.params.effects, [{
      effect: 'assign', values: {var: 'model_choices.row'}, ignoreCase: true, restriction: 'restricted',
      when: {'!!': {var: 'model_choices.row'}},
    }]);
    let outputs: string[] = [];
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const link = [...tree.linksState.links.values()].find((item) => item.matchInfo.spec.id === '::model:lookup::data')!;
      outputs = Object.keys(link.matchInfo.outputs).sort();
    });
    expectDeepEqual(outputs, ['bare', 'calc', 'city', 'cyl', 'metric', 'region', 'speed']);
  });
});
