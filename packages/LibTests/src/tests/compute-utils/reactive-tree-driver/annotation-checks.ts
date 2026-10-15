import {category, test, before, awaitCheck} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {Link} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/Link';
import {isAnnotationCheck} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {inspectLinks} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/inspection';
import {DriverLogger, isLinkLogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, errors} from '../../../test-utils';

const ANNOTATED = 'LibTests:TestAnnotatedInputs';
const NAMED = 'LibTests:TestNamedValidators';
const VALUES = 'LibTests:TestValueAnnotations';

const annotatedStep = (links: PipelineLinkConfigurationInput<string | string[]>[] = []): PipelineConfiguration =>
  ({id: 'p', type: 'static', steps: [{id: 's', nqName: ANNOTATED}], links});

const checkLinks = (tree: StateTree) =>
  [...tree.linksState.links.values()].filter((link) => isAnnotationCheck(link.matchInfo.spec));

const nodeUuid = (tree: StateTree, link: Link) => tree.nodeTree.getNode(link.prefix).getItem().uuid;

const linksOf = (tree: StateTree, uuid: string) =>
  checkLinks(tree).filter((link) => nodeUuid(tree, link) === uuid).map((link) => link.uuid);

const stepAt = (tree: StateTree, ...idx: number[]) =>
  tree.nodeTree.getNode(idx.map((i) => ({idx: i}))).getItem() as FuncCallNode;

const sameJson = (a: any, b: any) => JSON.stringify(a) === JSON.stringify(b);

category('ComputeUtils: Driver annotation checks', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Annotation checks follow the annotationChecks option', async () => {
    const pconf = await getProcessedConfig(annotatedStep());
    const validationOfV = (annotationChecks: boolean) => {
      let result: any;
      testScheduler.run(({cold}) => {
        const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks});
        StateTree.loadOrCreateCalls(tree, true).subscribe();
        tree.init().subscribe();
        cold('-a').subscribe(() => stepAt(tree, 0).getStateStore().setState('v', 20));
        cold('--a').subscribe(() => result = stepAt(tree, 0).validationInfo$.value.v);
      });
      return result;
    };
    expectDeepEqual(validationOfV(false) === undefined, true, {prefix: 'Without annotation checks'});
    expectDeepEqual(validationOfV(true), errors('Must be at most 10'), {prefix: 'With annotation checks'});
  });

  test('Annotation checks keep their link ids and kinds', async () => {
    const pconf = await getProcessedConfig(annotatedStep());
    let links: [string, string | undefined][] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      tree.init().subscribe();
      cold('-a').subscribe(() =>
        links = checkLinks(tree).map((link) => [link.matchInfo.spec.id, link.matchInfo.spec.annotation]));
    });
    const required = ['c', 'code', 'col', 'df', 'mode', 'mol', 'v'].map((io) => [`::${io}:required`, 'required']);
    const checks = ['code:validator', 'col:allowNulls', 'col:type', 'col:validators', 'mol:semType', 'v:max', 'v:min']
      .map((key) => [`::${key}`, 'check']);
    expectDeepEqual(links.sort(), [...required, ...checks].sort());
  });

  test('Annotation checks are hidden from links info and the link log', async () => {
    const pconf = await getProcessedConfig(annotatedStep());
    let info: string[] = [];
    let ids = new Set<string>();
    const logger = new DriverLogger();
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true, logger});
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        ids = new Set(checkLinks(tree).map((link) => link.matchInfo.spec.id));
        info = tree.linksState.getLinksInfo().map((link) => link.id);
      });
    });
    expectDeepEqual(ids.size > 0, true, {prefix: 'Annotation checks exist'});
    expectDeepEqual(info.filter((id) => ids.has(id)), [], {prefix: 'Links info'});
    const logged = logger.logs$.value.filter(isLinkLogItem).filter((item) => ids.has(item.id));
    expectDeepEqual(logged.filter((item) => item.type === 'linkAdded' || item.type === 'linkRemoved'), [],
      {prefix: 'Added and removed'});
    expectDeepEqual(logged.length > 0 && logged.every((item) => item.annotation === 'check' ||
      item.annotation === 'required'), true, {prefix: 'Runs carry the kind'});
  });

  test('Inspection lists annotation links with their kind', async () => {
    const pconf = await getProcessedConfig({
      id: 'p', type: 'static', steps: [{id: 'checked', nqName: ANNOTATED}, {id: 'chosen', nqName: VALUES}],
    });
    let inspected: Record<string, string | undefined> = {};
    let info: string[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        inspected = Object.fromEntries(inspectLinks(tree).map((link) => [link.id, link.annotation]));
        info = tree.linksState.getLinksInfo().map((link) => link.id);
      });
    });
    expectDeepEqual([inspected['::v:max'], inspected['::v:required'], inspected['::city:choices::meta']],
      ['check', 'required', 'rule']);
    expectDeepEqual(info.includes('::city:choices::meta'), true, {prefix: 'Rules stay in links info'});
    expectDeepEqual(info.some((id) => inspected[id] === 'check' || inspected[id] === 'required'), false,
      {prefix: 'Checks stay out of links info'});
  });

  test('Annotation checks run on nested and added steps', async () => {
    const pconf = await getProcessedConfig({
      id: 'root',
      type: 'static',
      steps: [
        {id: 'nested', type: 'static', steps: [{id: 'inner', nqName: ANNOTATED}]},
        {id: 'items', type: 'parallel', stepTypes: [{id: 'item', nqName: ANNOTATED}], initialSteps: []},
      ],
    });
    const results: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        const items = tree.nodeTree.getNode([{idx: 1}]).getItem();
        tree.addSubTree(items.uuid, 'item', 0).subscribe();
      });
      cold('100ms a').subscribe(() => {
        stepAt(tree, 0, 0).getStateStore().setState('v', 20);
        stepAt(tree, 1, 0).getStateStore().setState('v', 20);
      });
      cold('200ms a').subscribe(() =>
        results.push([stepAt(tree, 0, 0).validationInfo$.value.v, stepAt(tree, 1, 0).validationInfo$.value.v]));
    });
    expectDeepEqual(results, [[errors('Must be at most 10'), errors('Must be at most 10')]]);
  });

  test('Annotation checks are kept across tree mutations', async () => {
    const pconf = await getProcessedConfig({
      id: 'items', type: 'parallel', stepTypes: [{id: 'item', nqName: ANNOTATED}], initialSteps: ['item'],
    });
    const seen: Record<string, any> = {};
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      tree.init().subscribe();
      const rootUuid = () => tree.nodeTree.root.getItem().uuid;
      let first: FuncCallNode;
      cold('-a').subscribe(() => {
        first = stepAt(tree, 0);
        seen.initial = linksOf(tree, first.uuid);
      });
      cold('100ms a').subscribe(() => tree.addSubTree(rootUuid(), 'item', 1).subscribe());
      cold('200ms a').subscribe(() => {
        seen.afterAdd = linksOf(tree, first.uuid);
        seen.added = linksOf(tree, stepAt(tree, 1).uuid).length;
        tree.duplicateSubtree(first.uuid).subscribe();
      });
      cold('300ms a').subscribe(() => {
        seen.afterDuplicate = linksOf(tree, first.uuid);
        seen.totalAfterDuplicate = checkLinks(tree).length;
        tree.moveSubtree(first.uuid, 2).subscribe();
      });
      cold('400ms a').subscribe(() => {
        seen.moved = linksOf(tree, first.uuid).length;
        tree.removeSubtree(stepAt(tree, 0).uuid).subscribe();
      });
      cold('500ms a').subscribe(() => {
        seen.totalAfterRemove = checkLinks(tree).length;
        first.getStateStore().setState('v', 20);
      });
      cold('600ms a').subscribe(() => seen.validation = first.validationInfo$.value.v);
    });
    const perStep = seen.initial.length;
    expectDeepEqual(perStep > 0, true, {prefix: 'Checks of the first step'});
    expectDeepEqual(seen.afterAdd, seen.initial, {prefix: 'Kept after add'});
    expectDeepEqual(seen.added, perStep, {prefix: 'Added step'});
    expectDeepEqual(seen.afterDuplicate, seen.initial, {prefix: 'Kept after duplicate'});
    expectDeepEqual(seen.totalAfterDuplicate, 3 * perStep, {prefix: 'After duplicate'});
    expectDeepEqual(seen.moved, perStep, {prefix: 'Moved step'});
    expectDeepEqual(seen.totalAfterRemove, 2 * perStep, {prefix: 'After remove'});
    expectDeepEqual(seen.validation, errors('Must be at most 10'), {prefix: 'Moved step is checked'});
  });

  test('Annotation and config checks on one io keep their message order', async () => {
    const pconf = await getProcessedConfig(annotatedStep([
      {id: 'cfg', type: 'check', io: 's/v', check: {max: 5}, message: 'Config max'},
    ]));
    let result: any;
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      cold('-a').subscribe(() => stepAt(tree, 0).getStateStore().setState('v', 20));
      cold('--a').subscribe(() => result = stepAt(tree, 0).validationInfo$.value.v);
    });
    expectDeepEqual(result, errors('Config max', 'Must be at most 10'));
  });

  // the platform runs annotation validators through a real FuncCall, so this one runs in real time
  test('Annotation validators run through the step FuncCall', async () => {
    const tree = StateTree.fromPipelineConfig({
      config: await getProcessedConfig({id: 'p', type: 'static', steps: [{id: 's', nqName: NAMED}]}),
      annotationChecks: true,
    });
    await tree.init().toPromise();
    const node = stepAt(tree, 0);
    node.getStateStore().setState('x', 20);
    await awaitCheck(() => sameJson(node.validationInfo$.value.x, errors('too big')),
      `x is not too big, got ${JSON.stringify(node.validationInfo$.value.x)}`, 5000);
    node.getStateStore().setState('x', 1);
    await awaitCheck(() => node.validationInfo$.value.x == null || sameJson(node.validationInfo$.value.x, errors()),
      'x is still too big', 5000);
  });
});
