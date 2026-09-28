import * as grok from 'datagrok-api/grok';
import {v4 as uuidv4} from 'uuid';
import {BaseTree, NodePath, NodePathSegment} from '../data/BaseTree';
import {isFuncCallNode, StateTreeNode} from './StateTreeNodes';
import {LinkSpec, MatchInfo, matchNodeLink} from './link-matching';
import {DriverLogger, reportError} from '../data/Logger';
import {Link} from './Link';
import {parseLinkIO} from '../config/LinkSpec';
import {ruleMetaHandler, ruleValidatorHandler} from './rule-handlers';
import {CALL, CheckOptions, expandChecks, TABLE, TARGET, VALUE} from '../config/checks';

export class DependenciesData {
  nodes: Set<string> = new Set();
  links: Set<string> = new Set();
}

export interface IoDep {
  data?: string;
}

export type IoDeps = Record<string, IoDep>;

export function calculateStepsDependencies(state: BaseTree<StateTreeNode>, links: Link[]) {
  const deps: Map<string, DependenciesData> = new Map();
  for (const link of links) {
    const inputInfos = Object.values(link.matchInfo.inputs);
    for (const infosIn of inputInfos) {
      for (const infoIn of infosIn) {
        const inPathFull = [...link.prefix, ...infoIn.path];
        addOutputDeps(state, deps, link, inPathFull);
      }
    }
    if (inputInfos.length === 0)
      addOutputDeps(state, deps, link);
  }
  return deps;
}

function addOutputDeps(state: BaseTree<StateTreeNode>, deps: Map<string, DependenciesData>, link: Link, inPathFull?: NodePathSegment[]) {
  for (const infosOut of Object.values(link.matchInfo.outputs)) {
    for (const infoOut of infosOut) {
      const outPathFull = [...link.prefix, ...infoOut.path];
      const nodeOut = state.getNode(outPathFull);
      const depData = deps.get(nodeOut.getItem().uuid) ?? new DependenciesData();
      if (inPathFull && !BaseTree.isNodeAddressEq(inPathFull, outPathFull)) {
        const nodeIn = state.getNode(inPathFull);
        depData.nodes.add(nodeIn.getItem().uuid);
      }
      depData.links.add(link.uuid);
      deps.set(nodeOut.getItem().uuid, depData);
    }
  }
}

export function calculateIoDependencies(state: BaseTree<StateTreeNode>, links: Link[], logger?: DriverLogger) {
  const deps = new Map<string, IoDeps>();
  for (const link of links) {
    const linkId = link.uuid;
    for (const infosOut of Object.values(link.matchInfo.outputs)) {
      for (const infoOut of infosOut) {
        const stepPath = [...link.prefix, ...infoOut.path];
        const node = state.getNode(stepPath);
        const ioName = infoOut.ioName!;
        const depsData = deps.get(node.getItem().uuid) ?? {};
        const depType = link.matchInfo.spec.type ?? 'data';
        if (depsData[ioName] == null)
          depsData[ioName] = {};
        if (depType === 'data') {
          if (depsData[ioName][depType]) {
            const prevUuid = depsData[ioName][depType] as string;
            const prevSpecId = links.find((l) => l.uuid === prevUuid)?.matchInfo.spec.id ?? prevUuid;
            const newSpecId = link.matchInfo.spec.id;
            reportError(
              'warning',
              'ioDependencies',
              `Duplicate data link target: step ${JSON.stringify(stepPath)} io "${ioName}" is written by both "${prevSpecId}" and "${newSpecId}"`,
              logger,
              [prevSpecId, newSpecId],
            );
          }
          depsData[ioName][depType] = linkId;
        }
        deps.set(node.getItem().uuid, depsData);
      }
    }
  }
  return deps;
}

export function createDefaultValidators(state: BaseTree<StateTreeNode>, logger?: DriverLogger) {
  const defaultValidators = state.traverse(state.root, (acc, node, path) => {
    const item = node.getItem();
    if (!isFuncCallNode(item))
      return acc;
    const ios = item.config.io ?? [];
    const validators = ios.flatMap((io) => {
      if (io.direction === 'output')
        return [];
      const options: CheckOptions = {...io.checks, nullable: io.nullable};
      // the annotation's validators are run by the platform, through the step's FuncCall
      const annotationValidators = options.validators;
      delete options.validators;
      const tableIo = options.table == null ? undefined :
        ios.find((other) => other.id === options.table && other.direction === 'input');
      if (!tableIo)
        delete options.table;
      const expanded = expandChecks(options);
      if (annotationValidators?.length) {
        expanded.push({
          key: 'validators', family: 'validator', needsTable: false, needsCall: true, needsInputs: false,
          params: {
            when: {'!': {missing: [VALUE]}},
            sources: {verdicts: {validators: {input: VALUE, call: CALL}}},
            effects: [{effect: 'verdicts', targets: [TARGET], source: 'verdicts'}],
          },
        });
      }
      // a GrokScript expression sees every input of the step under its own name
      const reserved = new Set([VALUE, TABLE, TARGET, CALL]);
      const stepInputs = ios.filter((other) => other.direction === 'input' && !reserved.has(other.id));
      return expanded.map(({key, family, needsTable, needsCall, needsInputs, params}) => {
        const spec: LinkSpec = {
          id: `::${io.id}:${key}`,
          from: [
            ...parseLinkIO(`${VALUE}:${io.id}`, 'input'),
            ...(needsTable ? parseLinkIO(`${TABLE}:${tableIo!.id}`, 'input') : []),
            ...(needsCall ? parseLinkIO(`${CALL}(call,optional):.`, 'input') : []),
            ...(needsInputs ? stepInputs.flatMap((other) => parseLinkIO(`${other.id}:${other.id}`, 'input')) : []),
          ],
          to: parseLinkIO(`${TARGET}:${io.id}`, 'output'),
          type: family,
          handler: family === 'meta' ? ruleMetaHandler : ruleValidatorHandler,
          params,
        } as LinkSpec;
        const inputs: MatchInfo['inputs'] = {[VALUE]: [{path: [], ioName: io.id}]};
        if (needsTable)
          inputs[TABLE] = [{path: [], ioName: tableIo!.id}];
        if (needsCall)
          inputs[CALL] = [{path: []}];
        if (needsInputs) {
          for (const other of stepInputs)
            inputs[other.id] = [{path: [], ioName: other.id}];
        }
        const minfo: MatchInfo = {
          spec,
          inputs,
          outputs: {[TARGET]: [{path: [], ioName: io.id}]},
          actions: {},
          inputsUUID: new Map(),
          outputsUUID: new Map(),
          isDefaultValidator: true,
        };
        return new Link(path, minfo, 0, logger);
      });
    });
    return [...acc, ...validators];
  }, [] as Link[]);
  return defaultValidators;
}
