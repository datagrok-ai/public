import * as grok from 'datagrok-api/grok';
import {v4 as uuidv4} from 'uuid';
import {BaseTree, NodePath, NodePathSegment} from '../data/BaseTree';
import {isFuncCallNode, StateTreeNode} from './StateTreeNodes';
import {LinkSpec, MatchInfo, matchNodeLink} from './link-matching';
import {Link} from './Link';
import {parseLinkIO} from '../config/LinkSpec';

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

export function calculateIoDependencies(state: BaseTree<StateTreeNode>, links: Link[]) {
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
            const msg = `Duplicate deps path ${JSON.stringify(stepPath)} io ${ioName}`;
            console.error(msg);
            grok.shell.error(msg);
          }
          depsData[ioName][depType] = linkId;
        }
        deps.set(node.getItem().uuid, depsData);
      }
    }
  }
  return deps;
}

export function createDefaultValidators(state: BaseTree<StateTreeNode>) {
  const defaultValidators = state.traverse(state.root, (acc, node, path) => {
    const item = node.getItem();
    if (!isFuncCallNode(item))
      return acc;
    const validators = item.config.io?.map((io) => {
      if (io.nullable || io.direction === 'output')
        return;
      const spec: LinkSpec = {
        id: uuidv4(),
        from: parseLinkIO(`in:${io.id}`, io.direction),
        to: parseLinkIO(`out:${io.id}`, io.direction),
        type: 'validator',
        handler({controller}) {
          const val = controller.getFirst('in');
          if (val == null || val === '')
            controller.setValidation('out', ({errors: [{description: 'Missing value'}]}));
          else
            controller.setValidation('out', undefined);
        },
      };
      const minfo: MatchInfo = {
        spec,
        inputs: {
          'in': [{
            path: [],
            ioName: io.id,
          }],
        },
        outputs: {
          'out': [{
            path: [],
            ioName: io.id,
          }],
        },
        actions: {},
        inputsUUID: new Map(),
        outputsUUID: new Map(),
        isDefaultValidator: true,
      };
      return new Link(path, minfo, 0);
    }).filter((x) => !!x);
    return [...acc, ...(validators ?? [])];
  }, [] as Link[]);
  return defaultValidators;
}
