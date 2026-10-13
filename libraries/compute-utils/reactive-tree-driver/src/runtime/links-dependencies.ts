import {BaseTree, NodePathSegment} from '../data/BaseTree';
import {StateTreeNode} from './StateTreeNodes';
import {MatchedIO} from './link-matching';
import {DriverLogger, reportError} from '../data/Logger';
import {Link} from './Link';

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

const isUnlinked = (link: Link, alias: string) => !!link.matchInfo.spec.to?.find((io) => io.name === alias)?.unlinked;

const targetKey = (state: BaseTree<StateTreeNode>, link: Link, info: MatchedIO) =>
  `${state.getNode([...link.prefix, ...info.path]).getItem().uuid}/${info.ioName}`;

// ios written by data links of the config, keyed by node uuid and io name; `$linked` outputs yield to them
function dataLinkTargets(state: BaseTree<StateTreeNode>, links: Link[]) {
  const targets = new Set<string>();
  for (const link of links) {
    if ((link.matchInfo.spec.type ?? 'data') !== 'data')
      continue;
    for (const [alias, infos] of Object.entries(link.matchInfo.outputs)) {
      if (isUnlinked(link, alias))
        continue;
      for (const info of infos)
        targets.add(targetKey(state, link, info));
    }
  }
  return targets;
}

/** Drops from `$linked` outputs the ios that other data links write. */
export function pruneLinkedTargets(state: BaseTree<StateTreeNode>, links: Link[]) {
  const written = dataLinkTargets(state, links);
  for (const link of links) {
    const {outputs, outputsUUID} = link.matchInfo;
    for (const [alias, infos] of Object.entries(outputs)) {
      if (!isUnlinked(link, alias))
        continue;
      const kept = infos.map((info) => !written.has(targetKey(state, link, info)));
      if (kept.every((keep) => keep))
        continue;
      if (kept.some((keep) => keep)) {
        outputs[alias] = infos.filter((_, idx) => kept[idx]);
        outputsUUID.set(alias, (outputsUUID.get(alias) ?? []).filter((_, idx) => kept[idx]));
      } else {
        delete outputs[alias];
        outputsUUID.delete(alias);
      }
    }
  }
}
