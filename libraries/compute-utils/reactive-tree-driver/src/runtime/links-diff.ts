import {Link} from './Link';
import isEqual from 'lodash.isequal';

export interface LinksDiff {
  toRemove: Set<string>,
  toAdd: Set<string>,
}

// TODO: probably n*log(n) alg someday
export function getLinksDiff(oldLinks: Link[], newLinks: Link[]): LinksDiff {
  const toKeep = new Set<string>();
  const toRemove = new Set<string>();
  const toAdd = new Set<string>();
  const notToAdd = new Set<string>();

  for (const oldLink of oldLinks) {
    const foundNew = findLink(newLinks, oldLink);
    if (foundNew) {
      toKeep.add(oldLink.uuid);
      notToAdd.add(foundNew.uuid);
    } else
      toRemove.add(oldLink.uuid);
  }

  for (const newLink of newLinks) {
    if (!notToAdd.has(newLink.uuid))
      toAdd.add(newLink.uuid);
  }

  return {
    toRemove,
    toAdd,
  };
}

export function findLink(links: Link[], link: Link) {
  return links.find((l) => linkEqual(link, l));
}

export function linkEqual(link1: Link, link2: Link): boolean {
  if (!isEqual(link1.matchInfo.spec.id, link2.matchInfo.spec.id))
    return false;

  if (!isEqual(link1.prefix, link2.prefix))
    return false;

  if (!isEqual(link1.matchInfo.basePath, link2.matchInfo.basePath))
    return false;

  if (!isEqual(link1.matchInfo.basePathUUID, link2.matchInfo.basePathUUID))
    return false;

  if (!isEqual(link1.matchInfo.inputs, link2.matchInfo.inputs))
    return false;

  if (!isEqual(link1.matchInfo.inputsUUID, link2.matchInfo.inputsUUID))
    return false;

  if (!isEqual(link1.matchInfo.outputs, link2.matchInfo.outputs))
    return false;

  if (!isEqual(link1.matchInfo.outputsUUID, link2.matchInfo.outputsUUID))
    return false;

  return true;
}
