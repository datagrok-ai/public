import {BaseTree} from '../data/BaseTree';
import {Link} from './Link';
import {StateTreeNode} from './StateTreeNodes';
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

  const filteredOldLinks: Link[] = [];
  for (const link of oldLinks) {
    if (!link.matchInfo.isDefaultValidator)
      filteredOldLinks.push(link);
    else
      toRemove.add(link.uuid);
  }

  for (const oldLink of filteredOldLinks) {
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

  if (!isEqual(link1.matchInfo.inputs, link2.matchInfo.inputs))
    return false;

  if (!isEqual(link1.matchInfo.outputs, link2.matchInfo.outputs))
    return false;

  return true;
}