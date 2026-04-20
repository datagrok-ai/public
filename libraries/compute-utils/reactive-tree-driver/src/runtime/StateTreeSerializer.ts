import {Observable, BehaviorSubject, merge} from 'rxjs';
import {map, scan, debounceTime} from 'rxjs/operators';
import {BaseTree, TreeNode} from '../data/BaseTree';
import {PipelineSerializedState, PipelineState} from '../config/PipelineInstance';
import {ConsistencyInfo, FuncCallNode, FuncCallStateInfo, isFuncCallNode, StateTreeNode, StateTreeSerializationOptions} from './StateTreeNodes';
import {LinksState} from './LinksState';
import {ValidationResult} from '../data/common-types';

export function toStateRec(
  node: TreeNode<StateTreeNode>,
  options: StateTreeSerializationOptions = {},
  linksState?: LinksState,
): PipelineState {
  const item = node.getItem();
  const actions = linksState ? linksState.getNodeActionsData(item.uuid) : undefined;
  if (isFuncCallNode(item))
    return item.toState(options, actions);
  const selfState = item.toState(options, actions);
  const steps = node.getChildren().map((node) => {
    const item = toStateRec(node.item, options, linksState);
    return item;
  });
  const fullState = {...selfState, steps};
  const structureCheckResults = item.getStructureCheck(fullState);
  return {...fullState, structureCheckResults};
}

export function toSerializedStateRec(
  node: TreeNode<StateTreeNode>,
  options: StateTreeSerializationOptions = {},
): PipelineSerializedState {
  const item = node.getItem();
  if (isFuncCallNode(item))
    return item.toSerializedState(options);
  const state = item.toSerializedState(options);
  const steps = node.getChildren().map((node) => {
    const item = toSerializedStateRec(node.item, options);
    return item;
  });
  return {...state, steps};
}

export function getValidations(nodeTree: BaseTree<StateTreeNode>) {
  const entries = nodeTree.traverse(nodeTree.root, (acc, node) => {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return [...acc, [item.uuid, item.validationInfo$] as const];

    return acc;
  }, [] as (readonly [string, BehaviorSubject<Record<string, ValidationResult>>])[]);
  return Object.fromEntries(entries);
}

export function getConsistency(nodeTree: BaseTree<StateTreeNode>) {
  const entries = nodeTree.traverse(nodeTree.root, (acc, node) => {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return [...acc, [item.uuid, item.consistencyInfo$] as const];

    return acc;
  }, [] as (readonly [string, BehaviorSubject<Record<string, ConsistencyInfo>>])[]);
  return Object.fromEntries(entries);
}

export function getMeta(nodeTree: BaseTree<StateTreeNode>) {
  const entries = nodeTree.traverse(nodeTree.root, (acc, node) => {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return [...acc, [item.uuid, item.metaInfo$] as const];
    return acc;
  }, [] as (readonly [string, BehaviorSubject<Record<string, BehaviorSubject<any | undefined>>>])[]);
  return Object.fromEntries(entries);
}

export function getFuncCallStates(nodeTree: BaseTree<StateTreeNode>) {
  const entries = nodeTree.traverse(nodeTree.root, (acc, node) => {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return [...acc, [item.uuid, item.funcCallState$] as const];
    return acc;
  }, [] as (readonly [string, BehaviorSubject<FuncCallStateInfo | undefined>])[]);
  return Object.fromEntries(entries);
}

export function getNodesDescriptions(nodeTree: BaseTree<StateTreeNode>) {
  const entries = nodeTree.traverse(nodeTree.root, (acc, node) => {
    const item = node.getItem();
    const stateNames = item.nodeDescription.getStateNames();
    const stateChanges = stateNames.map((name) => item.nodeDescription.getStateChanges(name).pipe(
      map((val) => [name, val] as const),
    ));
    const descriptions$ = merge(...stateChanges).pipe(
      scan((acc, [name, val]) => {
        if (name === 'tags') {
          const tags = Object.values(val ?? {}).flat().filter((x) => x) as string[];
          return {...acc, [name]: tags};
        }
        return {...acc, [name]: val};
      }, {} as Record<string, string[] | string>),
      debounceTime(0),
    );
    return [...acc, [item.uuid, descriptions$] as const];
  }, [] as (readonly [string, Observable<Record<string, string | string[]> | undefined>])[]);
  return Object.fromEntries(entries);
}

export function getIOMutations(nodeTree: BaseTree<StateTreeNode>) {
  const allFlags = nodeTree.traverse(nodeTree.root, (acc, node) => {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return [...acc, item.instancesWrapper.getIOEditsFlag()];
    return acc;
  }, [] as Observable<boolean>[]);
  return merge(...allFlags);
}
