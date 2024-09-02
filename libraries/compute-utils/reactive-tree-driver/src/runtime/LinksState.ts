import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {NodePath} from '../data/BaseTree';
import {StateTree} from './StateTree';
import {isSequentialPipelineNode, isStaticPipelineNode} from './StateTreeNodes';
import {matchNodeLink} from './link-matching';
import {Link} from './Link';


export class LinksState {
  constructor() {}

  public update(state: StateTree, mutationPath: NodePath) {
    const nlinks = this.createLinks(state);
    // TODO
  }

  public createLinks(state: StateTree) {
    const links = state.traverse(state.getRoot(), (acc, node, path) => {
      const item = node.getItem();
      if (isStaticPipelineNode(item) || isSequentialPipelineNode(item)) {
        const {config} = item;
        const matchedLinks = (config.links ?? [])
          .map((link) => matchNodeLink(node, link))
          .filter((x) => !!x)
          .flat();
        const links = matchedLinks.map((minfo) => new Link(path, minfo));
        return [...acc, ...links];
      }
      return acc;
    }, [] as Link[]);
    return links;
  }
}
