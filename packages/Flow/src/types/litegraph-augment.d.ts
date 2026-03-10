/** Type augmentations for litegraph.js properties missing from official .d.ts */

import 'litegraph.js';

declare module 'litegraph.js' {
  interface LGraphNode {
    /** Array of widgets attached to the node (runtime property, missing from .d.ts) */
    widgets?: IWidget[];
  }

  interface LGraph {
    /** Callback when a node is removed from the graph (runtime, missing from .d.ts) */
    onNodeRemoved?: (node: LGraphNode) => void;
    /** Callback when a link is added */
    onLinkAdded?: (link: LLink) => void;
    /** Callback when a link is removed */
    onLinkRemoved?: (link: LLink) => void;
  }
}
