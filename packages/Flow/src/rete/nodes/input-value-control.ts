/** The inline value editor on an input node's body — a Rete control whose
 *  React side (`DgControlComponent` in node-component.tsx) mounts a real
 *  Datagrok input built by `buildInputValueEditor`.
 *
 *  The editor is built lazily on first render and cached for the node's
 *  lifetime: React re-renders (status changes, invalidation) re-attach the
 *  same element, so focus and caret survive re-renders mid-typing. */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {buildInputValueEditor, InputValueEditor} from '../../utils/input-values';

export class InputValueControl extends ClassicPreset.Control {
  private editor: InputValueEditor | null | undefined;

  constructor(readonly node: FlowNode) {
    super();
  }

  /** The cached DG editor root, or null when the type has no inline value. */
  element(): HTMLElement | null {
    if (this.editor === undefined) {
      this.editor = buildInputValueEditor(this.node,
        () => this.node.editorBridge?.notifyParamsChanged(this.node.id));
    }
    return this.editor?.root ?? null;
  }

  /** Reflect a store change made elsewhere (the context panel) into the DG
   *  input — programmatic, never reported as an edit. */
  sync(): void {
    this.editor?.sync();
  }
}
