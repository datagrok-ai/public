/** The inline value editor on an input node's body — built lazily and cached for the node's
 *  lifetime so React re-renders re-attach the same element and focus survives mid-typing. */

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
        () => this.node.editorBridge?.notifyParamsChanged(this.node.id), {host: 'node'});
    }
    return this.editor?.root ?? null;
  }

  /** Reflect a store change made elsewhere into the DG input — programmatic, never reported as an edit. */
  sync(): void {
    this.editor?.sync();
  }
}
