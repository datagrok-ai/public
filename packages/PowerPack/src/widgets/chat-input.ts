/* The platform chat's compose box: the u2 MessageInput with @-mentions of users. The Dart chat talks
   to it through DOM events on the root, with the markup in `data-markup`: the box fires `u2-msg-send`
   and empties once the chat answers `u2-msg-sent` (`u2-msg-failed` keeps the text); the chat fires
   `u2-msg-insert` to drop markup in (an entity dragged onto the box) and `u2-msg-submit` to send what
   the box holds, which an empty box answers with an empty `u2-msg-send`. */
import * as DG from 'datagrok-api/dg';
import {messageInput} from '@datagrok-libraries/u2/src/dg/index.js';
import '@datagrok-libraries/u2/css/message-input.css';

export function chatInput(value?: string, placeholder?: string, sendOnEnter?: boolean, draftKey?: string): DG.Widget {
  const post = (markup: string) => new Promise<void>((resolve, reject) => {
    const root = input.root;
    const settle = (e: Event) => {
      root.removeEventListener('u2-msg-sent', settle);
      root.removeEventListener('u2-msg-failed', settle);
      if (e.type === 'u2-msg-sent')
        resolve();
      else
        reject(new Error('the chat did not post the message'));
    };
    root.addEventListener('u2-msg-sent', settle);
    root.addEventListener('u2-msg-failed', settle);
    root.dataset.markup = markup;
    root.dispatchEvent(new Event('u2-msg-send'));
  });
  const input = messageInput({
    value: value ?? undefined,
    placeholder: placeholder ?? undefined,
    sendOn: sendOnEnter ? 'enter' : 'ctrlEnter',
    draftKey: draftKey ?? undefined,
    onSend: post,
  });
  input.root.addEventListener('u2-msg-insert', () => {
    input.value = `${input.value}${input.root.dataset.markup ?? ''}`;
    input.focus();
  });
  input.root.addEventListener('u2-msg-submit', () => {
    if (input.empty.peek())
      post('').catch(() => undefined);
    else
      input.send();
  });
  return new DG.Widget(input.root).own(() => input.dispose());
}
