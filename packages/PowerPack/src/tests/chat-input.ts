import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {chatInput} from '../widgets/chat-input';

// The event protocol the platform chat (xamgle chat.dart, ChatInput) speaks to the compose box.
category('Chat input', () => {
  function sends(root: HTMLElement, answer: 'u2-msg-sent' | 'u2-msg-failed'): string[] {
    const sent: string[] = [];
    root.addEventListener('u2-msg-send', () => {
      sent.push(root.dataset.markup!);
      root.dispatchEvent(new Event(answer));
    });
    return sent;
  }

  function text(root: HTMLElement): string {
    return root.querySelector('.u2-msg-editor')!.textContent!;
  }

  test('is the chatInput func', async () => {
    const funcs = DG.Func.find({package: 'PowerPack', tags: ['chatInput']});
    expect(funcs.map((f) => f.name).join(), 'chatMessageInput');
  });

  test('submit sends the markup and empties the box once sent', async () => {
    const w = chatInput('hello\nthere');
    const sent = sends(w.root, 'u2-msg-sent');
    w.root.dispatchEvent(new Event('u2-msg-submit'));
    await DG.delay(0);
    expect(sent.join('|'), 'hello\nthere');
    expect(text(w.root), '');
    w.detach();
  });

  test('a failed send keeps the text', async () => {
    const w = chatInput('keep me');
    sends(w.root, 'u2-msg-failed');
    w.root.dispatchEvent(new Event('u2-msg-submit'));
    await DG.delay(0);
    expect(text(w.root), 'keep me');
    w.detach();
  });

  test('an empty box answers submit with an empty send', async () => {
    const w = chatInput();
    const sent = sends(w.root, 'u2-msg-failed');
    w.root.dispatchEvent(new Event('u2-msg-submit'));
    await DG.delay(0);
    expect(sent.length, 1);
    expect(sent[0], '');
    w.detach();
  });

  test('insert appends markup', async () => {
    const w = chatInput('see ');
    const sent = sends(w.root, 'u2-msg-sent');
    w.root.dataset.markup = 'the table';
    w.root.dispatchEvent(new Event('u2-msg-insert'));
    w.root.dispatchEvent(new Event('u2-msg-submit'));
    await DG.delay(0);
    expect(sent[0], 'see the table');
    w.detach();
  });
});
