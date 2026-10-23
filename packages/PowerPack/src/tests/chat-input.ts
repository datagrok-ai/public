import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {ChatInputWidget} from '../widgets/chat-input';

// What the platform chat (xamgle chat.dart, ChatInput) calls on the compose box.
category('Chat input', () => {
  function text(w: ChatInputWidget): string {
    return w.root.querySelector('.u2-msg-editor')!.textContent!;
  }

  test('is the chatInput func', async () => {
    const funcs = DG.Func.find({package: 'PowerPack', tags: ['chatInput']});
    expect(funcs.map((f) => f.name).join(), 'chatMessageInput');
  });

  test('submit sends the markup and empties the box once taken', async () => {
    const w = new ChatInputWidget('hello\nthere');
    const sent: string[] = [];
    w.onSend = async (markup) => {sent.push(markup); return true;};
    expect(await w.submit(), true);
    expect(sent.join('|'), 'hello\nthere');
    expect(text(w), '');
    w.detach();
  });

  test('a refused send keeps the text', async () => {
    const w = new ChatInputWidget('keep me');
    w.onSend = async () => false;
    expect(await w.submit(), false);
    expect(text(w), 'keep me');
    w.detach();
  });

  test('an empty box sends nothing', async () => {
    const w = new ChatInputWidget();
    let calls = 0;
    w.onSend = async () => {calls++; return true;};
    expect(await w.submit(), false);
    expect(calls, 0);
    w.detach();
  });

  test('insert appends markup', async () => {
    const w = new ChatInputWidget('see ');
    const sent: string[] = [];
    w.onSend = async (markup) => {sent.push(markup); return true;};
    w.insert('the table');
    await w.submit();
    expect(sent[0], 'see the table');
    w.detach();
  });
});
