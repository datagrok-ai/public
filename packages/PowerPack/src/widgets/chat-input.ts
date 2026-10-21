/* The platform chat's compose box: the u2 MessageInput with @-mentions of users. The Dart chat
   (xamgle chat.dart, ChatInput, through jsu.JsChatInput) sets `onSend` and calls the methods below. */
import * as DG from 'datagrok-api/dg';
import {messageInput} from '@datagrok-libraries/u2/src/dg/index.js';
import type {MessageInput} from '@datagrok-libraries/u2';
import '@datagrok-libraries/u2/css/message-input.css';

export class ChatInputWidget extends DG.Widget {
  /** Posts the markup; false keeps it in the box. */
  onSend: (markup: string) => Promise<boolean> = async () => false;
  private readonly input: MessageInput;

  constructor(value?: string, placeholder?: string, sendOnEnter?: boolean, draftKey?: string) {
    super(document.createElement('div'));
    this.input = messageInput({
      value: value ?? undefined,
      placeholder: placeholder ?? undefined,
      sendOn: sendOnEnter ? 'enter' : 'ctrlEnter',
      draftKey: draftKey ?? undefined,
      onSend: async (markup) => {
        if (!await this.onSend(markup))
          throw new Error('the chat did not post the message');
      },
    });
    this.root.append(this.input.root);
    this.own(() => this.input.dispose());
  }

  /** Sends what the box holds; true once the chat has taken it. */
  async submit(): Promise<boolean> {
    if (this.input.empty.peek() || !await this.onSend(this.input.value))
      return false;
    this.input.clear();
    return true;
  }

  insert(markup: string): void {
    this.input.value = `${this.input.value}${markup}`;
    this.input.focus();
  }

  focus(): void {
    this.input.focus();
  }
}
