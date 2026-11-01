/* The platform chat's compose box: the u2 MessageInput with @-mentions of users. The Dart chat
   (xamgle chat.dart, ChatInput, through jsu.JsChatInput) sets `onSend` and calls the methods below. */
import * as DG from 'datagrok-api/dg';
import {messageInput} from '@datagrok-libraries/u2/dg';
import {Scope} from '@datagrok-libraries/u2';
import type {MessageInput} from '@datagrok-libraries/u2';

export class ChatInputWidget extends DG.Widget {
  /** Posts the markup; false keeps it in the box. */
  onSend: (markup: string) => Promise<boolean> = async () => false;
  private readonly input: MessageInput;

  constructor(value?: string, placeholder?: string, sendOnEnter?: boolean, draftKey?: string) {
    super(document.createElement('div'));
    // the Dart chat creates and disposes this widget, so the input lives on a tree of its own
    const owner = Scope.root('chat input');
    this.input = owner.build(() => messageInput({
      value: value ?? undefined,
      placeholder: placeholder ?? undefined,
      sendOn: sendOnEnter ? 'enter' : 'ctrlEnter',
      draftKey: draftKey ?? undefined,
      onSend: async (markup) => {
        if (!await this.onSend(markup))
          throw new Error('the chat did not post the message');
      },
    }));
    this.root.append(this.input.root);
    this.own(() => owner.dispose());
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
