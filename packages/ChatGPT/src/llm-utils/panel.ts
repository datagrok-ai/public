/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {OpenAI} from 'openai';
// @ts-ignore .... idk why it does not like it
import '../../css/ai.css';
import {ChatModel} from 'openai/resources/shared';
import {dartLike, fireAIAbortEvent} from '../utils';
import {ConversationStorage, StoredConversationWithContext, UIMessage} from './storage';

export type ModelOption = 'Fast' | 'Deep Research';
export const ModelType: {[type in ModelOption]: ChatModel} = {
  Fast: 'gpt-4o-mini',
  //'Deep Research': 'gpt-5.2' as ChatModel, // hell of a smart model but a bit expensive expensive
  'Deep Research': 'o4-mini', // good balance between speed, quality and $$$
} as const;

// in future might extend it with other types for response API
type MessageType = OpenAI.Chat.ChatCompletionMessageParam;

type AIPanelInputs = {
    prompt: string,
    model: ModelOption,
}

type DBAIPanelInputs = AIPanelInputs & {
    schemaName: string,
}

const actionButtionValues = {
  run: 'Run AI Prompt',
  stop: 'Stop AI Generation',
} as const;

export type UIMessageOptions = {
  result?: {
    finalResult?: string;
  }
}

export class AIPanel<T extends MessageType = OpenAI.Chat.ChatCompletionMessageParam, K extends AIPanelInputs = AIPanelInputs> {
  private root: HTMLElement;
  private inputArea: HTMLElement;
  protected header: HTMLElement;
  private outputArea: HTMLElement;
  private textArea: HTMLTextAreaElement;
  private runButton: HTMLElement;
  private modelInput: DG.InputBase<ModelOption>;
  private newChatButton: HTMLElement;
  private copyConversationButton: HTMLElement;
  private historyButton: HTMLElement;
  private tryAgainButton: HTMLElement;
  private _onRunRequest = new rxjs.Subject<{prevMessages: T[], currentPrompt: K}>();
  private _messages: T[] = [];
  private _uiMessages: UIMessage[] = [];
  protected get placeHolder() { return 'Type your prompt here...'; }
  public get onRunRequest() {
    return this._onRunRequest.asObservable();
  }
  private _onClearChatRequest: rxjs.Subject<void> = new rxjs.Subject<void>();
  public get onClearChatRequest(): rxjs.Observable<void> {
    return this._onClearChatRequest.asObservable();
  }
  private runButtonTooltip: typeof actionButtionValues[keyof typeof actionButtionValues] = actionButtionValues.run;
  public inputControlsDiv: HTMLElement;
  protected get contextId(): string {
    return this._contextID;// these should be overriden in subclasses
  }


  private currentConversationId: string | null = null;
  constructor(private _contextID: string = 'global-ai-panel') {
    this.root = ui.divV([], 'd4-ai-generation-panel');
    this.inputArea = ui.divV([], 'd4-ai-panel-input-area');
    this.outputArea = ui.divV([], 'd4-ai-panel-output-area');
    this.header = ui.divH([], 'd4-ai-panel-header');
    this.textArea = dartLike(document.createElement('textarea'))
      .set('className', 'd4-ai-input-textarea').set('placeholder', this.placeHolder).value;
    this.runButton = ui.iconFA('paper-plane', () => {
      if (this.runButtonTooltip === actionButtionValues.stop) {
        this.terminate();
        return;
      }
      this.handleRun();
    });
    this.textArea.addEventListener('keydown', (event: KeyboardEvent) => {
      if (event.key === 'Enter' && (!event.ctrlKey && !event.metaKey)) {
        event.preventDefault();
        event.stopImmediatePropagation();
        this.handleRun();
      }
    });
    ui.tooltip.bind(this.runButton, () => this.runButtonTooltip, 'left');
    this.tryAgainButton = ui.icons.sync(() => this.tryAgain(), 'Try Again');
    this.historyButton = ui.iconFA('history', () => this.showHistory(), 'Chat History...');
    this.copyConversationButton = ui.iconFA('copy', async () => {
      const success = await this.copyConversationToClipboard();
      if (success)
        grok.shell.info('Conversation copied to clipboard');
      else
        grok.shell.error('Failed to copy conversation to clipboard');
    }, 'Copy Conversation to Clipboard');

    this.newChatButton = ui.icons.add(async () => {
      ui.setUpdateIndicator(this.root, true, 'Saving conversation...');
      try {
        await this.saveCurrentConversation(); // will always be first user message
      } catch (error) {
        grok.shell.error('Failed to save conversation before clearing.');
        console.error('Failed to save conversation before clearing:', error);
      } finally {
        ui.setUpdateIndicator(this.root, false);
      }
      this.handleClear();
      this.currentConversationId = null;
    }, 'Start New Chat');
    this.modelInput = ui.input.choice('Model', {
      items: Object.keys(ModelType) as ModelOption[],
      value: 'Deep Research',
      nullable: false,
      tooltipText: 'Select AI model to use',
    }) as DG.InputBase<ModelOption>;
    this.hideContentIcons();
    this.inputControlsDiv = ui.divH([
      this.historyButton,
      this.modelInput.input
    ], 'd4-ai-panel-input-controls');
    this.runButton.style.color = 'var(--blue-1)';
    const inputControlsRight = ui.divH([this.tryAgainButton, this.runButton], 'd4-ai-panel-run-controls');
    inputControlsRight.style.marginLeft = 'auto';
    const controlsDiv = ui.divH([this.inputControlsDiv, inputControlsRight], 'd4-ai-panel-controls-container');
    this.inputArea.appendChild(this.textArea);
    this.inputArea.appendChild(controlsDiv);
    this.root.appendChild(this.header);
    this.root.appendChild(this.outputArea);
    this.root.appendChild(this.inputArea);

    const headerTitle = ui.h2('Chat', 'd4-ai-panel-header-title');
    dartLike(headerTitle.style).set('margin', '0px').set('userSelect', 'none').set('cursor', 'pointer');
    headerTitle.addEventListener('click', () => this.textArea.focus());
    this.header.appendChild(headerTitle);
    const rightHeaderDiv = ui.divH([], 'd4-ai-panel-header-right-buttons');
    rightHeaderDiv.appendChild(this.newChatButton);
    rightHeaderDiv.appendChild(this.copyConversationButton);
    rightHeaderDiv.appendChild(this.historyButton);
    this.header.appendChild(rightHeaderDiv);
  }

  private static _lastDockedPanel: DG.DockNode | null = null;
  show() {
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showContextPanel = false;
    // grok.shell.windows.context.root?.style
    // grok.shell.o = this.root;
    // grok.shell.setCurrentObject(this.root, true, true);
    if (AIPanel._lastDockedPanel) {
      if (AIPanel._lastDockedPanel.container.containerElement.contains(this.root) && document.contains(this.root))
        return;
      grok.shell.dockManager.close(AIPanel._lastDockedPanel);
    }
    AIPanel._lastDockedPanel = grok.shell.dockManager.dock(this.root, 'right', null, undefined, 0.25);
  }

  hide() {
    // save history before hiding
    grok.shell.dockManager.close(this.root);
  }

  formatConversation() {
    let previousFromAI = false;
    return this._uiMessages.map((msg) => {
      const prefix = msg.fromUser ? 'USER:' : 'Datagrok AI:';
      const separator = '-'.repeat(50);
      const isContinuation = !msg.fromUser && previousFromAI;
      previousFromAI = !msg.fromUser;
      const actualPrefix = isContinuation ? '\n' : `${separator}\n${prefix}\n`;
      return `${actualPrefix}${msg.text}\n`;
    }).join('\n');
  }

  async copyConversationToClipboard() {
    const formattedText = this.formatConversation();

    try {
      await navigator.clipboard.writeText(formattedText);
      return true;
    } catch (err) {
      console.error('Failed to copy:', err);
      return false;
    }
  }

  toggle() {
    document.contains(this.root) ? this.hide() : this.show();
  }

  dispose() {
    this.root.remove();
    this._messages = [];
    this._uiMessages = [];
  }


  public getCurrentInputs(): K {
    return {
      prompt: this.textArea.value,
      model: this.modelInput.value!,
    } as K;
  }

  private _aiMessagesAccordionPane: HTMLElement | null = null;
  protected appendMessage(aiMessage: T, uiMessage: {title: string, content: string, fromUser: boolean, onlyAddToMessages?: boolean, uiOnly?: boolean, finalMessage?: string}, loader?: HTMLElement) {
    if (!uiMessage.uiOnly)
      this._messages.push(aiMessage);
    if (uiMessage.onlyAddToMessages)
      return;
    // from this point we know that message is also in the ui.
    this._uiMessages.push({fromUser: !!uiMessage.fromUser, text: uiMessage.content, title: uiMessage.title});
    if (uiMessage.fromUser) {
      const userDiv = ui.div(ui.divText(uiMessage.content, 'd4-ai-user-prompt-divtext'), 'd4-ai-user-prompt-container');
      this.outputArea.appendChild(userDiv);
      this._aiMessagesAccordionPane = null; // reset accordion pane so that next AI message creates a new one
    } else {
      if (!this._aiMessagesAccordionPane) {
        const acord = ui.accordion();
        this._aiMessagesAccordionPane = ui.divV([], 'd4-ai-messages-accordion-pane');
        const pane = acord.addPane('Responses', () => this._aiMessagesAccordionPane!, true, undefined, false);
        pane.expanded = true;
        acord.root.style.width = 'calc(100% - 35px)';
        this.outputArea.appendChild(acord.root);
      }
      const markDown = ui.markdown(uiMessage.content);
      // if there is code block, make it copyable
      markDown.style.position = 'relative';
      if (markDown.querySelector('pre > code')) {
        const copyButton = ui.icons.copy(() => {}, 'Copy Code');
        //dartLike(copyButton.style).set('position', 'absolute').set('top', '5px').set('right', '5px').set('width', '20px').set('height', '20px').set('opacity', '0.6').set('cursor', 'pointer');
        copyButton.classList.add('d4-ai-copy-code-button');
        markDown.appendChild(copyButton);
        copyButton.addEventListener('click', () => {
          const codeElement = markDown.querySelector('pre > code');
          if (codeElement) {
            // make sure the header of the code is not overlapped with the button
            const header = markDown.children[0];
            if (header && header.tagName?.toLowerCase() !== 'pre')
              (header as HTMLElement).style.marginRight = '16px';

            navigator.clipboard.writeText(codeElement.textContent || '').then(() => {
              copyButton.classList.add('d4-ai-copy-code-button-copied');
              setTimeout(() => {
                copyButton.classList.remove('d4-ai-copy-code-button-copied');
              }, 600);
            }).catch((err) => {
              grok.shell.error('Failed to copy code to clipboard.');
              console.error('Failed to copy code to clipboard:', err);
            });
          }
        });
      }

      if (uiMessage.finalMessage && !uiMessage.fromUser) {
        // if it is a final message from AI, add thumbs up/down,
        const feedbackDiv = ui.divH([], 'd4-ai-panel-feedback-div');
        const thumbsUp = ui.iconFA('thumbs-up', () => {
          grok.shell.info('Thanks for your feedback!');
          handleFeedback(true);
        }, 'Helpful');
        const thumbsDown = ui.iconFA('thumbs-down', () => {
          grok.shell.info('Thanks for your feedback!');
          handleFeedback(false);
        }, 'Not Helpful');
        [thumbsUp, thumbsDown].forEach((el) => {
          dartLike(el.style).set('padding', '2px').set('borderRadius', '6px');
        });

        const that = this;
        function handleFeedback(helpful: boolean) {
          [thumbsUp, thumbsDown].forEach((el) => {
            el.style.backgroundColor = '';
          });
          (helpful ? thumbsUp : thumbsDown).style.backgroundColor = helpful ? 'rgba(0, 150, 30, 0.2)' : 'rgba(200, 0, 0, 0.2)';
          receiveFeedback(
            that._uiMessages[0]?.text ?? '',
            uiMessage.finalMessage!,
            helpful
          );
        }
        feedbackDiv.appendChild(thumbsUp);
        feedbackDiv.appendChild(thumbsDown);
        dartLike(feedbackDiv.style).set('gap', '8px').set('alignItems', 'center').set('width', '100%').set('paddingBottom', '8px').set('paddingLeft', '4px');
        markDown.appendChild(feedbackDiv);
      }
      dartLike(markDown.style).set('userSelect', 'text').set('maxWidth', '100%');
      const titleText = uiMessage.title ? [ui.h3(uiMessage.title, 'd4-ai-assistant-response-title')] : [];
      const assistantDiv = ui.divV([...titleText, markDown], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane.appendChild(assistantDiv);
    }
    // reapend loader to the end
    if (loader)
      this.outputArea.appendChild(loader);
    this.outputArea.scrollTop = this.outputArea.scrollHeight;
    this.showContentIcons();
  }

  public startChatSession() {
    const loader = ui.icons.loader();
    dartLike(loader.style).set('alignSelf', 'center').set('height', '20px').set('marginTop', '8px');
    this.runButton.classList.remove('fal', 'fa-paper-plane');
    this.runButton.classList.add('fas', 'fa-stop');
    this.runButton.style.color = 'orangered';
    this.runButtonTooltip = actionButtionValues.stop;
    return {
      addUserMessage: (message: T, content: string) => {
        this.appendMessage(message, {title: '', content: content, fromUser: true}, loader);
      },
      addAIMessage: (message: T, title: string, content: string, onlyAddToMessages?: boolean) => {
        this.appendMessage(message, {title: title, content: content, fromUser: false, onlyAddToMessages}, loader);
      },
      addUIMessage: (content: string, fromUser: boolean, messageOptions?: UIMessageOptions) => {
        this.appendMessage('' as any, {title: '', content: content, fromUser: fromUser, uiOnly: true, finalMessage: messageOptions?.result?.finalResult}, loader);
      },
      endSession: () => {
        this.runButton.classList.remove('fas', 'fa-stop');
        this.runButton.classList.add('fal', 'fa-paper-plane');
        this.runButton.style.color = 'var(--blue-1)';
        this.runButtonTooltip = actionButtionValues.run;
        this.saveCurrentConversation().catch((e) => console.error('Failed to save conversation before hiding panel:', e));
        loader.remove();
      }
    };
  }

  get isShown(): boolean {
    return document.contains(this.root);
  }

  private tryAgain() {
    if (this._messages.length === 0)
      return; // should never happen, but just in case
    const inputs = this.getCurrentInputs();
    inputs.prompt = 'Please try again. ';
    this._onRunRequest.next({
      prevMessages: this._messages,
      currentPrompt: inputs,
    });
  }

  private terminate() {
    fireAIAbortEvent();
  }

  private handleRun() {
    const inputs = this.getCurrentInputs();
    this.textArea.value = '';
    this._onRunRequest.next({
      prevMessages: this._messages,
      currentPrompt: inputs,
    });
  }

  protected showContentIcons() {
    this.newChatButton.style.display = 'flex';
    this.tryAgainButton.style.display = 'flex';
    this.copyConversationButton.style.display = 'flex';
  }

  protected hideContentIcons() {
    this.newChatButton.style.display = 'none';
    this.tryAgainButton.style.display = 'none';
    this.copyConversationButton.style.display = 'none';
  }

  private handleClear() {
    this._messages = [];
    this._uiMessages = [];
    this.outputArea.innerHTML = '';
    this._onClearChatRequest.next();
    this.hideContentIcons();
  }

  private async showHistory() {
    try {
      const conversations = await ConversationStorage.listConversations(
        this.contextId,
        10
      ) as StoredConversationWithContext<T>[];

      if (conversations.length === 0) {
        grok.shell.info('No conversation history found');
        return;
      }

      const menu = DG.Menu.popup();
      for (const conv of conversations) {
        if (this.currentConversationId && conv.id === this.currentConversationId)
          continue;
        const date = new Date(conv.timestamp).toLocaleString();
        const preview = conv.initialPrompt.substring(0, 40);
        const displayText = `${preview}${conv.initialPrompt.length > 40 ? '...' : ''} (${date})`;

        menu.item(displayText, () => this.loadConversation(conv.id));
      }

      menu.separator();
      menu.item('Clear All History', async () => {
        // if (await grok.shell.('Delete all conversation history?')) {
        //   await ConversationStorage.clearAll();
        //   grok.shell.info('History cleared');
        // }
        ui.dialog('Delete all conversation history').add(ui.divText('This action will permanently delete all saved conversations. Are you sure you want to proceed?'))
          .onOK(async () => {
            await ConversationStorage.clearAll();
            grok.shell.info('History cleared');
          }).show();
      });

      menu.show();
    } catch (error) {
      console.error('Failed to load history:', error);
      grok.shell.error('Failed to load conversation history');
    }
  }
  private async loadConversation(conversationId: string) {
    try {
      const conv = (await ConversationStorage.getConversation(conversationId)) as StoredConversationWithContext<T> | null;
      if (!conv) {
        grok.shell.error('Conversation not found');
        return;
      }

      this.currentConversationId = conversationId;
      this._messages = conv.messages;

      // Clear and rebuild UI from messages
      this.outputArea.innerHTML = '';
      this._uiMessages = [];
      conv.uiMessages.forEach((msg) => {
        this.appendMessage(null as any, {title: msg.title ?? '', content: msg.text, fromUser: msg.fromUser, uiOnly: true}); // no loader
      });
      //grok.shell.info(`Loaded conversation: ${conv.initialPrompt.substring(0, 50)}...`);
    } catch (error) {
      console.error('Failed to load conversation:', error);
      grok.shell.error('Failed to load conversation');
    }
  }

  public async saveCurrentConversation(initialPrompt?: string) {
    if (this._messages.length === 0) return;

    try {
      if (this.currentConversationId)
        await ConversationStorage.updateConversation(this.currentConversationId, this._messages, this._uiMessages);
      else {
        this.currentConversationId = await ConversationStorage.saveConversation(
          this._messages, this._uiMessages, initialPrompt ?? this._uiMessages[0]?.text ?? 'No Title', this.contextId
        );
      }
    } catch (error) {
      console.error('Failed to save conversation:', error);
    }
  }
}

export class DBAIPanel<T extends MessageType = OpenAI.Chat.ChatCompletionMessageParam> extends AIPanel<T, DBAIPanelInputs> {
  protected get placeHolder() { return 'Ask your database, like "Total sales by regions"'; }
  protected schemaInput: DG.InputBase<string>;
  constructor(schemas: string[], defaultSchema: string, connectionID: string) {
    super(connectionID); // context ID is connection ID
    this.schemaInput = ui.input.choice('Schema', {
      items: schemas,
      value: defaultSchema,
      nullable: false,
      tooltipText: 'Select the database schema to use for AI-assisted query generation.',
    }) as DG.InputBase<string>;
    this.inputControlsDiv.appendChild(this.schemaInput.input);
    ui.tooltip.bind(this.schemaInput.input, 'Select the database schema to use for AI-assisted query generation.');
  }

  public getCurrentInputs(): DBAIPanelInputs {
    const baseInputs = super.getCurrentInputs();
    return {
      ...baseInputs,
      schemaName: this.schemaInput.value!,
    };
  }
}


function receiveFeedback(userPrompt: string, aiResponse: string, helpful: boolean) {
  // not implemented yet
}
