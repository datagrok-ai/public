/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
// @ts-ignore .... idk why it does not like it
import '../../css/ai.css';
import {dartLike, fireAIAbortEvent, getAIPanelToggleSubscription} from '../utils';
import {ConversationStorage, StoredConversationWithContext} from './storage';
import {ModelOption, ModelType} from './LLM-client';
import {LanguageModelV3Message, LanguageModelV3Content} from '@ai-sdk/provider';

// in future might extend it with other types for response API
export type MessageType = LanguageModelV3Message;

type AIPanelInputs = {
    prompt: string,
    model: ModelOption,
}

type DBAIPanelInputs = AIPanelInputs & {
    catalogName: string,
}

export type ScriptingAIPanelInputs = AIPanelInputs & {
  language: DG.ScriptingLanguage
};


type TVAIPanelInputs = AIPanelInputs & {
    mode: 'agent' | 'ask',
}

const actionButtionValues = {
  run: 'Run AI Prompt',
  stop: 'Stop AI Generation',
} as const;

export type UIMessageOptions = {
  /** if set, will add feedback buttons to the message*/
  finalResult?: string,
  /** if set, will show a confirmation dialog with the given message before adding the message */
  confirm?: {confirmResult?: boolean, message?: string},
}

export interface UIMessage {
  fromUser: boolean;
  text: string;
  title?: string;
  messageOptions?: UIMessageOptions;
}

export type PanleMesageRet = {
  confirmPromie: Promise<boolean>,
}

export type AIPanelFuncs<T extends MessageType = LanguageModelV3Message> = {
  /** Adds @aiMsg to the message stack (from user) and @msg to the UI panel */
  addUserMessage: (aiMsg: T, msg: string) => void,
  /** Adds @aiMsg to the message stack (from AI) and @msg with @title to the UI panel*/
  addAIMessage: (aiMsg: T, title: string, msg: string) => void,
  /** Adds @aiMsg to the message stack without showing it to the ui  */
  addEngineMessage: (aiMsg: T) => void,
  /** Adds @msg to the UI without adding it to the AI message stack */
  addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => void
  /**Adds confirmation section to the panel and awaits result */
  addConfirmMessage: (msg?: string) => Promise<boolean>,

}

export class AIPanel<T extends MessageType = LanguageModelV3Message, K extends AIPanelInputs = AIPanelInputs> {
  private root: HTMLElement;
  private view: DG.View | DG.ViewBase;
  private inputArea: HTMLElement;
  protected header: HTMLElement;
  private outputArea: HTMLElement;
  private textArea: HTMLTextAreaElement;
  private textAreaDiv: HTMLElement;
  private runButton: HTMLElement;
  private modelInput: DG.InputBase<ModelOption>;
  private newChatButton: HTMLElement;
  private copyConversationButton: HTMLElement;
  private historyButton: HTMLElement;
  private tryAgainButton: HTMLElement;
  private micButton: HTMLElement;
  private recognition: SpeechRecognition | null = null;
  private isRecognizing: boolean = false;
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
  constructor(private _contextID: string = 'global-ai-panel', view: DG.View | DG.ViewBase) {
    this.view = view;
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
    this.micButton = ui.iconFA('microphone', () => this.toggleSpeechRecognition(), 'Voice Input');
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
      this.micButton,
      this.modelInput.input
    ], 'd4-ai-panel-input-controls');
    this.runButton.style.color = 'var(--blue-1)';
    const inputControlsRight = ui.divH([this.tryAgainButton, this.runButton], 'd4-ai-panel-run-controls');
    inputControlsRight.style.marginLeft = 'auto';
    const controlsDiv = ui.divH([this.inputControlsDiv, inputControlsRight], 'd4-ai-panel-controls-container');
    this.textAreaDiv = ui.divV([this.textArea], {classes: 'd4-ai-input-textarea-div', style: {position: 'relative'}});
    this.inputArea.appendChild(this.textAreaDiv);
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
    this.setupSubscriptions();
  }

  protected setupSubscriptions() {
    // do some subscriptions
    let wasShown = false;
    const sub = grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.v != this.view) {
        wasShown = this.isShown;
        if (wasShown)
          this.hide();
      } else {
        if (wasShown)
          this.show();
      }
    });

    const toggleSub = getAIPanelToggleSubscription().subscribe((rv) => {
      if (rv == this.view)
        this.toggle();
    });
    const that = this;
    function onKeyDownHandler(event: KeyboardEvent) {
      if (grok.shell.v === that.view && event.ctrlKey && event.key === 'i') {
        event.preventDefault();
        event.stopImmediatePropagation();
        that.toggle();
      }
    }
    document.addEventListener('keydown', onKeyDownHandler);

    const closeSub = grok.events.onViewRemoved.subscribe((view) => {
      if (view == this.view) {
        sub.unsubscribe();
        closeSub.unsubscribe();
        this.hide();
        this.dispose();
        toggleSub.unsubscribe();
        document.removeEventListener('keydown', onKeyDownHandler);
      }
    });
  }

  private static _lastDockedPanel: DG.DockNode | null = null;
  show() {
    grok.shell.windows.showHelp = false;
    // grok.shell.windows.context.root?.style
    // grok.shell.o = this.root;
    // grok.shell.setCurrentObject(this.root, true, true);
    if (AIPanel._lastDockedPanel) {
      if (AIPanel._lastDockedPanel.container.containerElement.contains(this.root) && document.contains(this.root))
        return;
      grok.shell.dockManager.close(AIPanel._lastDockedPanel);
    }
    const contextRoot = grok.shell.windows.context?.root;
    const contextNode = contextRoot && document.contains(contextRoot) ? grok.shell.dockManager.findNode(contextRoot) : null;
    AIPanel._lastDockedPanel = grok.shell.dockManager.dock(this.root, contextNode ? 'down' : 'right', contextNode, '', contextNode ? 0.5 : 0.25);
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
    this.stopRecognition();
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
  protected appendMessage(
    aiMessage: T, uiMessage: {
      title: string, content: string, fromUser: boolean, onlyAddToMessages?: boolean, uiOnly?: boolean, messageOptions?: UIMessageOptions
    }, loader?: HTMLElement
  ): PanleMesageRet | undefined {
    let ret: PanleMesageRet | undefined = undefined;
    if (!uiMessage.uiOnly)
      this._messages.push(aiMessage);
    if (uiMessage.onlyAddToMessages)
      return;
    // from this point we know that message is also in the ui.
    this._uiMessages.push({fromUser: !!uiMessage.fromUser, text: uiMessage.content, title: uiMessage.title, messageOptions: uiMessage.messageOptions});
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

      if (uiMessage?.messageOptions?.finalResult && !uiMessage.fromUser) {
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
            uiMessage?.messageOptions?.finalResult!,
            that.contextId,
            helpful
          );
        }
        feedbackDiv.appendChild(thumbsUp);
        feedbackDiv.appendChild(thumbsDown);
        dartLike(feedbackDiv.style).set('gap', '8px').set('alignItems', 'center').set('width', '100%').set('paddingBottom', '8px').set('paddingLeft', '4px');
        markDown.appendChild(feedbackDiv);
      }

      // if the action from AI requires confirmation, add confirm section with buttons
      if (uiMessage?.messageOptions?.confirm && !uiMessage.fromUser) {
        let resolve: (value: boolean) => void = () => {};
        const confirmPromise = new Promise<boolean>((res) => { resolve = res; });
        ret = {confirmPromie: confirmPromise};
        const confirmDiv = ui.divV([], 'd4-ai-panel-confirmation-div');
        const confirmText = ui.divText(uiMessage?.messageOptions?.confirm.message ?? 'Allow Datagrok to proceed?', 'd4-ai-panel-confirmation-text');
        const buttonsDiv = ui.divH([], 'd4-ai-panel-confirmation-buttons-div');
        const confirmButton = ui.iconFA('check', null, 'Confirm');
        const cancelButton = ui.iconFA('times', null, 'Cancel');
        dartLike(confirmButton.style).set('padding', '2px').set('borderRadius', '6px').set('marginRight', '2px').set('textAlign', 'center').set('width', '14px');
        dartLike(cancelButton.style).set('padding', '2px').set('borderRadius', '6px').set('textAlign', 'center').set('width', '14px');
        // if the result of confirmation is already set, do not subscribe to changes
        const handleConfirmResult = (confirmed: boolean) => {
          // make sure the reference of the result is also marked
          uiMessage?.messageOptions?.confirm && (uiMessage.messageOptions.confirm.confirmResult = confirmed);
          resolve(confirmed);
          [confirmButton, cancelButton].forEach((el) => {
            el.onclick = null;
            dartLike(el.style).set('backgroundColor', '').set('cursor', 'default').set('opacity', '0.6').set('pointerEvents', 'none');
          });
          (confirmed ? confirmButton : cancelButton).style.backgroundColor = confirmed ? 'rgba(0, 150, 30, 0.3)' : 'rgba(200, 0, 0, 0.3)';
        };
        confirmButton.onclick = () => handleConfirmResult(true);
        cancelButton.onclick = () => handleConfirmResult(false);
        buttonsDiv.appendChild(confirmButton);
        buttonsDiv.appendChild(cancelButton);
        confirmDiv.appendChild(confirmText);
        confirmDiv.appendChild(buttonsDiv);
        markDown.appendChild(confirmDiv);
        if (uiMessage?.messageOptions?.confirm.confirmResult != undefined)
          handleConfirmResult(uiMessage?.messageOptions?.confirm.confirmResult); // already have a result
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
    return ret;
  }

  public startChatSession(): {session: AIPanelFuncs<T>, endSession: () => void} {
    const loader = ui.icons.loader();
    dartLike(loader.style).set('alignSelf', 'center').set('height', '20px').set('marginTop', '8px');
    this.runButton.classList.remove('fal', 'fa-paper-plane');
    this.runButton.classList.add('fas', 'fa-stop');
    this.runButton.style.color = 'orangered';
    this.runButtonTooltip = actionButtionValues.stop;
    return {
      session: {
        addAIMessage: (aiMessage, title, content) => this.appendMessage(aiMessage, {title: title, content: content, fromUser: false}, loader),
        addEngineMessage: (aiMessage) => this.appendMessage(aiMessage, {title: '', content: '', fromUser: false, onlyAddToMessages: true}, loader),
        addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => this.appendMessage('' as any, {title: '', content: msg, fromUser: fromUser, uiOnly: true, messageOptions: messageOptions}, loader),
        addUserMessage: (aiMsg, content) => this.appendMessage(aiMsg, {title: '', content: content, fromUser: true}, loader),
        addConfirmMessage: (msg?: string) => this.appendMessage('' as any, {title: '', content: '', fromUser: false, uiOnly: true, messageOptions: {confirm: {message: msg}}}, loader)!.confirmPromie,
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
      prevMessages: [...this._messages],
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
        this.appendMessage(null as any, {title: msg.title ?? '', content: msg.text, fromUser: msg.fromUser, uiOnly: true, messageOptions: msg.messageOptions}); // no loader
      });
      this.afterConversationLoad(conv);
      //grok.shell.info(`Loaded conversation: ${conv.initialPrompt.substring(0, 50)}...`);
    } catch (error) {
      console.error('Failed to load conversation:', error);
      grok.shell.error('Failed to load conversation');
    }
  }

  protected afterConversationLoad(conversation: StoredConversationWithContext<T>) {
    // do nothing in base class
  }

  protected getConversationMeta(): any {
    return undefined;
  }

  public async saveCurrentConversation(initialPrompt?: string) {
    if (this._messages.length === 0) return;

    try {
      if (this.currentConversationId)
        await ConversationStorage.updateConversation(this.currentConversationId, this._messages, this._uiMessages, this.getConversationMeta());
      else {
        this.currentConversationId = await ConversationStorage.saveConversation(
          this._messages, this._uiMessages, initialPrompt ?? this._uiMessages[0]?.text ?? 'No Title', this.contextId, this.getConversationMeta()
        );
      }
    } catch (error) {
      console.error('Failed to save conversation:', error);
    }
  }

  private toggleSpeechRecognition() {
    if (this.isRecognizing)
      this.stopRecognition();
    else
      this.startRecognition();
  }

  private startRecognition() {
    try {
      this.recognition = new SpeechRecognition();
      this.recognition.continuous = false;
      this.recognition.interimResults = false;
      this.recognition.lang = 'en-US';

      // Update UI to show recording state
      this.isRecognizing = true;
      this.micButton.classList.remove('fa-microphone');
      this.micButton.classList.add('fa-stop');
      this.micButton.style.color = 'orangered';
      ui.setUpdateIndicator(this.textAreaDiv, true, 'Listening...');

      this.recognition.onresult = (event) => {
        const transcript = event.results[0][0].transcript;
        this.textArea.value = transcript;
        this.textArea.focus();
        ui.setUpdateIndicator(this.textAreaDiv, false);
      };

      this.recognition.onerror = (event) => {
        console.error('Speech recognition error:', event.error);
        ui.setUpdateIndicator(this.textAreaDiv, false);

        let errorMessage = 'Speech recognition error';
        switch (event.error) {
        case 'no-speech':
          errorMessage = 'No speech detected. Please try again.';
          break;
        case 'audio-capture':
          errorMessage = 'No microphone found or access denied';
          break;
        case 'not-allowed':
          errorMessage = 'Microphone access denied. Please enable microphone permissions.';
          break;
        case 'network':
          errorMessage = 'Network error during speech recognition';
          break;
        case 'aborted':
          // User stopped intentionally, no error message needed
          break;
        default:
          errorMessage = `Speech recognition error: ${event.error}`;
        }

        if (event.error !== 'aborted')
          grok.shell.error(errorMessage);


        this.stopRecognition();
      };

      this.recognition.onend = () => {
        this.isRecognizing = false;
        this.micButton.classList.remove('fa-stop');
        this.micButton.classList.add('fa-microphone');
        this.micButton.style.color = '';
        ui.setUpdateIndicator(this.textAreaDiv, false);
      };

      this.recognition.start();
    } catch (error) {
      console.error('Failed to start speech recognition:', error);
      grok.shell.error('Failed to start speech recognition');
      this.stopRecognition();
    }
  }

  private stopRecognition() {
    if (this.recognition) {
      try {
        this.recognition.stop();
      } catch (error) {
        console.error('Error stopping recognition:', error);
      }
      this.recognition = null;
    }

    this.isRecognizing = false;
    this.micButton.classList.remove('fa-stop');
    this.micButton.classList.add('fa-microphone');
    this.micButton.style.color = '';
    ui.setUpdateIndicator(this.textAreaDiv, false);
  }
}


function receiveFeedback(userPrompt: string, aiResponse: string, contextId: string, helpful: boolean) {
  // not implemented yet
}


export class DBAIPanel extends AIPanel<MessageType, DBAIPanelInputs> {
  protected get placeHolder() { return 'Ask your database, like "Total sales by regions"'; }
  protected catalogInput: DG.InputBase<string>;
  constructor(catalogs: string[], defaultCatalog: string, connectionID: string, view: DG.View | DG.ViewBase) {
    super(connectionID, view); // context ID is connection ID
    this.catalogInput = ui.input.choice('Catalog', {
      items: catalogs,
      value: defaultCatalog,
      nullable: false,
      tooltipText: 'Select the database catalog to use for AI-assisted query generation.',
    }) as DG.InputBase<string>;
    this.inputControlsDiv.appendChild(this.catalogInput.input);
    ui.tooltip.bind(this.catalogInput.input, 'Select the database catalog to use for AI-assisted query generation.');
  }

  public getCurrentInputs(): DBAIPanelInputs {
    const baseInputs = super.getCurrentInputs();
    return {
      ...baseInputs,
      catalogName: this.catalogInput.value!,
    };
  }
}

export class TVAIPanel extends AIPanel<MessageType, TVAIPanelInputs> {
  protected get placeHolder() { return 'Ask about your table data...'; }
  protected tableView: DG.TableView;
  constructor(view: DG.TableView) {
    super(view.dataFrame?.name ?? view.name ?? 'AI-Table-context', view); // context ID is table name
    this.tableView = view;
  }

  // meta of TVAIPanel will store the layout of the table view
  protected getConversationMeta() {
    return this.tableView.saveLayout().viewState;
  }

  protected afterConversationLoad(conversation: StoredConversationWithContext<MessageType>) {
    const currentViewers = Array.from(this.tableView.viewers);
    // only apply layout if the view is standard grid and there is layout meta
    if (!!conversation.meta && currentViewers.length === 1 && currentViewers[0].type === DG.VIEWER.GRID) {
      const layout = DG.ViewLayout.fromViewState(conversation.meta);
      this.tableView.loadLayout(layout, true);
    }
  }
}

export class ScriptingAIPanel extends AIPanel<MessageType, ScriptingAIPanelInputs> {
  protected get placeHolder() { return 'Ask AI to generate a script...'; }
  protected languageInput: DG.InputBase<string>;

  constructor(view: DG.View | DG.ViewBase) {
    super('scripting-ai-panel', view); // context ID is fixed for scripting panel
    this.languageInput = ui.input.choice('Language', {
      items: Object.values(DG.SCRIPT_LANGUAGE),
      value: DG.SCRIPT_LANGUAGE.JAVASCRIPT,
      nullable: false,
      tooltipText: 'Select scripting language for the generated script.',
    }) as DG.InputBase<string>;
    this.inputControlsDiv.appendChild(this.languageInput.input);
    ui.tooltip.bind(this.languageInput.input, 'Select scripting language for the generated script.');
  }

  public getCurrentInputs(): ScriptingAIPanelInputs {
    const baseInputs = super.getCurrentInputs();
    return {
      ...baseInputs,
      language: this.languageInput.value as DG.ScriptingLanguage,
    };
  }
}
