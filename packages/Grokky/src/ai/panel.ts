/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {delay} from 'rxjs/operators'
// @ts-ignore .... idk why it does not like it
import '../../css/ai.css';
import {dartLike, fireAIAbortEvent, getAIPanelToggleSubscription, createStyledMarkdown, isEnterKey, copyToClipboard} from '../utils';
import {buildViewContext, renderEntityBlocks} from '../claude/exec-blocks';
import {ConversationStorage, StoredConversationWithContext} from './storage';
import {ClaudeRuntimeClient} from '../claude/runtime-client';
import {resolveScopes, showSuggestionsMenu} from './prompt-suggestions';

export type MessageType = {role: string; content: any};

type AIPanelInputs = {
    prompt: string,
}

type DBAIPanelInputs = AIPanelInputs & {
    catalogName: string,
}

export type ScriptingAIPanelInputs = AIPanelInputs & {
  language: DG.ScriptingLanguage
};


type TVAIPanelInputs = AIPanelInputs;

export interface AskUserOption {
  label: string;
  description?: string;
}

export interface AskUserQuestion {
  question: string;
  header?: string;
  options: AskUserOption[];
  multiSelect?: boolean;
}

interface AskUserInput {
  questions: AskUserQuestion[];
}

interface AskUserResponse {
  questions: AskUserQuestion[];
  answers: Record<string, string>;
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
  /** if set on a user message, shows a small green check ("Handled natively") instead of a response block */
  handledNatively?: boolean,
  /** if set, renders the message as a centered system event (retry notice, workflow header, etc.) */
  system?: boolean,
}

export interface UIMessage {
  fromUser: boolean;
  text: string;
  title?: string;
  messageOptions?: UIMessageOptions;
}

export type PanelMessageRet = {
  confirmPromise: Promise<boolean>,
}

export type AIPanelFuncs<T extends MessageType = MessageType> = {
  /** Adds @aiMsg to the message stack (from user) and @msg to the UI panel */
  addUserMessage: (aiMsg: T, msg: string) => void,
  /** Adds @aiMsg to the message stack (from AI) and @msg with @title to the UI panel*/
  addAIMessage: (aiMsg: T, title: string, msg: string) => void,
  /** Adds @aiMsg to the message stack without showing it to the ui  */
  addEngineMessage: (aiMsg: T) => void,
  /** Adds @msg to the UI without adding it to the AI message stack */
  addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => void
  /** Marks the last user prompt as handled by Datagrok's built-in handler: shows a green check, no response block */
  markHandledNatively: () => void,
  /**Adds confirmation section to the panel and awaits result */
  addConfirmMessage: (msg?: string) => Promise<boolean>,
  /** Shows options to the user and returns their selection */
  showInputRequest?: (input: AskUserInput) => Promise<AskUserResponse | null>,
}

export interface StreamingPanel<T extends MessageType = MessageType> {
  sessionId: string;
  startChatSession(): {session: AIPanelFuncs<T>, endSession: () => void, loader: HTMLElement};
  prependViewContext(prompt: string, view: DG.ViewBase): string;
  prependEntityContext(prompt: string): string;
  updateStreaming(content: string, loader: HTMLElement): void;
  finalizeStreaming(displayContent: string, execContent: string, view: DG.ViewBase): Promise<void>;
  appendStreamedElement(el: HTMLElement): void;
  appendUiMessage(content: string): void;
  clearStreaming(): void;
  showInputRequest(input: any): Promise<any>;
  cancelInputRequest(): void;
  get rawRender(): boolean;
  get noPrompt(): boolean;
  enableNoPrompt(): void;
}

export class AIPanel<T extends MessageType = MessageType, K extends AIPanelInputs = AIPanelInputs> implements StreamingPanel<T> {
  private root: HTMLElement;
  protected view: DG.View | DG.ViewBase;
  private inputArea: HTMLElement;
  protected header: HTMLElement;
  protected outputArea: HTMLElement;
  protected textArea: HTMLTextAreaElement;
  private textAreaDiv: HTMLElement;
  private runButton: HTMLElement;
  private newChatButton: HTMLElement;
  private copyConversationButton: HTMLElement;
  private historyButton: HTMLElement;
  private tryAgainButton: HTMLElement;
  private micButton: HTMLElement;
  private rawRenderButton: HTMLElement;
  private wandButton: HTMLElement;
  private _rawRender: boolean = false;
  private _noPrompt: boolean = false;
  private recognition: SpeechRecognition | null = null;
  private isRecognizing: boolean = false;
  /** `Say "cancel" to stop` caption shown next to the loader while the AI is working in voice mode. */
  private _voiceCancelHint: HTMLElement | null = null;
  private _onRunRequest = new rxjs.Subject<{prevMessages: T[], currentPrompt: K}>();
  protected _messages: T[] = [];
  protected _uiMessages: UIMessage[] = [];
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
  protected attachedEntities: DG.Entity[] = [];
  protected attachmentsRow: HTMLElement = ui.divH([], {style: {flexWrap: 'wrap'}});
  private _pendingEntityContext = '';
  private _pendingAttachmentsForRender: DG.Entity[] = [];
  protected get contextId(): string {
    return this._contextID;// these should be overriden in subclasses
  }

  protected _streamingContainer: HTMLElement | null = null;
  protected _streamingMarkdownEl: HTMLElement | null = null;
  private _sessionId: string;
  private _pendingInputResolve: ((value: AskUserResponse | null) => void) | null = null;
  private _skillMenu: DG.Menu | null = null;
  private _inline: boolean = false;
  /** Index into {@link promptHistory} while cycling with Ctrl+[ / Ctrl+]; `null` means the live draft is shown. */
  private _promptHistoryIndex: number | null = null;
  /** The unsubmitted draft saved when the user starts cycling, restored when they cycle back past the newest entry. */
  private _promptDraft: string = '';

  get sessionId(): string { return this._sessionId; }


  private currentConversationId: string | null = null;
  constructor(private _contextID: string = 'global-ai-panel', view: DG.View | DG.ViewBase, opts: {inline?: boolean} = {}) {
    this.view = view;
    this._sessionId = `claude-${_contextID}-${crypto.randomUUID()}`;
    this._inline = !!opts.inline;
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
    }, 'Send');
    this.textArea.addEventListener('keydown', (event: KeyboardEvent) => {
      if (event.key === 'Escape' && this._skillMenu) {
        this._skillMenu.hide();
        this._skillMenu = null;
        return;
      }
      // Ctrl+[ / Ctrl+] cycle through this session's prompt history (keyCode fallback for Chrome ≤ 50).
      if (event.ctrlKey && (event.key === '[' || event.keyCode === 219)) {
        event.preventDefault();
        this.navigatePromptHistory(-1);
        return;
      }
      if (event.ctrlKey && (event.key === ']' || event.keyCode === 221)) {
        event.preventDefault();
        this.navigatePromptHistory(1);
        return;
      }
      if (isEnterKey(event) && (!event.ctrlKey && !event.metaKey)) {
        event.preventDefault();
        event.stopImmediatePropagation();
        this.handleRun();
      }
    });
    this.textArea.addEventListener('input', () => { this._promptHistoryIndex = null; this._updateSkillMenu(); });
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
    this.rawRenderButton = ui.iconFA('terminal', () => {
      this._rawRender = !this._rawRender;
      this.rawRenderButton.style.color = this._rawRender ? 'var(--blue-1)' : '';
      this.root.classList.toggle('d4-ai-raw-mode', this._rawRender);
    }, 'Toggle raw console');
    this.wandButton = ui.iconFA('magic', async (e) => {
      const scopes = await resolveScopes('panel', this.view);
      showSuggestionsMenu(scopes, (prompt) => this.runSuggestion(prompt), e);
    }, 'Prompt suggestions');
    this.wandButton.classList.add('grokky-search-wand');
    this.setWandVisible(true);
    this.hideContentIcons();
    this.inputControlsDiv = ui.divH([
      this.wandButton, this.micButton, this.rawRenderButton,
    ], 'd4-ai-panel-input-controls');
    this.runButton.style.color = 'var(--blue-1)';
    const sessionControls = ui.divH([this.copyConversationButton, this.historyButton, this.newChatButton], 'd4-ai-panel-run-controls');
    sessionControls.style.marginLeft = 'auto';
    const messageControls = ui.divH([this.tryAgainButton, this.runButton], 'd4-ai-panel-run-controls');
    const controlsDiv = ui.divH([this.inputControlsDiv, sessionControls, ui.div([], 'd4-ribbon-separator'), messageControls], 'd4-ai-panel-controls-container');
    this.textAreaDiv = ui.divV([this.attachmentsRow, this.textArea], {classes: 'd4-ai-input-textarea-div', style: {position: 'relative'}});
    ui.makeDroppable(this.textAreaDiv, {
      acceptDrop: (o) => o instanceof DG.Entity,
      doDrop: (args: any) => {
        if (args?.dragObject instanceof DG.Entity) {
          this.addEntityChip(args.dragObject);
          this.textArea.focus();
        }
      },
    });
    this.inputArea.appendChild(this.textAreaDiv);
    this.inputArea.appendChild(controlsDiv);
    this.root.appendChild(this.header);
    this.root.appendChild(this.outputArea);
    this.root.appendChild(this.inputArea);

    if (this._inline)
      this.root.classList.add('d4-ai-inline-mode');

    this.setupSubscriptions();
  }

  protected setupSubscriptions() {
    if (this._inline)
      return;
    let wasShown = false;
    // Defer the dock/undock until after the view switch has settled.
    const sub = grok.events.onCurrentViewChanged.pipe(delay(150)).subscribe(() => {
      if (grok.shell.v != this.view) {
        wasShown = this.isShown;
        if (wasShown)
          this.hide();
      } else if (wasShown)
        this.show();
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

  mountInto(parent: HTMLElement) {
    parent.appendChild(this.root);
  }

  show(focus: boolean = false) {
    // Opening the AI panel steals focus from the active element (e.g., query editor)
    const previouslyFocused = document.activeElement as HTMLElement | null;
    const aiContainer = grok.shell.windows.ai;
    if (!aiContainer.contains(this.root))
      aiContainer.appendChild(this.root);
    grok.shell.windows.showAI = true;
    this.renderEmptyState();
    if (focus)
      this.textArea.focus();
    else
      previouslyFocused?.focus();
  }

  hide() {
    if (grok.shell.windows.ai.contains(this.root))
      grok.shell.windows.showAI = false;
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
    return copyToClipboard(this.formatConversation());
  }

  toggle() {
    this.isShown ? this.hide() : this.show(true);
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
    } as K;
  }

  private addEntityChip(e: DG.Entity): void {
    if (this.attachedEntities.some((x) => x.id === e.id))
      return;
    this.attachedEntities.push(e);
    const icon = DG.ObjectHandler.forEntity(e)?.renderIcon(e.dart) ?? ui.iconFA('tag');
    const label = ui.label(e.friendlyName ?? e.name);
    const close = ui.iconFA('times', () => {
      this.attachedEntities = this.attachedEntities.filter((y) => y !== e);
      chip.remove();
    }, 'Remove');
    const chip = ui.divH([icon, label, close], 'grokky-entity-chip');
    this.attachmentsRow.appendChild(chip);
  }

  private clearAttachments(): void {
    this.attachedEntities = [];
    ui.empty(this.attachmentsRow);
  }

  private describeEntity(e: DG.Entity): string {
    const parts = [`type: ${e.entityType}`, `id: ${e.id}`];
    if (e.nqName)
      parts.push(`nqName: ${e.nqName}`);
    if (e instanceof DG.FileInfo)
      parts.push(`path: ${e.fullPath}`);
    return `- "${e.friendlyName ?? e.name}" (${parts.join(', ')})`;
  }

  public prependEntityContext(prompt: string): string {
    const ctx = this._pendingEntityContext;
    if (!ctx)
      return prompt;
    this._pendingEntityContext = '';
    return ctx + '\n---\n\n' + prompt;
  }

  protected _aiMessagesAccordionPane: HTMLElement | null = null;
  private _lastUserPromptContainer: HTMLElement | null = null;

  /** Ensures there is an open response block for the current turn. The block holds all consecutive
   * AI messages and exposes a hover-only minimize icon on its left that collapses it to one line. */
  protected ensureResponseBlock(): void {
    if (this._aiMessagesAccordionPane)
      return;
    this._aiMessagesAccordionPane = ui.divV([], 'd4-ai-messages-accordion-pane');
    const minimizeIcon = ui.iconFA('window-minimize', null, 'Minimize');
    minimizeIcon.classList.add('d4-ai-response-minimize-icon');
    const block = ui.divH([minimizeIcon, this._aiMessagesAccordionPane], 'd4-ai-response-block');
    let collapsed = false;
    minimizeIcon.onclick = () => {
      collapsed = !collapsed;
      block.classList.toggle('d4-ai-response-block-collapsed', collapsed);
      minimizeIcon.classList.toggle('fa-window-minimize', !collapsed);
      minimizeIcon.classList.toggle('fa-window-maximize', collapsed);
    };
    ui.tooltip.bind(minimizeIcon, () => collapsed ? 'Expand' : 'Minimize');
    this.outputArea.appendChild(block);
  }

  protected createStyledMarkdown(content: string): HTMLElement {
    return createStyledMarkdown(content);
  }

  private createHandledNativelyIcon(): HTMLElement {
    const icon = ui.iconFA('check', null, 'Handled natively');
    icon.classList.add('d4-ai-handled-natively-icon');
    return icon;
  }

  /** Marks the most recent user prompt as handled by Datagrok's built-in handler:
   * shows a small green check (tooltip "Handled natively") and skips the response block entirely. */
  public markPromptHandledNatively(): void {
    const last = this._uiMessages[this._uiMessages.length - 1];
    if (!last?.fromUser)
      return;
    last.messageOptions = {...last.messageOptions, handledNatively: true};
    if (this._lastUserPromptContainer && !this._lastUserPromptContainer.querySelector('.d4-ai-handled-natively-icon'))
      this._lastUserPromptContainer.insertBefore(this.createHandledNativelyIcon(), this._lastUserPromptContainer.firstChild);
  }

  protected appendFeedbackButtons(markDown: HTMLElement, onFeedback?: (helpful: boolean) => void): void {
    const feedbackDiv = ui.divH([], 'd4-ai-panel-feedback-div');
    const thumbsUp = ui.iconFA('thumbs-up', () => {
      grok.shell.info('Thanks for your feedback!');
      handleFeedback(true);
    }, 'Helpful');
    const thumbsDown = ui.iconFA('thumbs-down', () => {
      grok.shell.info('Thanks for your feedback!');
      handleFeedback(false);
    }, 'Not Helpful');
    const copyMsg = ui.iconFA('copy', () => {
      copyToClipboard(markDown.textContent || '').then((ok) =>
        ok ? grok.shell.info('Message copied to clipboard') : grok.shell.error('Failed to copy message'));
    }, 'Copy Message');
    [copyMsg, thumbsUp, thumbsDown].forEach((el) => dartLike(el.style).set('padding', '2px').set('borderRadius', '6px'));
    function handleFeedback(helpful: boolean) {
      [thumbsUp, thumbsDown].forEach((el) => { el.style.backgroundColor = ''; });
      (helpful ? thumbsUp : thumbsDown).style.backgroundColor = helpful ? 'rgba(0, 150, 30, 0.2)' : 'rgba(200, 0, 0, 0.2)';
      onFeedback?.(helpful);
    }
    feedbackDiv.appendChild(copyMsg);
    feedbackDiv.appendChild(thumbsUp);
    feedbackDiv.appendChild(thumbsDown);
    dartLike(feedbackDiv.style).set('alignItems', 'center').set('width', '100%').set('paddingBottom', '8px').set('paddingLeft', '4px');
    markDown.appendChild(feedbackDiv);
  }

  protected appendMessage(
    aiMessage: T, uiMessage: {
      title: string, content: string, fromUser: boolean, onlyAddToMessages?: boolean, uiOnly?: boolean, messageOptions?: UIMessageOptions
    }, loader?: HTMLElement
  ): PanelMessageRet | undefined {
    let ret: PanelMessageRet | undefined = undefined;
    const emptyState = this.outputArea.querySelector('.grokky-empty-state');
    if (emptyState) {
      emptyState.remove();
      this.setWandVisible(true);
    }
    if (!uiMessage.uiOnly)
      this._messages.push(aiMessage);
    if (uiMessage.onlyAddToMessages)
      return;
    // from this point we know that message is also in the ui.
    this._uiMessages.push({fromUser: !!uiMessage.fromUser, text: uiMessage.content, title: uiMessage.title, messageOptions: uiMessage.messageOptions});
    if (uiMessage.fromUser) {
      const promptText = ui.divText(uiMessage.content, 'd4-ai-user-prompt-divtext');
      const userDiv = ui.div(
        uiMessage.messageOptions?.handledNatively ? [this.createHandledNativelyIcon(), promptText] : [promptText],
        'd4-ai-user-prompt-container');
      this.outputArea.appendChild(userDiv);
      if (this._pendingAttachmentsForRender.length) {
        const chips = this._pendingAttachmentsForRender.map((e) => {
          const icon = DG.ObjectHandler.forEntity(e)?.renderIcon(e.dart) ?? ui.iconFA('tag');
          const label = ui.label(e.friendlyName ?? e.name);
          return ui.divH([icon, label], 'grokky-entity-chip');
        });
        const chipsRow = ui.divH(chips, {style: {flexWrap: 'wrap', justifyContent: 'flex-end'}});
        this.outputArea.appendChild(chipsRow);
        this._pendingAttachmentsForRender = [];
      }
      this._lastUserPromptContainer = userDiv;
      this._aiMessagesAccordionPane = null; // reset accordion pane so that next AI message creates a new one
    } else if (uiMessage.messageOptions?.system) {
      this.outputArea.appendChild(ui.divText(uiMessage.content, 'grokky-system-message'));
      this._aiMessagesAccordionPane = null;
    } else {
      this.ensureResponseBlock();
      const markDown = this.createStyledMarkdown(uiMessage.content);

      if (uiMessage?.messageOptions?.finalResult && !uiMessage.fromUser) {
        this.appendFeedbackButtons(markDown, (helpful) => {
          receiveFeedback(
            this._uiMessages[0]?.text ?? '',
            uiMessage?.messageOptions?.finalResult!,
            this.contextId,
            helpful
          );
        });
      }

      if (uiMessage?.messageOptions?.confirm && !uiMessage.fromUser) {
        let resolve: (value: boolean) => void = () => {};
        const confirmPromise = new Promise<boolean>((res) => { resolve = res; });
        ret = {confirmPromise: confirmPromise};
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

      const titleText = uiMessage.title ? [ui.h3(uiMessage.title, 'd4-ai-assistant-response-title')] : [];
      const assistantDiv = ui.divV([...titleText, markDown], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane!.appendChild(assistantDiv);
    }
    // reapend loader to the end
    if (loader)
      this.outputArea.appendChild(loader);
    this.outputArea.scrollTop = this.outputArea.scrollHeight;
    this.showContentIcons();
    return ret;
  }

  public startChatSession(): {session: AIPanelFuncs<T>, endSession: () => void, loader: HTMLElement} {
    const spinner = ui.icons.loader();
    dartLike(spinner.style).set('height', '20px');
    // Speech is easy to mishear, so while the AI works in voice mode show the word that aborts the run.
    const cancelHint = ui.divText('Say "cancel" to stop', 'd4-ai-voice-cancel-hint');
    cancelHint.style.display = this.isRecognizing ? '' : 'none';
    this._voiceCancelHint = cancelHint;
    const loader = ui.divH([spinner, cancelHint], 'd4-ai-loader');
    dartLike(loader.style).set('alignSelf', 'center').set('marginTop', '8px');
    this.runButton.classList.remove('fal', 'fa-paper-plane');
    this.runButton.classList.add('fas', 'fa-stop');
    this.runButton.style.color = 'orangered';
    this.runButtonTooltip = actionButtionValues.stop;
    return {
      loader,
      session: {
        addAIMessage: (aiMessage, title, content) => this.appendMessage(aiMessage, {title: title, content: content, fromUser: false}, loader),
        addEngineMessage: (aiMessage) => this.appendMessage(aiMessage, {title: '', content: '', fromUser: false, onlyAddToMessages: true}, loader),
        addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => this.appendMessage('' as any, {title: '', content: msg, fromUser: fromUser, uiOnly: true, messageOptions: messageOptions}, loader),
        markHandledNatively: () => this.markPromptHandledNatively(),
        addUserMessage: (aiMsg, content) => this.appendMessage(aiMsg, {title: '', content: content, fromUser: true}, loader),
        addConfirmMessage: (msg?: string) => this.appendMessage('' as any, {title: '', content: '', fromUser: false, uiOnly: true, messageOptions: {confirm: {message: msg}}}, loader)!.confirmPromise,
        showInputRequest: (input: AskUserInput) => this.showInputRequest(input),
      },
      endSession: () => {
        this.runButton.classList.remove('fas', 'fa-stop');
        this.runButton.classList.add('fal', 'fa-paper-plane');
        this.runButton.style.color = 'var(--blue-1)';
        this.runButtonTooltip = actionButtionValues.run;
        this._voiceCancelHint = null;
        this.saveCurrentConversation().catch((e) => console.error('Failed to save conversation before hiding panel:', e));
        loader.remove();
      }
    };
  }

  get isShown(): boolean {
    return grok.shell.windows.showAI && grok.shell.windows.ai.contains(this.root);
  }

  get rawRender(): boolean { return this._rawRender; }
  get noPrompt(): boolean { return this._noPrompt; }

  public appendArtifact(node: HTMLElement): void {
    this.outputArea.appendChild(ui.divV([node], 'd4-ai-assistant-response-container'));
    this.showContentIcons();
  }

  enableNoPrompt(): void {
    this._noPrompt = true;
    if (!this._rawRender) {
      this._rawRender = true;
      this.rawRenderButton.style.color = 'var(--blue-1)';
      this.root.classList.add('d4-ai-raw-mode');
    }
    this.resetSession();
    this.handleClear();
  }

  resetSession(): void {
    this._sessionId = `claude-${crypto.randomUUID()}`;
    this._streamingContainer = null;
    this._streamingMarkdownEl = null;
  }

  prependViewContext(prompt: string, view: DG.ViewBase): string {
    const ctx = buildViewContext(view);
    if (!ctx)
      return prompt;
    return ctx + '\n---\n\n' + prompt;
  }

  private createStreamingEl(content: string): HTMLElement {
    if (this._rawRender) {
      const pre = document.createElement('pre');
      pre.textContent = content;
      pre.style.cssText = 'white-space:pre-wrap;user-select:text;margin:0';
      return pre;
    }
    const md = ui.markdown(content);
    dartLike(md.style).set('userSelect', 'text').set('maxWidth', '100%');
    return md;
  }

  updateStreaming(content: string, loader: HTMLElement): void {
    if (!this._streamingContainer) {
      loader.style.display = 'none';
      this.ensureResponseBlock();
      this._streamingMarkdownEl = this.createStreamingEl(content);
      this._streamingContainer = ui.divV([this._streamingMarkdownEl], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane!.appendChild(this._streamingContainer);
    } else {
      const el = this.createStreamingEl(content);
      this._streamingMarkdownEl!.replaceWith(el);
      this._streamingMarkdownEl = el;
    }
    this.outputArea.scrollTop = this.outputArea.scrollHeight;
  }

  async finalizeStreaming(displayContent: string, _execContent: string, _view: DG.ViewBase): Promise<void> {
    if (this._rawRender) {
      this._streamingMarkdownEl = null;
      this._streamingContainer = null;
      this._uiMessages.push({fromUser: false, text: displayContent, messageOptions: {finalResult: displayContent}});
      return;
    }
    this.renderFinalContent(displayContent);
  }

  public appendStreamedElement(el: HTMLElement): void {
    this.ensureResponseBlock();
    this._aiMessagesAccordionPane!.appendChild(ui.divV([el], 'd4-ai-assistant-response-container'));
  }

  public appendUiMessage(content: string): void {
    this.appendMessage('' as any, {title: '', fromUser: false, uiOnly: true, content, messageOptions: {system: true}});
  }

  protected renderFinalContent(content: string): void {
    const markDown = this.createStyledMarkdown(content);
    renderEntityBlocks(markDown);
    this.appendFeedbackButtons(markDown);

    if (this._streamingMarkdownEl) {
      this._streamingMarkdownEl.replaceWith(markDown);
      this._streamingMarkdownEl = null;
      this._streamingContainer = null;
    } else {
      this.ensureResponseBlock();
      this._aiMessagesAccordionPane!.appendChild(ui.divV([markDown], 'd4-ai-assistant-response-container'));
    }

    this._uiMessages.push({fromUser: false, text: content, messageOptions: {finalResult: content}});
  }

  clearStreaming(): void {
    this._streamingContainer?.remove();
    this._streamingContainer = null;
    this._streamingMarkdownEl = null;
  }

  showInputRequest(input: AskUserInput): Promise<AskUserResponse | null> {
    const questions = input.questions ?? [];
    let resolved = false;
    const choiceInputs: DG.InputBase<string>[] = [];

    for (const q of questions) {
      const items = q.options.map((o) => o.label);
      choiceInputs.push(ui.input.choice(q.header ?? '', {items, value: items[0], nullable: false}) as DG.InputBase<string>);
    }
    const form = ui.divV(questions.map((q, i) => ui.divV([ui.divText(q.question), choiceInputs[i].root])));

    return new Promise<AskUserResponse | null>((resolve) => {
      const doResolve = (value: AskUserResponse | null) => {
        if (resolved)
          return;
        resolved = true;
        this._pendingInputResolve = null;
        for (const inp of choiceInputs)
          (inp.input as HTMLSelectElement).disabled = true;
        resolve(value);
      };

      form.appendChild(ui.button('Submit', () => {
        const answers: Record<string, string> = {};
        for (let i = 0; i < questions.length; i++)
          answers[questions[i].question] = choiceInputs[i].value!;
        doResolve({questions, answers});
      }));

      this._pendingInputResolve = doResolve;
      this.ensureResponseBlock();
      const wrapper = ui.divV([form], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane!.appendChild(wrapper);
      this.outputArea.scrollTop = this.outputArea.scrollHeight;
    });
  }

  cancelInputRequest(): void {
    if (this._pendingInputResolve) {
      this._pendingInputResolve(null);
      this._pendingInputResolve = null;
    }
  }

  private tryAgain() {
    if (this._messages.length === 0)
      return; // should never happen, but just in case
    if (this.runButtonTooltip === actionButtionValues.stop)
      return;
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

  /** Prompts the user submitted in this session, oldest first — the navigable input history. */
  private get promptHistory(): string[] {
    return this._uiMessages.filter((m) => m.fromUser).map((m) => m.text);
  }

  /** Replaces the input with the previous (-1) or next (+1) prompt from this session's history. */
  private navigatePromptHistory(direction: -1 | 1): void {
    const hist = this.promptHistory;
    if (hist.length === 0)
      return;
    if (this._promptHistoryIndex === null) {
      if (direction === 1)
        return; // already at the live draft — nothing newer
      this._promptDraft = this.textArea.value;
      this._promptHistoryIndex = hist.length;
    }
    const next = this._promptHistoryIndex + direction;
    if (next < 0)
      return; // already at the oldest entry
    if (next >= hist.length) {
      this._promptHistoryIndex = null;
      this.textArea.value = this._promptDraft;
    } else {
      this._promptHistoryIndex = next;
      this.textArea.value = hist[next];
    }
    this.textArea.selectionStart = this.textArea.selectionEnd = this.textArea.value.length;
    this._updateSkillMenu();
  }

  protected handleRun() {
    if (this._pendingInputResolve)
      return;
    if (this.runButtonTooltip === actionButtionValues.stop)
      return;
    if (!this.textArea.value.trim())
      return;
    const inputs = this.getCurrentInputs();
    this.textArea.value = '';
    this._pendingEntityContext = this.attachedEntities.length ?
      'Attached Datagrok entities (use MCP tools to fetch full details by id/nqName/path):\n' +
        this.attachedEntities.map((e) => this.describeEntity(e)).join('\n') :
      '';
    this._pendingAttachmentsForRender = [...this.attachedEntities];
    this.clearAttachments();
    this._promptHistoryIndex = null;
    this._onRunRequest.next({
      prevMessages: this._messages,
      currentPrompt: inputs,
    });
  }

  private _hideSkillMenu(): void {
    if (this._skillMenu) {
      this._skillMenu.hide();
      this._skillMenu = null;
    }
  }

  private _updateSkillMenu(): void {
    const text = this.textArea.value;
    if (!text.startsWith('/'))
      return this._hideSkillMenu();

    const query = text.slice(1).toLowerCase();
    const allNames = ClaudeRuntimeClient.getInstance().getSkillNames();
    if (!allNames.length)
      return this._hideSkillMenu();

    const filtered = allNames.filter((name) => name.toLowerCase().includes(query));
    if (!filtered.length)
      return this._hideSkillMenu();

    this._hideSkillMenu();
    this._skillMenu = DG.Menu.popup();
    this._skillMenu.header('Skills');
    for (const name of filtered) {
      this._skillMenu.item(name, () => {
        this.textArea.value = `/${name} `;
        this.textArea.focus();
        this._skillMenu = null;
      });
    }
    this._skillMenu.show({element: this.textAreaDiv, y: 0});
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

  protected handleClear() {
    this._messages = [];
    this._uiMessages = [];
    this._promptHistoryIndex = null;
    this._lastUserPromptContainer = null;
    this.outputArea.innerHTML = '';
    this._onClearChatRequest.next();
    this.hideContentIcons();
    this.renderEmptyState();
  }

  protected runSuggestion(prompt: string): void {
    this.textArea.value = prompt;
    this.handleRun();
  }

  protected shouldShowEmptyState(): boolean {
    return this.view instanceof DG.TableView &&
      this._uiMessages.length === 0 &&
      !this.outputArea.querySelector('.grokky-empty-state');
  }

  private setWandVisible(visible: boolean): void {
    this.wandButton.style.display = visible && this.view instanceof DG.TableView ? '' : 'none';
  }

  protected async renderEmptyState(): Promise<void> {
    if (!this.shouldShowEmptyState())
      return;
    const scopes = await resolveScopes('panel', this.view);
    if (!this.shouldShowEmptyState())
      return;

    const blocks = scopes.map((s) => {
      const icon = ui.iconFA(s.icon ?? 'circle');
      icon.classList.add('grokky-scope-icon');
      if (s.key)
        icon.classList.add(`grokky-scope-${s.key}`);
      const header = ui.h3(ui.span([icon, s.label]));
      const cards = s.suggestions.slice(0, 2).map((sg) => {
        const card = ui.card(ui.divText(sg.label ?? sg.prompt));
        card.onclick = () => this.runSuggestion(sg.prompt);
        return card;
      });
      return ui.divV([header, ui.divH(cards)]);
    });

    const root = ui.panel([ui.h2('What can I help you with?'), ...blocks], 'grokky-empty-state');
    this.outputArea.appendChild(root);
    this.setWandVisible(false);
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
      this._promptHistoryIndex = null;
      this._lastUserPromptContainer = null;
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

  /** The voice escape hatch only matters while the mic is live — keep the loader caption in sync. */
  private syncVoiceCancelHint() {
    if (this._voiceCancelHint)
      this._voiceCancelHint.style.display = this.isRecognizing ? '' : 'none';
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
      this.syncVoiceCancelHint();

      this.recognition.onresult = (event) => {
        const transcript = event.results[0][0].transcript;
        ui.setUpdateIndicator(this.textAreaDiv, false);

        const command = transcript.trim().replace(/[.,;:!?]+$/, '').toLowerCase();
        const isCancelWord = command === 'stop' || command === 'cancel';

        // While a prompt is being processed, "stop"/"cancel" aborts the AI run instead of being sent
        // as a new prompt — speech is easy to mishear, so keep this escape hatch reliable. Anything
        // else said while busy is ignored (a new prompt can't be submitted yet anyway).
        if (this.runButtonTooltip === actionButtionValues.stop) {
          if (isCancelWord)
            this.terminate();
          return;
        }

        if (isCancelWord) {
          this.stopRecognition();
          this.textArea.focus();
          return;
        }

        this.textArea.value = transcript;
        this.handleRun();
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
        if (this.isRecognizing) {
          this.startRecognition();
          return;
        }
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
    this.isRecognizing = false;
    this.syncVoiceCancelHint();

    if (this.recognition) {
      try {
        this.recognition.stop();
      } catch (error) {
        console.error('Error stopping recognition:', error);
      }
      this.recognition = null;
    }
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
  private setAndRunFunc: (query: string) => void;

  constructor(catalogs: string[], defaultCatalog: string, connectionID: string, view: DG.View | DG.ViewBase, setAndRunFunc: (query: string) => void) {
    super(connectionID, view); // context ID is connection ID
    this.setAndRunFunc = setAndRunFunc;
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

  async finalizeStreaming(displayContent: string, execContent: string, _view: DG.ViewBase): Promise<void> {
    this.renderFinalContent(displayContent);
    // Extract SQL from fenced code blocks and inject into query editor
    const sqlMatch = /```(?:sql)?\n([\s\S]*?)```/.exec(execContent);
    if (sqlMatch) {
      const sql = sqlMatch[1].trimEnd().replace(/;+$/, '');
      this.setAndRunFunc(sql);
    }
  }
}

export class TVAIPanel extends AIPanel<MessageType, TVAIPanelInputs> {
  protected get placeHolder() { return 'Ask Claude about your data...'; }
  protected tableView: DG.TableView;

  constructor(view: DG.TableView) {
    super(view.dataFrame?.name ?? view.name ?? 'AI-Table-context', view);
    this.tableView = view;
  }

  protected getConversationMeta() {
    return {viewState: this.tableView.saveLayout().viewState, sessionId: this.sessionId};
  }

  protected afterConversationLoad(conversation: StoredConversationWithContext<MessageType>) {
    if (conversation.meta?.sessionId)
      (this as any)._sessionId = conversation.meta.sessionId;
    const viewState = conversation.meta?.viewState ?? conversation.meta;
    const currentViewers = Array.from(this.tableView.viewers);
    if (!!viewState && currentViewers.length === 1 && currentViewers[0].type === DG.VIEWER.GRID) {
      const layout = DG.ViewLayout.fromViewState(viewState);
      this.tableView.loadLayout(layout, true);
    }
  }
}

export class ShellAIPanel extends AIPanel {
  protected get placeHolder() { return 'Ask AI anything...'; }

  constructor() {
    super('shell-ai-panel', null as any);
  }

  protected setupSubscriptions(): void {
    // Shell panel is not tied to a view — no view-change tracking needed
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

  async finalizeStreaming(displayContent: string, execContent: string, _view: DG.ViewBase): Promise<void> {
    this.renderFinalContent(displayContent);
    // Extract code from datagrok-exec blocks and set on the script editor
    const codeMatch = /```datagrok-exec\n([\s\S]*?)```/.exec(execContent);
    if (codeMatch) {
      ui.setUpdateIndicator(this.view.root, true, 'Updating script...');
      const indicator = this.view.root.querySelector('.d4-update-shadow') as HTMLElement;
      if (indicator)
        indicator.style.zIndex = '1000';
      (this.view as DG.ScriptView).code = codeMatch[1].trimEnd();
      ui.setUpdateIndicator(this.view.root, false);
    }
  }
}
