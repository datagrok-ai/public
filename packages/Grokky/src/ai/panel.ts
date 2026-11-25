/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
// @ts-ignore .... idk why it does not like it
import '../../css/ai.css';
import {dartLike, fireAIAbortEvent, createStyledMarkdown, codeCard, tryMarkdown, plainTextEl,
  copyToClipboard} from '../utils';
import {buildWorkspaceContext, executeSingleBlock} from '../claude/exec-blocks';
import {ConversationStorage, StoredConversationWithContext} from './storage';
import {ClaudeRuntimeClient} from '../claude/runtime-client';
import {resolveScopes, showSuggestionsMenu, runSuggestionAction, Suggestion, Block, ChoiceOption} from './prompt-suggestions';
import {AiFeature} from '../claude/ai-feature';

export type MessageType = {role: string; content: any};

type AIPanelInputs = {
    prompt: string,
    isAnswer?: boolean,
}


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

interface QuestionSection {
  strip: HTMLElement;
  body: HTMLElement;
}

export type StripStyle = 'success' | 'muted';

const actionButtionValues = {
  run: 'Run AI Prompt',
  stop: 'Stop AI Generation',
} as const;

/** Silence after which a spoken phrase is sent as the prompt. */
const VOICE_SUBMIT_DELAY_MS = 1500;

const micTooltips = {
  default: 'Voice Input',
  accessDenied: 'Microphone access denied. Please enable microphone permissions.',
  noDevice: 'No microphone found or access denied',
  notSupported: 'Voice input is not supported in this browser',
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
  execCode?: string[];
}

type ExecInfo = {codes: string[], ranThisSession: boolean, view?: DG.ViewBase};
const ExecState = {idle: 'idle', ran: 'ran', error: 'error'} as const;
type ExecState = typeof ExecState[keyof typeof ExecState];
const EXEC_TOOLTIPS: Record<ExecState, string> = {
  idle: 'Show code — not run this session',
  ran: 'Show code — already executed',
  error: 'Show code — last re-run failed',
};

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
  showInputRequest?: (input: AskUserInput) => void,
}

export interface StreamingPanel<T extends MessageType = MessageType> {
  sessionId: string;
  startChatSession(): {session: AIPanelFuncs<T>, endSession: () => void};
  prependViewContext(prompt: string, view: DG.ViewBase): string;
  prependEntityContext(prompt: string): string;
  updateStreaming(content: string): void;
  setStatus(text: string): void;
  showError(text: string): void;
  finalizeStreaming(displayContent: string, execCodes: string[], view: DG.ViewBase): Promise<void>;
  appendStreamedElement(el: HTMLElement): void;
  appendUiMessage(content: string): void;
  clearStreaming(): void;
  showInputRequest(input: AskUserInput): void;
  dismissQuestion(): void;
  restoredTranscript(): string;
  get rawRender(): boolean;
  get noPrompt(): boolean;
  enableNoPrompt(): void;
  pushNativeContext(prompt: string): void;
  flushNativeContext(): string;
}

export class AIPanel<T extends MessageType = MessageType, K extends AIPanelInputs = AIPanelInputs> implements StreamingPanel<T> {
  readonly root: HTMLElement;
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
  private _pendingNativeContext: string[] = [];
  private recognition: SpeechRecognition | null = null;
  private isRecognizing: boolean = false;
  /** The final results of the phrase being spoken, not yet submitted. */
  private voiceTranscript: string = '';
  private voiceSubmitTimer?: number;
  /** `Say "cancel" to stop` caption shown next to the loader while the AI is working in voice mode. */
  private _voiceCancelHint: HTMLElement | null = null;
  private micAccessDenied: boolean = false;
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
  private get isBusy(): boolean { return this.runButtonTooltip === actionButtionValues.stop; }

  private syncConversationButtons(): void {
    const buttons: [HTMLElement, string][] = [[this.newChatButton, 'Start New Chat'], [this.historyButton, 'Chat History...']];
    for (const [b, tooltip] of buttons) {
      ui.setDisabled(b, this.isBusy, this.isBusy ? 'Wait for response or press Stop' : null);
      if (!this.isBusy)
        ui.tooltip.bind(b, tooltip);
    }
  }
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
  private statusBar: HTMLElement = ui.divH([], 'grokky-status-bar');
  private statusLabel: HTMLElement = ui.divText('', 'grokky-status-label');
  private statusRepeats = 0;
  private followOutput = true;
  private _sessionId: string;
  private _inline: boolean = false;
  private _feature: AiFeature;
  /** Guards async {@link renderEmptyState} against out-of-order completion on rapid view switches. */
  private _emptyStateSeq = 0;
  /** Index into {@link promptHistory} while recalling prompts with the arrow keys; `null` means the live draft is shown. */
  private _promptHistoryIndex: number | null = null;
  /** The unsubmitted draft saved when the user starts cycling, restored when they cycle back past the newest entry. */
  private _promptDraft: string = '';
  private _dismissQuestion: (() => void) | null = null;

  get sessionId(): string { return this._sessionId; }


  private currentConversationId: string | null = null;
  constructor(private _contextID: string = 'global-ai-panel', view: DG.View | DG.ViewBase, opts: {inline?: boolean, feature?: AiFeature} = {}) {
    this.view = view;
    this._feature = opts.feature ?? AiFeature.Chat;
    this._sessionId = `${this._feature}-${_contextID}-${crypto.randomUUID()}`;
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
      if (event.isComposing)
        return;
      if (event.key === 'Escape' && this._promptHistoryIndex !== null) {
        this.textArea.value = this._promptDraft;
        this._promptHistoryIndex = null;
        return;
      }
      const up = event.key === 'ArrowUp';
      if ((up || event.key === 'ArrowDown') && !event.shiftKey && !event.altKey && !event.ctrlKey && !event.metaKey &&
        (this._promptHistoryIndex !== null || (up && this.textArea.selectionEnd === 0)) &&
        this.navigatePromptHistory(up ? -1 : 1))
        event.preventDefault();
      if (event.key === 'Enter' && (!event.ctrlKey && !event.metaKey)) {
        event.preventDefault();
        event.stopImmediatePropagation();
        this.handleRun();
      }
    });
    this.textArea.addEventListener('input', () => this._promptHistoryIndex = null);
    this.textArea.addEventListener('blur', (e) => {
      if (this.isBusy && (e as any).sourceCapabilities === null)
        Promise.resolve().then(() => this.textArea.focus());
    });
    ui.tooltip.bind(this.runButton, () => this.runButtonTooltip, 'left');
    this.tryAgainButton = ui.icons.sync(() => this.tryAgain(), 'Try Again');
    this.historyButton = ui.iconFA('history', () => this.showHistory());
    this.micButton = ui.iconFA('microphone', () => this.toggleSpeechRecognition(), micTooltips.default);
    if (typeof SpeechRecognition === 'undefined')
      this.setMicDisabled(micTooltips.notSupported);
    else
      this.checkMicPermission();
    this.copyConversationButton = ui.iconFA('copy', async () => {
      const success = await this.copyConversationToClipboard();
      if (success)
        grok.shell.info('Conversation copied to clipboard');
      else
        grok.shell.error('Failed to copy conversation to clipboard');
    }, 'Copy Conversation to Clipboard');

    this.newChatButton = ui.icons.add(async () => {
      if (this.isBusy)
        return;
      ui.setUpdateIndicator(this.root, true, 'Saving conversation...');
      try {
        await this.saveCurrentConversation(); // will always be first user message
      } catch (error) {
        grok.shell.error('Failed to save conversation before clearing.');
        console.error('Failed to save conversation before clearing:', error);
      } finally {
        ui.setUpdateIndicator(this.root, false);
      }
      this.resetSession();
      this.handleClear();
      this.textArea.value = '';
      this.clearAttachments();
      this.currentConversationId = null;
    });
    this.syncConversationButtons();
    this.rawRenderButton = ui.iconFA('terminal', () => {
      this._rawRender = !this._rawRender;
      this.rawRenderButton.style.color = this._rawRender ? 'var(--blue-1)' : '';
      this.root.classList.toggle('d4-ai-raw-mode', this._rawRender);
    }, 'Toggle raw console');
    this.wandButton = ui.iconFA('magic', async (e) => {
      const scopes = await resolveScopes('panel', this.view ?? grok.shell.v);
      showSuggestionsMenu(scopes, (s) => this.runSuggestion(s), e);
    }, 'Prompt suggestions');
    this.wandButton.classList.add('grokky-search-wand');
    this.setWandVisible(true);
    this.hideContentIcons();
    this.inputControlsDiv = ui.divH([
      this.wandButton, this.micButton, this.rawRenderButton,
    ], 'd4-ai-panel-input-controls');
    this.runButton.classList.add('grokky-send');
    const sessionControls = ui.divH([this.copyConversationButton, this.historyButton, this.newChatButton], 'd4-ai-panel-run-controls');
    sessionControls.style.marginLeft = 'auto';
    const messageControls = ui.divH([this.tryAgainButton, this.runButton], 'd4-ai-panel-run-controls');
    const controlsDiv = ui.divH([this.inputControlsDiv, sessionControls, ui.div([], 'd4-ribbon-separator'), messageControls], 'd4-ai-panel-controls-container');
    controlsDiv.addEventListener('mousedown', (e) => {
      if (document.activeElement === this.textArea)
        e.preventDefault();
    });
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
    this.root.appendChild(this.statusBar);
    this.root.appendChild(this.inputArea);
    this.outputArea.addEventListener('scroll', () => {
      const o = this.outputArea;
      this.followOutput = o.scrollHeight - o.clientHeight - o.scrollTop < 30;
    });

    if (this._inline)
      this.root.classList.add('d4-ai-inline-mode');
  }

  mountInto(parent: HTMLElement) {
    parent.appendChild(this.root);
  }

  get contextView(): DG.View | DG.ViewBase { return this.view; }
  setContextView(view: DG.View | DG.ViewBase): void { this.view = view; }

  activate(focus: boolean = false): void {
    this.renderEmptyState();
    if (focus)
      this.textArea.focus();
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

  protected appendFeedbackButtons(markDown: HTMLElement, onFeedback?: (helpful: boolean) => void, exec?: ExecInfo): void {
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
    const execIcons = exec?.codes.length ? this.buildExecIcons(markDown, exec) : [];
    feedbackDiv.append(copyMsg, thumbsUp, thumbsDown, ...execIcons);
    dartLike(feedbackDiv.style).set('alignItems', 'center').set('width', '100%').set('paddingBottom', '8px').set('paddingLeft', '4px');
    markDown.appendChild(feedbackDiv);
  }

  private buildExecIcons(markDown: HTMLElement, {codes, ranThisSession, view = this.view}: ExecInfo): HTMLElement[] {
    const pre = ui.element('pre', 'grokky-exec-code');
    pre.textContent = codes.map((c) => c.trim()).join('\n\n');
    const codeBlock = codeCard(pre, 'JavaScript');
    codeBlock.classList.add('grokky-hidden');
    markDown.appendChild(codeBlock);

    let state: ExecState = ranThisSession ? ExecState.ran : ExecState.idle;
    const setState = (next: ExecState) => {
      state = next;
      codeToggle.classList.remove('grokky-exec-ran', 'grokky-exec-error');
      if (next !== ExecState.idle)
        codeToggle.classList.add(`grokky-exec-${next}`);
    };
    const codeToggle = ui.iconFA('code', () => codeBlock.classList.toggle('grokky-hidden'), null);
    codeToggle.classList.add('grokky-exec-icon');
    ui.tooltip.bind(codeToggle, () => EXEC_TOOLTIPS[state]);
    setState(state);

    let rerunning = false;
    const rerunIcon = ui.iconFA('redo', async () => {
      if (rerunning) return;
      rerunning = true;
      rerunIcon.classList.add('fa-spin');
      let rerunFailed = false;
      for (const code of codes) {
        const {error} = await executeSingleBlock(code, grok.shell.v ?? view, 0);
        if (error) {
          rerunFailed = true;
          break;
        }
      }
      setState(rerunFailed ? ExecState.error : ExecState.ran);
      rerunIcon.classList.remove('fa-spin');
      rerunning = false;
    }, 'Re-run this code');
    for (const el of [codeToggle, rerunIcon])
      dartLike(el.style).set('padding', '2px').set('borderRadius', '6px');
    return [ui.div([], 'grokky-exec-sep'), codeToggle, rerunIcon];
  }

  protected appendMessage(
    aiMessage: T, uiMessage: {
      title: string, content: string, fromUser: boolean, onlyAddToMessages?: boolean, uiOnly?: boolean, messageOptions?: UIMessageOptions, execCode?: string[]
    }
  ): PanelMessageRet | undefined {
    let ret: PanelMessageRet | undefined = undefined;
    this.removeChoiceBlocks();
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
    this._uiMessages.push({fromUser: !!uiMessage.fromUser, text: uiMessage.content, title: uiMessage.title, messageOptions: uiMessage.messageOptions, execCode: uiMessage.execCode});
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
        }, uiMessage.execCode?.length ? {codes: uiMessage.execCode, ranThisSession: false} : undefined);
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
    this.outputArea.scrollTop = this.outputArea.scrollHeight;
    this.showContentIcons();
    return ret;
  }

  public startChatSession(): {session: AIPanelFuncs<T>, endSession: () => void} {
    const spinner = ui.icons.loader();
    // Speech is easy to mishear, so while the AI works in voice mode show the word that aborts the run.
    const cancelHint = ui.divText('Say "cancel" to stop', 'd4-ai-voice-cancel-hint');
    cancelHint.style.display = this.isRecognizing ? '' : 'none';
    this._voiceCancelHint = cancelHint;
    this.statusRepeats = 0;
    this.dismissQuestion();
    this.setStatus('Thinking…');
    const loader = ui.divH([spinner, this.statusLabel, cancelHint], 'd4-ai-loader');
    ui.empty(this.statusBar);
    this.statusBar.appendChild(loader);
    this.runButton.classList.remove('fal', 'fa-paper-plane');
    this.runButton.classList.add('fas', 'fa-stop');
    this.runButtonTooltip = actionButtionValues.stop;
    this.syncConversationButtons();
    return {
      session: {
        addAIMessage: (aiMessage, title, content) => this.appendMessage(aiMessage, {title: title, content: content, fromUser: false}),
        addEngineMessage: (aiMessage) => this.appendMessage(aiMessage, {title: '', content: '', fromUser: false, onlyAddToMessages: true}),
        addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => this.appendMessage('' as any, {title: '', content: msg, fromUser: fromUser, uiOnly: true, messageOptions: messageOptions}),
        markHandledNatively: () => this.markPromptHandledNatively(),
        addUserMessage: (aiMsg, content) => this.appendMessage(aiMsg, {title: '', content: content, fromUser: true}),
        addConfirmMessage: (msg?: string) => this.appendMessage('' as any, {title: '', content: '', fromUser: false, uiOnly: true, messageOptions: {confirm: {message: msg}}})!.confirmPromise,
        showInputRequest: (input: AskUserInput) => this.showInputRequest(input),
      },
      endSession: () => {
        this.runButton.classList.remove('fas', 'fa-stop');
        this.runButton.classList.add('fal', 'fa-paper-plane');
        this.runButtonTooltip = actionButtionValues.run;
        this.syncConversationButtons();
        this._voiceCancelHint = null;
        this.saveCurrentConversation().catch((e) => console.error('Failed to save conversation before hiding panel:', e));
        loader.remove();
      }
    };
  }

  get rawRender(): boolean { return this._rawRender; }
  get noPrompt(): boolean { return this._noPrompt; }

  public appendArtifact(node: HTMLElement): void {
    this.outputArea.appendChild(ui.divV([node], 'd4-ai-assistant-response-container'));
    this.showContentIcons();
  }

  /** Removes any pending inline choice blocks — they're transient and shouldn't linger once
   * the user moves on (submits a prompt, gets a new message, or opens another picker). */
  private removeChoiceBlocks(): void {
    const blocks = this.outputArea.querySelectorAll('.grokky-choice-container');
    for (let i = 0; i < blocks.length; i++)
      blocks[i].remove();
  }

  /** Renders an inline choice block that reads like an assistant reply (used by suggestion actions).
   * Picking an option removes the block; the option's handler adds its own follow-up message. */
  public addChoice(prompt: string | null, options: ChoiceOption[]): void {
    const emptyState = this.outputArea.querySelector('.grokky-empty-state');
    if (emptyState) {
      emptyState.remove();
      this.setWandVisible(true);
    }
    this.removeChoiceBlocks();
    const container = ui.divV([], 'd4-ai-assistant-response-container grokky-choice-container');
    const cards = options.map((o) => {
      const card = ui.div(ui.divText(o.label), 'grokky-choice-card');
      card.onclick = () => {
        container.remove();
        o.onSelect();
      };
      return card;
    });
    const children: HTMLElement[] = [];
    if (prompt)
      children.push(this.createStyledMarkdown(prompt));
    children.push(ui.divH(cards, 'grokky-choice-row'));
    container.appendChild(ui.divV(children, 'grokky-choice-block'));
    this.outputArea.appendChild(container);
    this.showContentIcons();
  }

  /** Appends a short assistant-style note (e.g. an action's result confirmation) to the transcript. */
  public addNote(markdown: string): void {
    this.appendMessage('' as any, {title: '', content: markdown, fromUser: false, uiOnly: true});
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

  pushNativeContext(prompt: string): void {
    this._pendingNativeContext.push(prompt);
  }

  /** Serializes the restored messages into a compact transcript the model can act on
   * ("reproduce what we did") — includes executed code blocks recorded as engine messages. */
  restoredTranscript(): string {
    const messages = this._messages[this._messages.length - 1]?.role === 'user' ?
      this._messages.slice(0, -1) : this._messages;
    const parts: string[] = [];
    for (const m of messages) {
      const c: any = (m as any).content;
      const text = typeof c === 'string' ? c :
        Array.isArray(c) ? c.map((x: any) => x?.text ?? '').filter((x: string) => x).join('\n') : '';
      if (!text.trim())
        continue;
      const role = (m as any).role === 'user' ? 'USER' : 'ASSISTANT';
      parts.push(`${role}: ${text.length > 1500 ? text.slice(0, 1500) + ' …[truncated]' : text}`);
    }
    let out = parts.join('\n');
    if (out.length > 9000)
      out = out.slice(0, 4500) + '\n[... middle of the conversation truncated ...]\n' + out.slice(-4500);
    return out ?
      '[Conversation restored from saved history — you have no memory of it. ' +
      'The transcript below is what happened earlier; treat it as this conversation\'s history. ' +
      'ASSISTANT entries starting with "[executed datagrok_exec]" are code that actually ran.]\n' +
      out : '';
  }

  flushNativeContext(): string {
    if (this._pendingNativeContext.length === 0)
      return '';
    const items = `- "${this._pendingNativeContext.join('"\n- "')}"`;
    this._pendingNativeContext = [];
    return `[The user previously asked:]\n${items}\n\n`;
  }

  resetSession(): void {
    this._sessionId = `${this._feature}-${crypto.randomUUID()}`;
    this._streamingContainer = null;
    this._streamingMarkdownEl = null;
  }

  prependViewContext(prompt: string, _view: DG.ViewBase): string {
    const ctx = buildWorkspaceContext();
    if (!ctx)
      return prompt;
    return ctx + '\n---\n\n' + prompt;
  }

  private createStreamingEl(content: string): HTMLElement | null {
    if (this._rawRender) {
      const pre = document.createElement('pre');
      pre.textContent = content;
      pre.style.cssText = 'white-space:pre-wrap;user-select:text;margin:0';
      return pre;
    }
    const md = tryMarkdown(content);
    if (md)
      dartLike(md.style).set('userSelect', 'text').set('maxWidth', '100%');
    return md;
  }

  updateStreaming(content: string): void {
    if (!this._streamingContainer) {
      this.ensureResponseBlock();
      this._streamingMarkdownEl = this.createStreamingEl(content) ?? plainTextEl(content);
      this._streamingContainer = ui.divV([this._streamingMarkdownEl], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane!.appendChild(this._streamingContainer);
    } else {
      const el = this.createStreamingEl(content);
      if (el) {
        this._streamingMarkdownEl!.replaceWith(el);
        this._streamingMarkdownEl = el;
      }
    }
    this.scrollIfFollowing();
  }

  setStatus(text: string): void {
    if (this.statusRepeats && this.statusLabel.dataset.status === text) {
      this.statusLabel.textContent = `${text} \u00d7${++this.statusRepeats}`;
      return;
    }
    this.statusRepeats = 1;
    this.statusLabel.dataset.status = text;
    this.statusLabel.textContent = text;
    this.statusLabel.classList.remove('grokky-status-in');
    void this.statusLabel.offsetWidth;
    this.statusLabel.classList.add('grokky-status-in');
  }

  showError(text: string): void {
    const error = ui.divText(text, 'grokky-status-error');
    ui.tooltip.bind(error, text);
    ui.empty(this.statusBar);
    this.statusBar.appendChild(error);
  }

  private scrollIfFollowing(): void {
    if (this.followOutput)
      this.outputArea.scrollTop = this.outputArea.scrollHeight;
  }

  async finalizeStreaming(displayContent: string, execCodes: string[], view: DG.ViewBase): Promise<void> {
    if (this._rawRender) {
      this._streamingMarkdownEl = null;
      this._streamingContainer = null;
      this._uiMessages.push({fromUser: false, text: displayContent, execCode: execCodes.length ? execCodes : undefined, messageOptions: {finalResult: displayContent}});
      return;
    }
    this.renderFinalContent(displayContent, execCodes, view);
  }

  public appendStreamedElement(el: HTMLElement): void {
    this.ensureResponseBlock();
    this._aiMessagesAccordionPane!.appendChild(ui.divV([el], 'd4-ai-assistant-response-container'));
    this._streamingContainer = null;
    this._streamingMarkdownEl = null;
    this.scrollIfFollowing();
  }

  public appendUiMessage(content: string): void {
    this.appendMessage('' as any, {title: '', fromUser: false, uiOnly: true, content, messageOptions: {system: true}});
  }

  protected renderFinalContent(content: string, execCodes: string[], view: DG.ViewBase): void {
    const markDown = this.createStyledMarkdown(content);
    this.appendFeedbackButtons(markDown, undefined,
      execCodes.length ? {codes: execCodes, ranThisSession: true, view} : undefined);

    if (this._streamingMarkdownEl) {
      this._streamingMarkdownEl.replaceWith(markDown);
      this._streamingMarkdownEl = null;
      this._streamingContainer = null;
    } else if (content || execCodes.length) {
      this.ensureResponseBlock();
      this._aiMessagesAccordionPane!.appendChild(ui.divV([markDown], 'd4-ai-assistant-response-container'));
    }

    this._uiMessages.push({fromUser: false, text: content, execCode: execCodes.length ? execCodes : undefined, messageOptions: {finalResult: content}});
  }

  clearStreaming(): void {
    this._streamingContainer?.remove();
    this._streamingContainer = null;
    this._streamingMarkdownEl = null;
  }

  showInputRequest(input: AskUserInput): void {
    const questions = input.questions ?? [];
    const answers = questions.map((q) => q.options[0]?.label ?? '');
    const sections = questions.map((q, i): QuestionSection => {
      const pills: HTMLElement[] = q.options.map((o) => {
        const pill = ui.divText(o.label, 'grokky-choice-card');
        if (o.description)
          ui.tooltip.bind(pill, o.description);
        pill.onclick = () => {
          answers[i] = o.label;
          pills.forEach((p) => p.classList.toggle('grokky-choice-selected', p === pill));
        };
        return pill;
      });
      pills[0]?.classList.add('grokky-choice-selected');
      return {
        strip: ui.divH([ui.divText(q.question)], 'grokky-strip grokky-strip-question'),
        body: ui.div([ui.divH(pills, 'grokky-choice-row')], 'grokky-confirm-body'),
      };
    });
    const close = (style: StripStyle, icon: string, answerText: (questionIdx: number) => string): void => {
      this._dismissQuestion = null;
      footer.remove();
      sections.forEach((s, i) => {
        s.body.remove();
        s.strip.classList.replace('grokky-strip-question', `grokky-strip-${style}`);
        s.strip.append(ui.divH([ui.iconFA(icon), ui.divText(answerText(i))], 'grokky-strip-answer'));
      });
    };
    this.dismissQuestion();
    this._dismissQuestion = () => close('muted', 'minus', () => 'Not answered');
    const footer = ui.divH([ui.bigButton('Submit', () => {
      if (this.isBusy)
        return;
      close('success', 'check', (i) => answers[i]);
      const prompt = questions.length === 1 ? answers[0] :
        questions.map((q, i) => `${q.question} ${answers[i]}`).join('\n');
      this._onRunRequest.next({prevMessages: this._messages, currentPrompt: {prompt, isAnswer: true} as K});
      this.textArea.focus();
    })], 'grokky-question-footer');

    const card = ui.divV([...sections.flatMap((s) => [s.strip, s.body]), footer],
      'grokky-confirm-card grokky-question-card');
    this.ensureResponseBlock();
    this._aiMessagesAccordionPane!.appendChild(ui.divV([card], 'd4-ai-assistant-response-container'));
    this.outputArea.scrollTop = this.outputArea.scrollHeight;
  }

  dismissQuestion(): void {
    this._dismissQuestion?.();
  }

  private tryAgain() {
    if (this._messages.length === 0)
      return; // should never happen, but just in case
    // The prompt itself is resent: a request stopped before the agent picked it up is not in its history.
    const last = this.promptHistory.pop();
    if (!last || this.isBusy)
      return;
    const inputs = this.getCurrentInputs();
    inputs.prompt = last;
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

  /** Replaces the input with the previous (-1) or next (+1) prompt from this session's history; false if there is none. */
  private navigatePromptHistory(direction: -1 | 1): boolean {
    const hist = this.promptHistory;
    const current = this._promptHistoryIndex ?? hist.length;
    const next = current + direction;
    if (next < 0 || next > hist.length)
      return false;
    if (this._promptHistoryIndex === null)
      this._promptDraft = this.textArea.value;
    this._promptHistoryIndex = next === hist.length ? null : next;
    this.textArea.value = next === hist.length ? this._promptDraft : hist[next];
    this.textArea.selectionStart = this.textArea.selectionEnd = this.textArea.value.length;
    return true;
  }

  protected handleRun() {
    if (this.isBusy)
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
    this._pendingNativeContext = [];
    this._promptHistoryIndex = null;
    this._lastUserPromptContainer = null;
    this.outputArea.innerHTML = '';
    this._onClearChatRequest.next();
    this.hideContentIcons();
    this.renderEmptyState();
  }

  protected runSuggestion(s: Suggestion): void {
    if (!runSuggestionAction(s.action, this)) {
      this.textArea.value = s.prompt ?? '';
      this.handleRun();
    }
    this.textArea.focus();
  }

  /** The view the panel currently works against: the pinned one for owned panels, else the live current view. */
  protected get liveView(): DG.View | DG.ViewBase | null {
    return this.view ?? grok.shell.v;
  }

  protected shouldShowEmptyState(): boolean {
    // Global suggestions (Anywhere, Code generation) apply on any view, so show the empty state
    // and wand everywhere — not just on table views (view/column blocks self-add when applicable).
    return this._uiMessages.length === 0;
  }

  private setWandVisible(visible: boolean): void {
    this.wandButton.style.display = visible ? '' : 'none';
  }

  /** Rebuilds the suggestion cards for the current view; removes them when they no longer apply. */
  protected async renderEmptyState(): Promise<void> {
    const seq = ++this._emptyStateSeq;
    const removeExisting = () => this.outputArea.querySelector('.grokky-empty-state')?.remove();
    if (!this.shouldShowEmptyState()) {
      removeExisting();
      this.setWandVisible(true);
      return;
    }
    const scopes = await resolveScopes('panel', this.liveView);
    if (seq !== this._emptyStateSeq || !this.shouldShowEmptyState())
      return;
    removeExisting();

    const accordion = ui.accordion();
    const defaultOpenCategories = 1;
    scopes.forEach((s, i) => this.addSuggestionCategory(accordion, s, i < defaultOpenCategories));

    const root = ui.panel([ui.h2('I\'m Grokky, your AI assistant. What can I help you with?'), accordion.root], 'grokky-empty-state');
    this.outputArea.appendChild(root);
    this.setWandVisible(false);
  }

  private addSuggestionCategory(accordion: DG.Accordion, s: Block, expanded: boolean): void {
    const pane = accordion.addPane(s.label, () => {
      const rows = s.suggestions.map((sg) => {
        const row = ui.divText(sg.label ?? sg.prompt ?? '', 'grokky-suggestion-item');
        row.onclick = () => this.runSuggestion(sg);
        if (sg.prompt && sg.label && sg.label !== sg.prompt)
          ui.tooltip.bind(row, () => sg.prompt ?? '');
        return row;
      });
      return ui.divV(rows);
    }, expanded);

    const icon = ui.iconFA(s.icon ?? 'circle');
    icon.classList.add('grokky-scope-icon');
    if (s.color)
      icon.style.color = s.color;
    pane.root.querySelector('.d4-accordion-pane-header')?.prepend(icon);
  }

  private async showHistory() {
    if (this.isBusy)
      return;
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
      const deleting = new Map<string, number>();
      const finish = (id: string, remove: boolean) => {
        clearTimeout(deleting.get(id));
        deleting.delete(id);
        if (remove)
          ConversationStorage.deleteConversation(id).catch(() => grok.shell.error('Failed to delete chat'));
      };
      for (const conv of conversations) {
        if (this.currentConversationId && conv.id === this.currentConversationId)
          continue;
        const date = new Date(conv.timestamp).toLocaleString();
        const preview = conv.initialPrompt.substring(0, 40).replace(/\|/g, '/'); // '|' would make a submenu
        const displayText = `${preview}${conv.initialPrompt.length > 40 ? '...' : ''} (${date})`;

        menu.item(displayText, () => {
          finish(conv.id, false);
          this.loadConversation(conv.id);
        });
        const row = menu.root.querySelector('.d4-menu-item:last-child')!;
        const toggle = (e: MouseEvent) => {
          e.stopPropagation();
          finish(conv.id, false);
          if (row.classList.toggle('grokky-history-deleted')) {
            deleting.set(conv.id, window.setTimeout(() => {
              row.remove();
              finish(conv.id, true);
            }, 5000));
          }
        };
        const trash = ui.icons.delete(toggle, 'Delete chat');
        const undo = ui.icons.undo(toggle, 'Undo');
        trash.classList.add('grokky-history-icon', 'grokky-history-delete');
        undo.classList.add('grokky-history-icon', 'grokky-history-undo', 'grok-icon-blue');
        row.append(trash, undo);
      }
      menu.onClose.subscribe(() => Array.from(deleting.keys()).forEach((id) => finish(id, true)));

      menu.separator();
      menu.item('Clear All History', async () => {
        // if (await grok.shell.('Delete all conversation history?')) {
        //   await ConversationStorage.clearAll();
        //   grok.shell.info('History cleared');
        // }
        ui.dialog('Delete all conversation history').add(ui.divText('This action will permanently delete all saved conversations. Are you sure you want to proceed?'))
          .onOK(async () => {
            await ConversationStorage.clearAll();
            this.currentConversationId = null;
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
        this.appendMessage(null as any, {title: msg.title ?? '', content: msg.text, fromUser: msg.fromUser, uiOnly: true, messageOptions: msg.messageOptions, execCode: msg.execCode}); // no loader
      });
      this.resetSession();
      // Chat history is per browser, so only the user who saved the chat resumes its runtime session.
      const meta = conv.meta;
      if (meta?.sessionId && meta.userId === grok.shell.user.id) {
        this._sessionId = meta.sessionId;
        ClaudeRuntimeClient.getInstance().expectResume(meta.sessionId);
      }
      this.afterConversationLoad(conv);
      this.textArea.focus();
      //grok.shell.info(`Loaded conversation: ${conv.initialPrompt.substring(0, 50)}...`);
    } catch (error) {
      console.error('Failed to load conversation:', error);
      grok.shell.error('Failed to load conversation');
    }
  }

  protected afterConversationLoad(conversation: StoredConversationWithContext<T>) {
    // do nothing in base class
  }

  protected getConversationMeta(): {sessionId: string, userId: string} {
    return {sessionId: this._sessionId, userId: grok.shell.user.id};
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
    if (this.micAccessDenied)
      return;
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

  private async checkMicPermission(): Promise<void> {
    if (!navigator.permissions?.query)
      return;

    try {
      const status = await navigator.permissions.query({name: 'microphone' as PermissionName});
      this.applyMicPermissionStatus(status.state);
      status.onchange = () => this.applyMicPermissionStatus(status.state);
    } catch (error) {
      console.error('Failed to query microphone permission:', error);
    }
  }

  private applyMicPermissionStatus(state: PermissionState): void {
    this.setMicDisabled(state === 'denied' ? micTooltips.accessDenied : null);
  }

  private setMicDisabled(reason: string | null): void {
    this.micAccessDenied = reason != null;
    ui.setDisabled(this.micButton, this.micAccessDenied, reason ?? micTooltips.default);
  }

  private startRecognition() {
    this.isRecognizing = true;
    this.voiceTranscript = '';
    this.micButton.classList.remove('fa-microphone');
    this.micButton.classList.add('fa-stop');
    this.micButton.style.color = 'orangered';
    this.syncVoiceCancelHint();
    this.listen();
  }

  /** Starts a recognition session. Chrome ends a continuous session on its own after a while or on
   * a network hiccup; `onend` then calls this again, keeping the phrase heard so far. */
  private listen() {
    try {
      // Continuous mode, so a pause between words does not end the phrase: Chrome's own end-of-speech
      // detection cut it after the first word for mics that gate the signal to silence between words.
      // The phrase is submitted after VOICE_SUBMIT_DELAY_MS without new speech instead.
      this.recognition = new SpeechRecognition();
      this.recognition.continuous = true;
      this.recognition.interimResults = true;
      this.recognition.lang = 'en-US';
      if (!this.voiceTranscript)
        ui.setUpdateIndicator(this.textAreaDiv, true, 'Listening...');

      this.recognition.onresult = (event) => {
        ui.setUpdateIndicator(this.textAreaDiv, false);
        let interim = '';
        for (let i = event.resultIndex; i < event.results.length; i++) {
          const text = event.results[i][0].transcript;
          if (!event.results[i].isFinal) {
            interim += text;
            continue;
          }
          const command = text.trim().replace(/[.,;:!?]+$/, '').toLowerCase();
          if (command === 'stop' || command === 'cancel') {
            this.onVoiceCancel();
            return;
          }
          this.voiceTranscript += text;
        }

        // Speech while the AI works is not a prompt (one can't be submitted yet); only "stop"/"cancel"
        // counts then, to abort the run.
        if (this.isBusy) {
          this.voiceTranscript = '';
          return;
        }

        this._promptHistoryIndex = null;
        this.textArea.value = (this.voiceTranscript + interim).trim();
        window.clearTimeout(this.voiceSubmitTimer);
        this.voiceSubmitTimer = window.setTimeout(() => this.submitVoicePrompt(), VOICE_SUBMIT_DELAY_MS);
      };

      this.recognition.onerror = (event) => {
        if (event.error === 'no-speech')
          return; // silence timeout — onend restarts the mic, no user action needed
        if (event.error !== 'aborted')
          console.error('Speech recognition error:', event.error);
        ui.setUpdateIndicator(this.textAreaDiv, false);

        let errorMessage = 'Speech recognition error';
        switch (event.error) {
        case 'audio-capture':
          errorMessage = micTooltips.noDevice;
          this.setMicDisabled(micTooltips.noDevice);
          break;
        case 'not-allowed':
          errorMessage = micTooltips.accessDenied;
          this.setMicDisabled(micTooltips.accessDenied);
          break;
        case 'network':
          errorMessage = 'Network error during speech recognition';
          break;
        case 'aborted':
          // Something else took the microphone, e.g. voice input in another tab: no error message needed
          break;
        default:
          errorMessage = `Speech recognition error: ${event.error}`;
        }

        if (event.error !== 'aborted')
          grok.shell.error(errorMessage);


        this.stopRecognition();
      };

      this.recognition.onend = () => {
        if (this.isRecognizing)
          this.listen();
      };

      this.recognition.start();
    } catch (error) {
      console.error('Failed to start speech recognition:', error);
      grok.shell.error('Failed to start speech recognition');
      this.stopRecognition();
    }
  }

  /** "stop"/"cancel" aborts the AI run while it works, and otherwise drops the phrase and the mic. */
  private onVoiceCancel() {
    window.clearTimeout(this.voiceSubmitTimer);
    this.voiceTranscript = '';
    if (this.isBusy) {
      this.terminate();
      return;
    }
    this._promptHistoryIndex = null;
    this.textArea.value = '';
    this.stopRecognition();
    this.textArea.focus();
  }

  private submitVoicePrompt() {
    this.voiceTranscript = '';
    this.handleRun();
  }

  /** Stops listening; a phrase heard so far stays in the prompt box, unsent, to edit or send by hand. */
  private stopRecognition() {
    this.isRecognizing = false;
    this.voiceTranscript = '';
    window.clearTimeout(this.voiceSubmitTimer);
    this.syncVoiceCancelHint();

    if (this.recognition) {
      // A stopped session can still report late; it must not log or switch off a newer one.
      this.recognition.onresult = this.recognition.onerror = this.recognition.onend = null;
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
