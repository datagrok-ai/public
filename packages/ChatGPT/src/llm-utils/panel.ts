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

export type ModelOption = 'Fast' | 'Deep Research';
export const ModelType: {[type in ModelOption]: ChatModel} = {
  Fast: 'gpt-4o-mini',
  //   ['Deep Research']: 'gpt-5.1', // hell of a smart model but expensive
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
  stop: 'Click to stop AI Generation',
} as const;

export class AIPanel<T extends MessageType = OpenAI.Chat.ChatCompletionMessageParam, K extends AIPanelInputs = AIPanelInputs> {
  private root: HTMLElement;
  private inputArea: HTMLElement;
  private outputArea: HTMLElement;
  private textArea: HTMLTextAreaElement;
  private runButton: HTMLElement;
  private modelInput: DG.InputBase<ModelOption>;
  private clearButton: HTMLElement;
  private tryAgainButton: HTMLElement;
  private _onRunRequest = new rxjs.Subject<{prevMessages: T[], currentPrompt: K}>();
  private _messages: T[] = [];
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
  constructor() {
    this.root = ui.divV([], 'd4-ai-generation-panel');
    this.inputArea = ui.divV([], 'd4-ai-panel-input-area');
    this.outputArea = ui.divV([], 'd4-ai-panel-output-area');
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
    ui.tooltip.bind(this.runButton, () => this.runButtonTooltip, 'top');
    this.tryAgainButton = ui.icons.sync(() => this.tryAgain(), 'Try Again');

    this.clearButton = ui.icons.delete(() => this.handleClear(), 'Clear Chat');
    this.modelInput = ui.input.choice('Model', {
      items: Object.keys(ModelType) as ModelOption[],
      value: 'Deep Research',
      nullable: false,
      tooltipText: 'Select AI model to use',
    }) as DG.InputBase<ModelOption>;
    this.hideContentIcons();
    this.inputControlsDiv = ui.divH([
      this.modelInput.input
    ], 'd4-ai-panel-input-controls');
    this.runButton.style.color = 'var(--blue-1)';
    const inputControlsRight = ui.divH([this.clearButton, this.tryAgainButton, this.runButton], 'd4-ai-panel-run-controls');
    inputControlsRight.style.marginLeft = 'auto';
    const controlsDiv = ui.divH([this.inputControlsDiv, inputControlsRight], 'd4-ai-panel-controls-container');
    this.inputArea.appendChild(this.textArea);
    this.inputArea.appendChild(controlsDiv);
    this.root.appendChild(this.outputArea);
    this.root.appendChild(this.inputArea);
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
    AIPanel._lastDockedPanel = grok.shell.dockManager.dock(this.root, 'right', null, 'Chat', 0.25);
  }

  hide() {
    grok.shell.dockManager.close(this.root);
  }

  toggle() {
    document.contains(this.root) ? this.hide() : this.show();
  }

  dispose() {
    this.root.remove();
    this._messages = [];
  }


  public getCurrentInputs(): K {
    return {
      prompt: this.textArea.value,
      model: this.modelInput.value!,
    } as K;
  }

  private _aiMessagesAccordionPane: HTMLElement | null = null;
  protected appendMessage(aiMessage: T, uiMessage: {title: string, content: string, fromUser: boolean, onlyAddToMessages?: boolean, uiOnly?: boolean}, loader: HTMLElement) {
    if (!uiMessage.uiOnly)
      this._messages.push(aiMessage);
    if (uiMessage.onlyAddToMessages)
      return;
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
      dartLike(markDown.style).set('userSelect', 'text').set('maxWidth', '100%');
      const titleText = uiMessage.title ? [ui.h3(uiMessage.title, 'd4-ai-assistant-response-title')] : [];
      const assistantDiv = ui.divV([...titleText, markDown], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane.appendChild(assistantDiv);
    }
    // reapend loader to the end
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
      addUIMessage: (content: string, fromUser: boolean) => {
        this.appendMessage('' as any, {title: '', content: content, fromUser: fromUser, uiOnly: true}, loader);
      },
      endSession: () => {
        this.runButton.classList.remove('fas', 'fa-stop');
        this.runButton.classList.add('fal', 'fa-paper-plane');
        this.runButton.style.color = 'var(--blue-1)';
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
    this.clearButton.style.display = 'flex';
    this.tryAgainButton.style.display = 'flex';
  }

  protected hideContentIcons() {
    this.clearButton.style.display = 'none';
    this.tryAgainButton.style.display = 'none';
  }

  private handleClear() {
    this._messages = [];
    this.outputArea.innerHTML = '';
    this._onClearChatRequest.next();
    this.hideContentIcons();
  }
}

export class DBAIPanel<T extends MessageType = OpenAI.Chat.ChatCompletionMessageParam> extends AIPanel<T, DBAIPanelInputs> {
  protected get placeHolder() { return 'Type your database query prompt here...'; }
  protected schemaInput: DG.InputBase<string>;
  constructor(schemas: string[], defaultSchema: string) {
    super();
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
