import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ClaudeRuntimeClient} from './claude-runtime-client';
import {AIPanel} from '../llm-utils/panel';
import {fireAIPanelToggleEvent, dartLike} from '../utils';

class ClaudeAIPanel extends AIPanel {
  protected get placeHolder() { return 'Ask Claude about your data...'; }
  private _streamingContainer: HTMLElement | null = null;
  private _streamingMarkdownEl: HTMLElement | null = null;

  constructor(view: DG.TableView) {
    super(`claude-${view.dataFrame?.name ?? view.name ?? 'view'}`, view);
  }

  /** Create or update the streaming response using the same visual structure as appendMessage. */
  updateStreaming(content: string, loader: HTMLElement): void {
    if (!this._streamingContainer) {
      loader.style.display = 'none';
      // Create the same accordion + response structure as appendMessage()
      if (!this._aiMessagesAccordionPane) {
        const acord = ui.accordion();
        this._aiMessagesAccordionPane = ui.divV([], 'd4-ai-messages-accordion-pane');
        const pane = acord.addPane('Responses', () => this._aiMessagesAccordionPane!, true, undefined, false);
        pane.expanded = true;
        acord.root.style.width = 'calc(100% - 35px)';
        this.outputArea.appendChild(acord.root);
      }
      this._streamingMarkdownEl = ui.markdown(content);
      dartLike(this._streamingMarkdownEl.style).set('userSelect', 'text').set('maxWidth', '100%');
      this._streamingContainer = ui.divV([this._streamingMarkdownEl], 'd4-ai-assistant-response-container');
      this._aiMessagesAccordionPane.appendChild(this._streamingContainer);
    } else {
      const md = ui.markdown(content);
      dartLike(md.style).set('userSelect', 'text').set('maxWidth', '100%');
      this._streamingMarkdownEl!.replaceWith(md);
      this._streamingMarkdownEl = md;
    }
    this.outputArea.scrollTop = this.outputArea.scrollHeight;
  }

  /** Finalize the streaming response: update with final content, add copy/feedback buttons, register in history. */
  finalizeStreaming(content: string): void {
    if (!this._streamingContainer || !this._streamingMarkdownEl)
      return;

    // Final markdown render
    const markDown = ui.markdown(content);
    markDown.style.position = 'relative';
    dartLike(markDown.style).set('userSelect', 'text').set('maxWidth', '100%');

    // Add copy button for code blocks
    if (markDown.querySelector('pre > code')) {
      const copyButton = ui.icons.copy(() => {}, 'Copy Code');
      copyButton.classList.add('d4-ai-copy-code-button');
      markDown.appendChild(copyButton);
      copyButton.addEventListener('click', () => {
        const codeElement = markDown.querySelector('pre > code');
        if (codeElement) {
          const header = markDown.children[0];
          if (header && header.tagName?.toLowerCase() !== 'pre')
            (header as HTMLElement).style.marginRight = '16px';
          navigator.clipboard.writeText(codeElement.textContent || '').then(() => {
            copyButton.classList.add('d4-ai-copy-code-button-copied');
            setTimeout(() => copyButton.classList.remove('d4-ai-copy-code-button-copied'), 600);
          }).catch(() => grok.shell.error('Failed to copy code to clipboard.'));
        }
      });
    }

    // Add feedback buttons
    const feedbackDiv = ui.divH([], 'd4-ai-panel-feedback-div');
    const thumbsUp = ui.iconFA('thumbs-up', () => {
      grok.shell.info('Thanks for your feedback!');
      handleFeedback(true);
    }, 'Helpful');
    const thumbsDown = ui.iconFA('thumbs-down', () => {
      grok.shell.info('Thanks for your feedback!');
      handleFeedback(false);
    }, 'Not Helpful');
    [thumbsUp, thumbsDown].forEach((el) => dartLike(el.style).set('padding', '2px').set('borderRadius', '6px'));

    function handleFeedback(helpful: boolean) {
      [thumbsUp, thumbsDown].forEach((el) => { el.style.backgroundColor = ''; });
      (helpful ? thumbsUp : thumbsDown).style.backgroundColor =
        helpful ? 'rgba(0, 150, 30, 0.2)' : 'rgba(200, 0, 0, 0.2)';
    }
    feedbackDiv.appendChild(thumbsUp);
    feedbackDiv.appendChild(thumbsDown);
    dartLike(feedbackDiv.style).set('gap', '8px').set('alignItems', 'center').set('width', '100%').set('paddingBottom', '8px').set('paddingLeft', '4px');
    markDown.appendChild(feedbackDiv);

    this._streamingMarkdownEl.replaceWith(markDown);
    this._streamingMarkdownEl = null;
    this._streamingContainer = null;

    // Register in UI messages for history/copy support
    this._uiMessages.push({fromUser: false, text: content, messageOptions: {finalResult: content}});
  }

  /** Remove the streaming element on error/cancellation. */
  clearStreaming(): void {
    this._streamingContainer?.remove();
    this._streamingContainer = null;
    this._streamingMarkdownEl = null;
  }
}

export async function setupClaudeAIPanelUI() {
  if (!grok.ai.config.configured)
    return;

  const handleView = (tableView: DG.TableView) => {
    const iconFse = ui.iconSvg('ai.svg', () => fireAIPanelToggleEvent(tableView), 'Ask Claude \n Ctrl+Shift+I');
    iconFse.style.width = iconFse.style.height = '18px';
    tableView.setRibbonPanels([...tableView.getRibbonPanels(), [iconFse]]);

    const panel = new ClaudeAIPanel(tableView);
    panel.hide();

    panel.onRunRequest.subscribe(async (args) => {
      const {session, endSession, loader} = panel.startChatSession();
      const sessionId = `claude-${tableView.dataFrame?.name ?? 'view'}-${crypto.randomUUID()}`;
      let accumulated = '';
      const subs: {unsubscribe: () => void}[] = [];

      const cleanup = () => {
        subs.forEach((s) => s.unsubscribe());
        subs.length = 0;
      };

      try {
        const client = ClaudeRuntimeClient.getInstance();
        if (!client.connected)
          await client.connect(tableView);
        else
          client.updateContext(tableView);

        // Add user message to UI (uiOnly — Claude SDK manages its own message history)
        session.addUiMessage(args.currentPrompt.prompt, true);

        // Subscribe to streaming chunks — render in proper accordion structure
        subs.push(client.onChunk.subscribe((evt) => {
          if (evt.sessionId !== sessionId) return;
          accumulated += evt.content;
          panel.updateStreaming(accumulated, loader);
        }));

        // Subscribe to tool use notifications
        subs.push(client.onToolUse.subscribe((evt) => {
          if (evt.sessionId !== sessionId) return;
          const statusText = evt.status === 'completed' ? 'completed' : 'running';
          panel.updateStreaming(accumulated + `\n\n*Tool: ${evt.tool}* (${statusText})`, loader);
        }));

        // Subscribe to final message — finalize in place, no flash
        subs.push(client.onFinal.subscribe((evt) => {
          if (evt.sessionId !== sessionId) return;
          panel.finalizeStreaming(evt.content);
          endSession();
          cleanup();
        }));

        // Subscribe to errors
        subs.push(client.onError.subscribe((evt) => {
          if (evt.sessionId !== sessionId) return;
          panel.clearStreaming();
          grok.shell.error(`Claude: ${evt.message}`);
          endSession();
          cleanup();
        }));

        // Determine context and send
        const context: {viewType?: string, connectionId?: string} = {
          viewType: DG.VIEW_TYPE.TABLE_VIEW,
        };
        client.send(sessionId, args.currentPrompt.prompt, context);
      }
      catch (e: any) {
        panel.clearStreaming();
        grok.shell.error(`Claude runtime: ${e.message}`);
        console.error('Claude runtime error:', e);
        endSession();
        cleanup();
      }
    });
  };

  // Handle already opened views
  Array.from(grok.shell.tableViews).filter((v) => v.dataFrame != null).forEach((view) => {
    handleView(view as DG.TableView);
  });

  grok.events.onViewAdded.subscribe((view) => {
    if (view.type === DG.VIEW_TYPE.TABLE_VIEW)
      handleView(view as DG.TableView);
  });
}
