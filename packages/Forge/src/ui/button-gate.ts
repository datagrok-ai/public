import * as ui from 'datagrok-api/ui';

/** Keeps [button] disabled while [reason] returns a text, shown as the button's tooltip; [enabledTooltip] is its
 * tooltip otherwise (none without it). */
export class ButtonGate {
  private readonly button: HTMLElement;
  private readonly reason: () => string | null;
  private readonly enabledTooltip: string | undefined;
  private readonly reasonText = ui.divText('');
  private shown: string | null = null;
  private overlay: Element | null = null;

  constructor(button: HTMLElement, reason: () => string | null, enabledTooltip?: string) {
    this.button = button;
    this.reason = reason;
    this.enabledTooltip = enabledTooltip;
    if (enabledTooltip !== undefined)
      ui.tooltip.bind(button, enabledTooltip);
  }

  /** Reads the reason again. A disabled button shows its tooltip through a body-level overlay that lives until the
   * button is enabled or out of the document; its content is [reasonText], so a new reason only changes the text and
   * the overlay is built again only when the last one is gone. */
  update(): void {
    const reason = this.reason();
    if (reason !== null && reason !== this.shown) {
      this.reasonText.textContent = reason;
      if (this.overlay?.isConnected)
        ui.setDisabled(this.button, true);
      else {
        // The overlay shows an element as it is, but a function as its source: the element is typed as text for it.
        ui.setDisabled(this.button, true, this.reasonText as unknown as string);
        this.overlay = document.querySelector('.d4-tooltip-overlays')?.lastElementChild ?? null;
      }
    } else if (reason === null && this.shown !== null) {
      ui.setDisabled(this.button, false);
      // The disabled reason is the button's own tooltip too: the last binding wins.
      ui.tooltip.bind(this.button, () => this.enabledTooltip ?? null);
    }
    this.shown = reason;
  }
}
