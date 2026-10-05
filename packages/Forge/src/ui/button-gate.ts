import * as ui from 'datagrok-api/ui';

/** Keeps [button] disabled while [reason] returns a text, shown as the button's tooltip. */
export class ButtonGate {
  private readonly button: HTMLElement;
  private readonly reason: () => string | null;
  private hasTooltip = false;

  constructor(button: HTMLElement, reason: () => string | null) {
    this.button = button;
    this.reason = reason;
  }

  /** Reads the reason again. A disabled button shows its tooltip through an overlay, built once per disabled period
   * on the attached button. */
  update(): void {
    const isBlocked = this.reason() !== null;
    if (isBlocked && !this.hasTooltip)
      ui.setDisabled(this.button, true, this.reason);
    else if (!isBlocked)
      ui.setDisabled(this.button, false);
    this.hasTooltip = isBlocked;
  }
}
