import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import '../../css/forge.css';

/** A block that shows or hides its content under a header with a chevron, a caption and an optional summary
 * (the input categories of Diff Studio). The state is not remembered. */
export class CollapsibleGroup {
  readonly root: HTMLDivElement;
  readonly body: HTMLDivElement;
  private readonly expandIcon: HTMLElement;
  private readonly collapseIcon: HTMLElement;
  private readonly summary = ui.div([], 'forge-group-summary');
  private readonly onExpand: (() => void) | undefined;

  /** [onExpand] runs when the group goes from collapsed to expanded. */
  constructor(caption: string, content: HTMLElement[], isExpanded: boolean = true, onExpand?: () => void) {
    this.expandIcon = ui.iconFA('chevron-right', null, `Show ${caption}`);
    this.collapseIcon = ui.iconFA('chevron-down', null, `Hide ${caption}`);
    // The two chevrons differ in width: a fixed box keeps every caption at the same x in either state.
    const chevron = ui.div([this.expandIcon, this.collapseIcon], 'forge-group-chevron');
    const header = ui.divH([chevron, ui.label(caption), this.summary], 'forge-group-header');
    header.addEventListener('click', () => this.setExpanded(!this.isExpanded));
    this.body = ui.div(content, 'forge-group-body');
    this.root = ui.div([header, this.body], 'forge-group');
    this.setExpanded(isExpanded);
    this.onExpand = onExpand;
  }

  get isExpanded(): boolean {
    return this.body.style.display !== 'none';
  }

  setExpanded(isExpanded: boolean): void {
    const isOpening = isExpanded && !this.isExpanded;
    ui.setDisplay(this.body, isExpanded);
    ui.setDisplay(this.expandIcon, !isExpanded);
    ui.setDisplay(this.collapseIcon, isExpanded);
    if (isOpening)
      this.onExpand?.();
  }

  /** The text after the caption, in the error color when [isInvalid]. */
  setSummary(text: string, isInvalid: boolean): void {
    this.summary.textContent = text;
    this.summary.classList.toggle('forge-group-invalid', isInvalid);
  }

  /** Expands the group whenever one of [inputs] is validated with an error; the caller keeps the subscriptions. */
  expandOnError(inputs: DG.InputBase[]): rxjs.Subscription[] {
    return inputs.map((input) => input.onValidated.subscribe((messages) => {
      if (messages.length > 0)
        this.setExpanded(true);
    }));
  }
}
