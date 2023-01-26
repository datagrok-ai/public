/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// external modules dependencies
import $ from 'cash-dom';

// inner dependencies
import '../../css/colored-text-input.css';


/** Class for colorizing input in the textarea of DG.InputBase.  */
export class ColoredTextInput {
  constructor(
    textInputBase: DG.InputBase<string>,
    painter: (str: string) => HTMLSpanElement[],
    /** Resize, no scrolls  */
    resizeable: boolean = true
  ) {
    this.textInputBase = textInputBase;
    this.painter = painter;
    $(this.root).addClass('colored-text-input');
    if (resizeable) {
      // make input field automatically resizeable
      this.textInputBase.onInput(
        () => {
          // necessary for the field to be squeezable, not only expandable
          $(this.textArea).css('height', 0);
          $(this.textArea).css('height', (this.textArea!.scrollHeight) + 'px');
        }
      );
    }
    this.highlights = ui.div([]);
    this.root.appendChild(this.highlights);

    this.textInputBase.onInput(() => this.colorize());
  }

  private textInputBase: DG.InputBase<string>;
  private highlights: HTMLDivElement;
  /** Divide input value into an array of spans, each with its own text color, use -webkit-text-fill-color. */
  private painter: (str: string) => HTMLSpanElement[];

  get textArea() {
    return this.textInputBase.root.getElementsByTagName('textarea').item(0);
  };

  get root() { return this.textInputBase.root; };

  private colorize() {
    const spans = this.painter(this.textInputBase.value);
    this.highlights.innerHTML = '';
    spans.forEach((span: HTMLSpanElement) => this.highlights.appendChild(span));
  }
}
