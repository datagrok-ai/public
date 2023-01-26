/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/rich-text-input.css';

import $ from 'cash-dom';

export class RichTextInput {
  constructor(
    textInputBase: DG.InputBase<string>,
    colorizer: (str: string) => HTMLSpanElement[],
    resizeable: boolean = true
  ) {
    this.textInputBase = textInputBase;
    this.colorizer = colorizer;
    $(this.root).addClass('rich-text');
    if (resizeable) {
      // make input field automatically resizeable
      this.textInputBase.onInput(
        () => {
          $(this.textArea).css('height', 0);
          $(this.textArea).css('height', (this.textArea!.scrollHeight) + 'px');
        }
      );
    }
    this.highlights = ui.div([] );
    this.root.appendChild(this.highlights);

    this.textInputBase.onInput(() => this.colorize());
  }

  private textInputBase: DG.InputBase<string>;
  private highlights: HTMLDivElement;
  private colorizer: (str: string) => HTMLSpanElement[];

  get textArea() {
    return this.textInputBase.root.getElementsByTagName('textarea').item(0);
  };

  get root() { return this.textInputBase.root; };

  private colorize() {
    const spans = this.colorizer(this.textInputBase.value);
    this.highlights.innerHTML = '';
    spans.forEach((span: HTMLSpanElement) => this.highlights.appendChild(span));
  }
}

export function demoColorizer(input: string): HTMLSpanElement[] {
  const colors = ['red', 'blueviolet', 'chartreuse',
    'aquamarine', 'darkcyan', 'gold', 'green', 'aqua', 'orange',
    'blue'];
  const spans: HTMLSpanElement[] = [];
  for (let i = 0; i < input.length; ++i) {
    const span = ui.span([input.at(i)]);
    // $(span).css('-webkit-text-fill-color', colors.at(Math.round(Math.random() * colors.length))!);
    $(span).css('-webkit-text-fill-color', colors.at(i % colors.length)!);
    spans.push(span);
  }
  return spans;
}
