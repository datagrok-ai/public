import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Observable } from 'rxjs';
import { first } from 'rxjs/operators';


/** A base class for tutorials */
export abstract class Tutorial extends DG.Widget {
  abstract get name(): string;
  abstract get description(): string;

  demoTable: string = 'demog';

  nextLink: HTMLAnchorElement = ui.link('next',
    '',
    'Go to the next tutorial', {
    classes: 'grok-tutorial-next',
    style: { display: 'none' },
  });
  header: HTMLHeadingElement = ui.h1('');
  subheader: HTMLHeadingElement = ui.h2('');
  activity: HTMLDivElement = ui.div([], 'grok-tutorial-description');

  static DATA_STORAGE_KEY: string = 'tutorials';

  constructor() {
    super(ui.div());
    this.root.append(this.header);
    this.root.append(this.subheader);
    this.root.append(this.activity);
    this.root.append(this.nextLink);
  }

  protected abstract _run(): Promise<void>;

  async run(): Promise<void> {
    if (this.demoTable) {
      const v = grok.shell.addTableView(await grok.data.getDemoTable(this.demoTable));
      v.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, 'tutorial-widget', 0.3);
    }

    await this._run();
    grok.shell.info('Completed!');
  }

  title(text: string): void {
    this.activity.append(ui.h1(text));
  }

  describe(text: string): void {
    this.activity.append(ui.divText(text, 'grok-tutorial-description-entry'));
    this._scroll();
  }

  async action(instructions: string, completed: Observable<any>, hint?: HTMLElement, hintSub?: DG.StreamSubscription): Promise<void> {
    hint?.classList.add('grok-tutorial-target-hint');
    const instructionDiv = ui.divText(instructions, 'grok-tutorial-entry-instruction');
    const entry = ui.div([instructionDiv], 'grok-tutorial-entry');
    this.activity.append(entry);
    this._scroll();
    await this.firstEvent(completed);
    instructionDiv.classList.add('grok-tutorial-entry-success');
    this.root.append(ui.divText('Done!', 'grok-tutorial-entry-success'));
    hintSub?.cancel();
    hint?.classList?.remove('grok-tutorial-target-hint');
  }

  protected _scroll(): void {
    this.root.scrollTop = this.root.scrollHeight;
  }

  firstEvent(eventStream: Observable<any>): Promise<void> {
    return new Promise<void>((resolve, reject) => {
      eventStream.pipe(first()).subscribe((_) => resolve());
    });
  }
}


/** A collection of tutorials */
export class Track {
  tutorials: Tutorial[];
  name: string;

  constructor(name: string, ...tutorials: Tutorial[]) {
    this.name = name;
    this.tutorials = tutorials;
  }
}
