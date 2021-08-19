import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { Observable } from 'rxjs';
import { filter, first } from 'rxjs/operators';


/** A base class for tutorials */
export abstract class Tutorial extends DG.Widget {
  abstract get name(): string;
  abstract get description(): string;

  demoTable: string = 'demog.csv';
  get t(): DG.DataFrame {
    return grok.shell.t;
  }

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

  async action(instructions: string, completed: Observable<any>,
    hint?: HTMLElement | null, hintSub?: DG.StreamSubscription | null): Promise<void> {
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

  async openPlot(name: string, check: (viewer: DG.Viewer) => boolean): Promise<DG.Viewer> {
    // TODO: Expand toolbox / accordion API coverage
    const getViewerIcon = (el: HTMLElement) => $(el).find(`i.svg-${name.replace(' ', '-')}`).get()[0];
    const view: DG.View = grok.shell.v;
    let viewer: DG.Viewer;

    await this.action(`Open ${name}`,
      grok.events.onViewerAdded.pipe(filter((data) => {
        const found = check(data.args.viewer);
        if (found) {
          viewer = data.args.viewer;
        }
        return found;
      })),
      view.type === 'TableView' ? getViewerIcon((<DG.TableView>view).toolboxPage.accordion.root) : null,
    );

    return viewer!;
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
