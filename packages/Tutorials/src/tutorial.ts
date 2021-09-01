import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { Observable } from 'rxjs';
import { filter, first } from 'rxjs/operators';
import { _package } from './package';


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
  header: HTMLHeadingElement = ui.h2('');
  subheader: HTMLHeadingElement = ui.h3('');
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
    }

    await this._run();
    
    this.title('Congratulations!');
    this.describe('You have successfully completed this tutorial.');
  }

  title(text: string): void {
    this.activity.append(ui.h2(text));
  }

  describe(text: string): void {
    const div = ui.div([], 'grok-tutorial-description-entry');
    div.innerHTML = text;
    this.activity.append(div);
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
    hintSub?.cancel();
    hint?.classList?.remove('grok-tutorial-target-hint');
  }

  protected _scroll(): void {
    this.root.scrollTop = this.root.scrollHeight;
  }

  firstEvent(eventStream: Observable<any>): Promise<void> {
    return new Promise<void>((resolve, reject) => {
      eventStream.pipe(first()).subscribe((_: any) => resolve());
    });
  }

  async openPlot(name: string, check: (viewer: DG.Viewer) => boolean): Promise<DG.Viewer> {
    // TODO: Expand toolbox / accordion API coverage
    const getViewerIcon = (el: HTMLElement) => $(el).find(`i.svg-${name.replace(' ', '-')}`).get()[0];
    const view: DG.View = grok.shell.v;
    let viewer: DG.Viewer;

    await this.action(`Open ${name}`,
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
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

  async dlgInputAction(dlg: DG.Dialog, instructions: string, caption: string, value: string) {
    const inp = dlg.inputs.filter((input: DG.InputBase) => input.caption == caption)[0];
    await this.action(instructions,
      new Observable((subscriber: any) => inp.onChanged(() => {
        if (inp.stringValue === value) subscriber.next(inp.stringValue);
      })),
      inp.root,
    );
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

export class TutorialRunner {
  root: HTMLDivElement = ui.div([], 'grok-tutorial,grok-welcome-panel');

  async run(t: Tutorial): Promise<void> {
    $(this.root).empty();
    this.root.append(t.root);
    await t.run();
  }

  constructor(track: Track, onStartTutorial?: (t: Tutorial) => Promise<void>) {
    this.root.append(ui.h2(`${track.name} Tutorials`));
    this.root.append(ui.divV(track.tutorials.map((t) => {
      const el = new TutorialCard(t).root;
      el.addEventListener('click', () => {
        if (onStartTutorial == null) {
          this.run(t);
        } else {
          onStartTutorial(t);
        }
      });
      return el;
    })));

  }

}

class TutorialCard {
  root: HTMLDivElement = ui.div();
  tutorial: Tutorial;

  constructor(tutorial: Tutorial) {
    this.tutorial = tutorial;

    this.root = ui.divH([
      ui.image(
        `${_package.webRoot}images/${tutorial.name.toLowerCase().replace(/ /g, '-')}.png`,
        100, 100),
      ui.divV([
        ui.h2(tutorial.name),
        ui.divText(tutorial.description, { id: 'description' }),
      ])
    ], 'grok-tutorial-card');
  }
}
