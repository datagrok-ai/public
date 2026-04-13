import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {u2} from '@datagrok-libraries/utils/src/u2';

import {_package} from './package';
import {getMyModelFiles, getRecentModelsTable, getEquationsFromFile} from './utils';
import {TITLE, MODEL_HINT, TEMPLATE_TITLES, EXAMPLE_TITLES} from './ui-constants';
import {DiffStudio, EDITOR_STATE, STATE_BY_TITLE} from './app';

import '../css/app-styles.css';

export class DiffStudioHub {
  name = 'Diff Studio Hub';
  root = ui.divV([]);
  view: DG.View;

  constructor() {
    this.view = DG.View.fromRoot(this.root);
    this.view.name = this.name;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;
  }

  async render(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      ui.empty(this.root);
      this.root.append(
        this.buildHeader(),
        this.buildButtons(),
      );

      await this.buildMyModels();
      await this.buildRecent();
      this.buildTemplates();
      this.buildLibrary();
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private buildHeader(): HTMLElement {
    const header = u2.appHeader({
      iconPath: `${_package.webRoot}/package.png`,
      learnMoreUrl: 'https://datagrok.ai/help/compute/diff-studio',
      description:
        '- Solve initial value problems for ordinary differential equations\n' +
        '- Enter formulas in a declarative form or modify templates\n' +
        '- Run computations and explore your model interactively\n' +
        '- Use sensitivity analysis and parameter fitting\n',
    });
    header.style.marginBottom = '20px';
    return header;
  }

  private buildButtons(): HTMLElement {
    const refreshBtn = ui.iconFA('sync', () => this.render(), 'Refresh view');
    const createBtn = ui.button('Create model', () => {}, 'Create new model');
    const uploadBtn = ui.button('Upload', () => {}, 'Upload model from local files');
    return ui.divH([refreshBtn, createBtn, uploadBtn], {style: {alignItems: 'center', gap: '8px'}});
  }

  private async buildMyModels(): Promise<void> {
    try {
      const myModelFiles = await getMyModelFiles();
      if (myModelFiles.length === 0)
        return;

      const cards = myModelFiles.map((file) => {
        const name = file.name;
        const description = file.fullPath;
        return this.buildModelCard(name, description, async () => {
          const equations = await getEquationsFromFile(file.fullPath);
          if (equations) {
            const solver = new DiffStudio();
            await solver.runSolverApp(equations, EDITOR_STATE.FROM_FILE);
          }
        });
      });

      this.root.append(ui.h1('My Models'), this.buildCardsContainer(cards));
    } catch (e) {
      // silently skip if loading fails
    }
  }

  private async buildRecent(): Promise<void> {
    try {
      const recentDf = await getRecentModelsTable();
      if (recentDf.rowCount === 0)
        return;

      const infoCol = recentDf.col(TITLE.INFO);
      const isCustomCol = recentDf.col(TITLE.IS_CUST);
      if (!infoCol || !isCustomCol)
        return;

      const cards: HTMLElement[] = [];
      for (let i = 0; i < recentDf.rowCount; ++i) {
        const info: string = infoCol.get(i);
        const isCustom: boolean = isCustomCol.get(i);

        if (isCustom) {
          const idx = info.lastIndexOf('/');
          const name = info.slice(idx + 1);
          cards.push(this.buildModelCard(name, info, async () => {
            const equations = await getEquationsFromFile(info);
            if (equations) {
              const solver = new DiffStudio();
              await solver.runSolverApp(equations, EDITOR_STATE.FROM_FILE);
            }
          }));
        } else {
          const title = info as TITLE;
          const description = MODEL_HINT.get(title) ?? '';
          cards.push(this.buildModelCard(info, description, async () => {
            const solver = new DiffStudio();
            await solver.runSolverApp(undefined, STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE);
          }));
        }
      }

      this.root.append(ui.h1('Recent'), this.buildCardsContainer(cards));
    } catch (e) {
      // silently skip if loading fails
    }
  }

  private buildTemplates(): void {
    const cards = TEMPLATE_TITLES.map((title) => {
      const description = MODEL_HINT.get(title) ?? '';
      return this.buildModelCard(title, description, async () => {
        const solver = new DiffStudio();
        await solver.runSolverApp(undefined, STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE);
      });
    });

    this.root.append(ui.h1('Templates'), this.buildCardsContainer(cards));
  }

  private buildLibrary(): void {
    const cards = EXAMPLE_TITLES.map((title) => {
      const description = MODEL_HINT.get(title) ?? '';
      return this.buildModelCard(title, description, async () => {
        const solver = new DiffStudio();
        await solver.runSolverApp(undefined, STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE);
      });
    });

    this.root.append(ui.h1('Library'), this.buildCardsContainer(cards));
  }

  private buildModelCard(name: string, description: string, onClick: () => void): HTMLElement {
    const icon = ui.iconImage('diff-studio', `${_package.webRoot}/package.png`);
    icon.style.width = '32px';
    icon.style.height = '32px';
    const label = ui.label(name);
    label.classList.add('d4-link-label');
    const desc = ui.divText(description, 'diff-studio-hub-card-description');
    const card = ui.divV([ui.divH([icon, label], {style: {alignItems: 'center', gap: '6px'}}), desc]);
    card.classList.add('diff-studio-hub-card');
    card.addEventListener('dblclick', onClick);
    return card;
  }

  private buildCardsContainer(cards: HTMLElement[]): HTMLElement {
    const container = ui.div(cards);
    container.classList.add('diff-studio-hub-grid');
    return container;
  }

  async show(): Promise<void> {
    await this.render();
    grok.shell.addPreview(this.view);
  }
}
