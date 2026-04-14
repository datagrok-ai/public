import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {u2} from '@datagrok-libraries/utils/src/u2';

import {_package} from './package';
import {getRecentModelsTable, getEquationsFromFile} from './utils';
import {TITLE, MODEL_HINT, TEMPLATE_TITLES, EXAMPLE_TITLES, MODEL_ICON, MISC, PATH, LINK} from './ui-constants';
import {DiffStudio, EDITOR_STATE, STATE_BY_TITLE} from './app';

import '../css/app-styles.css';

export class DiffStudioHub {
  name = 'Diff Studio Hub';
  root = ui.divV([]);
  view: DG.View;

  constructor() {
    this.view = DG.View.fromRoot(this.root);
    this.view.name = this.name;
    this.view.helpUrl = LINK.DIF_STUDIO_REL;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;
  }

  async render(): Promise<void> {
    this.renderHeader();
    await this.renderRest();
  }

  renderHeader(): void {
    ui.empty(this.root);
    this.root.append(this.buildHeader());
  }

  async renderRest(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      this.root.append(this.buildButtons());
      this.buildTemplates();
      this.buildLibrary();
      await this.buildRecent();
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
    header.classList.add('diff-studio-hub-header');
    return header;
  }

  private buildButtons(): HTMLElement {
    const refreshBtn = ui.iconFA('sync', () => this.render(), 'Refresh view');
    const createBtn = ui.button('Create', async () => {
      const solver = new DiffStudio();
      await solver.runSolverApp(undefined, EDITOR_STATE.BASIC_TEMPLATE);
      solver.showEditor();
    }, 'Create new model');
    const uploadBtn = ui.button('Upload', () => this.uploadFromLocalFile(), 'Upload model from local files');
    const buttons = ui.divH([refreshBtn, createBtn, uploadBtn]);
    buttons.classList.add('diff-studio-hub-buttons');
    return buttons;
  }

  private uploadFromLocalFile(): void {
    const fileInp = document.createElement('input');
    fileInp.type = 'file';
    fileInp.accept = `.${MISC.MODEL_FILE_EXT}`;
    fileInp.onchange = () => {
      const file = fileInp.files?.[0];
      if (!file)
        return;
      const reader = new FileReader();
      reader.addEventListener('load', async () => {
        const equations = reader.result as string;
        const solver = new DiffStudio();
        await solver.runSolverApp(equations, EDITOR_STATE.FROM_FILE);
        solver.showEditor();
      }, false);
      reader.readAsText(file);
    };
    fileInp.click();
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
          const iconPath = this.getIconUrl(title);
          cards.push(this.buildModelCard(info, description, async () => {
            const solver = new DiffStudio();
            await solver.runSolverApp(undefined, STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE);
          }, iconPath));
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
      const iconPath = this.getIconUrl(title);
      const state = STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE;
      const run = async () => {
        const solver = new DiffStudio();
        await solver.runSolverApp(undefined, state);
      };
      const card = this.buildModelCard(title, description, run, iconPath);
      this.addModelContextMenu(card, TITLE.TEMPL, state, run);
      return card;
    });

    this.root.append(ui.h1('Templates'), this.buildCardsContainer(cards));
  }

  private buildLibrary(): void {
    const cards = EXAMPLE_TITLES.map((title) => {
      const description = MODEL_HINT.get(title) ?? '';
      const iconPath = this.getIconUrl(title);
      const state = STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE;
      const run = async () => {
        const solver = new DiffStudio();
        await solver.runSolverApp(undefined, state);
      };
      const card = this.buildModelCard(title, description, run, iconPath);
      this.addModelContextMenu(card, TITLE.LIBRARY, state, run);
      return card;
    });

    this.root.append(ui.h1('Library'), this.buildCardsContainer(cards));
  }

  private addModelContextMenu(card: HTMLElement, section: TITLE, state: EDITOR_STATE, run: () => void): void {
    card.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      DG.Menu.popup()
        .item('Run', run, null, {description: 'Run model'})
        .item('Copy link', async () => {
          const url = `${window.location.origin}${PATH.APPS_DS}/${section}/${state}`;
          await navigator.clipboard.writeText(url);
          grok.shell.info('Model link copied to clipboard');
        }, null, {description: 'Copy the model link to clipboard'})
        .item('Help', () => window.open(LINK.DIF_STUDIO, '_blank'), null,
          {description: 'Open help in a new tab'})
        .show({x: ev.clientX, y: ev.clientY, causedBy: ev});
    });
  }

  private getIconUrl(title: TITLE): string {
    const iconFile = MODEL_ICON.get(title);
    return iconFile ? `${_package.webRoot}files/${iconFile}` : `${_package.webRoot}/package.png`;
  }

  private buildModelCard(name: string, description: string, onClick: () => void, iconUrl?: string): HTMLElement {
    const imgUrl = iconUrl ?? `${_package.webRoot}/files/icons/default.png`;
    const icon = ui.iconImage('diff-studio', imgUrl);
    icon.classList.add('diff-studio-hub-card-icon');
    const label = ui.label(name);
    const header = ui.divH([icon, label]);
    header.classList.add('diff-studio-hub-card-header');
    const card = ui.divV([header]);
    card.classList.add('diff-studio-hub-card');
    card.addEventListener('dblclick', onClick);

    ui.tooltip.bind(card, () => {
      const tooltipIcon = ui.iconImage('diff-studio', imgUrl);
      tooltipIcon.classList.add('diff-studio-hub-tooltip-icon');
      const tooltipName = ui.label(name);
      tooltipName.classList.add('diff-studio-hub-tooltip-name');
      const tooltipHeader = ui.divH([tooltipIcon, tooltipName]);
      tooltipHeader.classList.add('diff-studio-hub-tooltip-header');
      const tooltipDesc = ui.divText(description);
      return ui.divV([tooltipHeader, tooltipDesc]);
    });

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
