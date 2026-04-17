// Centralized hub for Diff Studio: templates, library, recent & custom models

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {u2} from '@datagrok-libraries/utils/src/u2';

import {_package} from './package';
import {getRecentModelsTable, getEquationsFromFile, extractIvpNameAndDescription} from './utils';
import {TITLE, MODEL_HINT, TEMPLATE_TITLES, EXAMPLE_TITLES, MODEL_ICON, MISC, PATH, LINK,
  LIBRARY_CHANGED_EVENT} from './ui-constants';
import {DiffStudio, EDITOR_STATE, STATE_BY_TITLE, getLink} from './app';

//@ts-ignore
import '../css/app-styles.css';

/** Diff Studio hub view with templates, library, and recent models */
export class DiffStudioHub {
  name = 'Diff Studio';
  root = ui.divV([]);
  view: DG.View;

  constructor() {
    this.view = DG.View.fromRoot(this.root);
    this.view.name = this.name;
    this.view.helpUrl = LINK.DIF_STUDIO_REL;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;
  }

  /** Render the full hub view */
  async render(): Promise<void> {
    this.renderHeader();
    await this.renderRest();
  }

  /** Render the hub header */
  renderHeader(): void {
    ui.empty(this.root);
    this.root.append(this.buildHeader());
  }

  /** Render buttons and all model sections */
  async renderRest(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      this.root.append(this.buildButtons());
      this.buildTemplates();
      await this.buildLibrary();
      await this.buildRecent();
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  /** Return the hub header element */
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

  /** Return the hub action buttons (refresh, create, upload) */
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

  /** Upload a model from a local `.ivp` file and run it */
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

  /** Build the "Recent" section from the recent models table */
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
          if (!(await grok.dapi.files.exists(info)))
            continue;
          const idx = info.lastIndexOf('/');
          const name = info.slice(idx + 1);
          const run = () => this.openFileAsPreview(info);
          const card = this.buildModelCard(name, info, run);
          const linkUrl = `${window.location.origin}/${PATH.FILE}/${info.replace(':', '.')}`;
          this.addRecentContextMenu(card, run, linkUrl);
          cards.push(card);
        } else {
          const title = info as TITLE;
          const description = MODEL_HINT.get(title) ?? '';
          const iconPath = this.getIconUrl(title);
          const state = STATE_BY_TITLE.get(title) ?? EDITOR_STATE.BASIC_TEMPLATE;
          const run = async () => {
            const solver = new DiffStudio();
            await solver.runSolverApp(undefined, state);
          };
          const card = this.buildModelCard(info, description, run, iconPath);
          const section = TEMPLATE_TITLES.includes(title) ? TITLE.TEMPL : TITLE.LIBRARY;
          const linkUrl = `${window.location.origin}${PATH.APPS_DS}/${section}/${state}${PATH.PARAM}`;
          this.addRecentContextMenu(card, run, linkUrl);
          cards.push(card);
        }
      }

      this.root.append(ui.h1('Recent'), this.buildCardsContainer(cards));
    } catch (e) {
      // silently skip if loading fails
    }
  } // buildRecent

  /** Build the "Templates" section */
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

  /** Build the "Library" section (built-in use cases + custom entries) */
  private async buildLibrary(): Promise<void> {
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

    const externalCards = await this.buildExternalModelCards();
    cards.push(...externalCards);

    this.root.append(ui.h1('Library'), this.buildCardsContainer(cards));
  }

  /** Build cards for custom Library models listed in `external-models.json` */
  private async buildExternalModelCards(): Promise<HTMLElement[]> {
    const cards: HTMLElement[] = [];
    try {
      const manifestPath = `${PATH.LIBRARY_FOLDER}/${PATH.EXTERNAL_MODELS_JSON}`;
      if (!(await grok.dapi.files.exists(manifestPath)))
        return cards;

      const text = await grok.dapi.files.readAsText(manifestPath);
      const parsed = JSON.parse(text);
      if (!Array.isArray(parsed?.models))
        return cards;

      for (const entry of parsed.models) {
        const modelPath: string | undefined = entry?.path;
        const iconFile: string | undefined = entry?.icon;
        const helpUrl: string | undefined = entry?.help;
        if (!modelPath)
          continue;

        const content = await getEquationsFromFile(modelPath);
        if (content === null)
          continue;

        const {name, description} = extractIvpNameAndDescription(content);
        const fallbackName = modelPath.slice(modelPath.lastIndexOf('/') + 1);
        const displayName = name || fallbackName;
        const iconUrl = `${_package.webRoot}files/icons/${iconFile || 'default.png'}`;

        const run = () => this.openFileAsPreview(modelPath);

        const card = this.buildModelCard(displayName, description, run, iconUrl);
        this.addCustomModelContextMenu(card, modelPath, helpUrl, run);
        cards.push(card);
      }
    } catch {
      // silently skip if the manifest is missing or malformed
    }
    return cards;
  } // buildExternalModelCards

  /** Open an IVP file as a Diff Studio preview */
  private async openFileAsPreview(path: string): Promise<void> {
    if (!(await grok.dapi.files.exists(path))) {
      grok.shell.warning(`File not found: ${path}`);
      return;
    }
    const folderPath = path.slice(0, path.lastIndexOf('/') + 1);
    const fileList = await grok.dapi.files.list(folderPath);
    const file = fileList.find((f) => f.nqName === path);
    if (!file) {
      grok.shell.warning(`File not found: ${path}`);
      return;
    }
    const solver = new DiffStudio(false, true, true);
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showBrowse = true;
    grok.shell.addView(await solver.getFilePreview(file, path));
  }

  /** Attach a context menu to a recent-model card */
  private addRecentContextMenu(card: HTMLElement, run: () => void, linkUrl: string): void {
    card.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      DG.Menu.popup()
        .item('Run', run, null, {description: 'Run model'})
        .item('Copy link', async () => {
          await navigator.clipboard.writeText(linkUrl);
          grok.shell.info('Model link copied to clipboard');
        }, null, {description: 'Copy the model link to clipboard'})
        .show({x: ev.clientX, y: ev.clientY, causedBy: ev});
    });
  }

  /** Open the settings dialog for a custom Library model */
  private openCustomModelSettings(modelPath: string, currentHelpUrl: string | undefined): void {
    const linkInput = ui.input.string('Help link', {
      value: currentHelpUrl ?? '',
      nullable: false,
      tooltipText: 'Set the link to the help page',
    });

    const dlg = ui.dialog({title: 'Model settings'})
      .add(linkInput)
      .addButton('Save', async () => {
        const value = linkInput.value?.trim();
        if (!value)
          return;
        try {
          const manifestPath = `${PATH.LIBRARY_FOLDER}/${PATH.EXTERNAL_MODELS_JSON}`;
          const text = await grok.dapi.files.readAsText(manifestPath);
          const parsed = JSON.parse(text);
          const entry = parsed?.models?.find?.((m: {path?: string}) => m?.path === modelPath);
          if (!entry) {
            grok.shell.error('Model entry not found in external-models.json');
            return;
          }
          entry.help = value;
          await grok.dapi.files.writeAsText(manifestPath, JSON.stringify(parsed, null, 2));
          dlg.close();
          grok.events.fireCustomEvent(LIBRARY_CHANGED_EVENT, null);
          await this.render();
        } catch (e) {
          grok.shell.error(`Failed to save settings: ${e instanceof Error ? e.message : 'platform issue'}`);
        }
      }, undefined, 'Save settings')
      .show();

    const saveBtn = dlg.getButton('Save');
    saveBtn.disabled = !linkInput.value?.trim();
    linkInput.onChanged.subscribe(() => {
      saveBtn.disabled = !linkInput.value?.trim();
    });
  } // openCustomModelSettings

  /** Attach a context menu to a custom Library model card */
  private addCustomModelContextMenu(card: HTMLElement, modelPath: string, helpUrl: string | undefined,
    run: () => void): void {
    card.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      DG.Menu.popup()
        .item('Run', run, null, {description: 'Run model'})
        .item('Copy link', async () => {
          const url = `${window.location.origin}/${PATH.FILE}/${modelPath.replace(':', '.')}`;
          await navigator.clipboard.writeText(url);
          grok.shell.info('Model link copied to clipboard');
        }, null, {description: 'Copy the model link to clipboard'})
        .item('Help', () => {
          if (helpUrl)
            window.open(helpUrl, '_blank');
          else
            grok.shell.warning('Help link is not defined. To set help, right click the model and select "Settings..."');
        }, null, {description: 'Open help in a new tab'})
        .item('Settings...', () => this.openCustomModelSettings(modelPath, helpUrl), null,
          {description: 'Configure model properties'})
        .show({x: ev.clientX, y: ev.clientY, causedBy: ev});
    });
  }

  /** Attach a context menu to a built-in template or library model card */
  private addModelContextMenu(card: HTMLElement, section: TITLE, state: EDITOR_STATE, run: () => void): void {
    card.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      DG.Menu.popup()
        .item('Run', run, null, {description: 'Run model'})
        .item('Copy link', async () => {
          const url = `${window.location.origin}${PATH.APPS_DS}/${section}/${state}${PATH.PARAM}`;
          await navigator.clipboard.writeText(url);
          grok.shell.info('Model link copied to clipboard');
        }, null, {description: 'Copy the model link to clipboard'})
        .item('Help', () => {
          grok.shell.windows.help.visible = true;
          grok.shell.windows.help.showHelp(getLink(state));
        }, null, {description: 'Open help for this model'})
        .show({x: ev.clientX, y: ev.clientY, causedBy: ev});
    });
  }

  /** Return the icon URL for a built-in model */
  private getIconUrl(title: TITLE): string {
    const iconFile = MODEL_ICON.get(title);
    return iconFile ? `${_package.webRoot}files/${iconFile}` : `${_package.webRoot}/package.png`;
  }

  /** Build a model card (icon + name + tooltip, double-click to run) */
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
  } // buildModelCard

  /** Wrap a list of cards into a grid container */
  private buildCardsContainer(cards: HTMLElement[]): HTMLElement {
    const container = ui.div(cards);
    container.classList.add('diff-studio-hub-grid');
    return container;
  }

  /** Render the hub and show it as a preview */
  async show(): Promise<void> {
    await this.render();
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showBrowse = true;
    grok.shell.addPreview(this.view);
  }
}
