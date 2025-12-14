/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from './app/hit-triage-app';
import {HitDesignApp} from './app/hit-design-app';
import {GasteigerPngRenderer} from './pngRenderers';
import {loadCampaigns, timeoutOneTimeEventListener} from './app/utils';
import {AppName, CampaignsType} from './app';
import {PeptiHitApp} from './app/pepti-hit-app';
import {CampaignJsonName, PeptiHitHelmColName} from './app/consts';
import {htPackageSettingsEditorWidget} from './packageSettingsEditor';

export * from './package.g';

export class HTPackage extends DG.Package {
  molToSmilesLruCache = new DG.LruCache<string, string>(2000);
  campaignsCache: {[key in AppName]?: {[name: string]: CampaignsType[key]}} = {};

  convertToSmiles(mol: string): string {
    if (!mol)
      return '';
    return this.molToSmilesLruCache.getOrCreate(mol,
      (mol) => grok.chem.convert(mol, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles),
    );
  }

  async loadCampaigns<T extends AppName>(appName: T, deletedCampaigns?: string[]) {
    if (!this.campaignsCache[appName])
      this.campaignsCache[appName] = (await loadCampaigns(appName, deletedCampaigns ?? [])) as typeof this.campaignsCache[T];

    return this.campaignsCache[appName]! as { [name: string]: CampaignsType[T] };
  }

  async saveCampaignJson<T extends AppName>(appName: T, campaign: CampaignsType[T]) {
    const path = `${appName}/campaigns/${campaign.name}/${CampaignJsonName}`;
    await this.files.writeAsText(path, JSON.stringify(campaign));
    (this.campaignsCache[appName] as { [name: string]: CampaignsType[T] })[campaign.name] = campaign as any;
  }
}

export const _package = new HTPackage();

async function hitAppTB(treeNode: DG.TreeViewGroup, name: AppName) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  const camps = (await _package.loadCampaigns(name, []));
  let prevTable: DG.TableView | null = null;
  let firstTemplate = true;

  for (const [_, camp] of Object.entries(camps)) {
    const templateName = camp.templateName ?? camp.template?.name;
    const templateGroup = templateName ? treeNode.getOrCreateGroup(templateName, null, firstTemplate) : treeNode;
    firstTemplate = false;
    const node = templateGroup.item(camp.friendlyName ?? camp.name);

    node.root.addEventListener('dblclick', (e) => {
      e.stopImmediatePropagation();
      e.stopPropagation();
      e.preventDefault();
    });

    DG.debounce(node.onSelected, 200).subscribe(async (_) => {
      try {
        const savePath = 'ingest' in camp ? camp.ingest.query : camp.savePath;
        if (!savePath || !(await grok.dapi.files.exists(savePath))) {
          grok.shell.error('File not found: ' + savePath);
          return;
        }
        const df = await grok.dapi.files.readCsv(savePath);
        if (!df)
          return;
        df.name = camp.friendlyName ?? camp.name;
        const semtypeInfo = camp.columnSemTypes;
        if (semtypeInfo) {
          for (const [colName, semType] of Object.entries(semtypeInfo)) {
            const col = df.columns.byName(colName);
            if (col) {
              col.semType = semType;
              if (semType === DG.SEMTYPE.MACROMOLECULE && colName === PeptiHitHelmColName) {
                col.setTag('units', 'helm');
                col.setTag('.alphabetIsMultichar', 'true');
                col.setTag('cell.renderer', 'helm');
              }
            }
          }
        }
        const activeElement = document.activeElement;
        let clicked = false;
        function resetActiveElement() {
          setTimeout(() => {
            if (!clicked && activeElement && document.activeElement !== activeElement)
              activeElement && 'focus' in activeElement && (activeElement as HTMLElement).focus();
          });
        }
        window.addEventListener('focusout', resetActiveElement);
        try {
          if (prevTable)
            prevTable.close();
        } catch (e) {}
        try {
          prevTable = grok.shell.addTableView(df);
          const layout = camp.layout;
          if (layout)
            prevTable.loadLayout(DG.ViewLayout.fromViewState(layout));
        } catch (e) {
          console.error(e);
        }
        function removeListeners() {
          window.removeEventListener('focusout', resetActiveElement);
        }
        function clickListener() {
          removeListeners();
          clicked = true;
        }
        window.addEventListener('click', clickListener, {once: true});

        setTimeout(() => {
          removeListeners();
          window.removeEventListener('click', clickListener);
        }, 5000);
      } catch (e) {
        console.error(e);
      }
    });
  }
  loaderItem.remove();
}

export class PackageFunctions {
  @grok.decorators.func()
  static async hitTriageAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    await hitAppTB(treeNode, 'Hit Triage');
  }


  @grok.decorators.func()
  static async hitDesignAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    await hitAppTB(treeNode, 'Hit Design');
  }


  @grok.decorators.func()
  static async peptiHitAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    await hitAppTB(treeNode, 'PeptiHit');
  }


  @grok.decorators.app({
    browsePath: 'Chem',
    name: 'Hit Triage',
  })
  static async hitTriageApp(): Promise<DG.ViewBase> {
    const c = grok.functions.getCurrentCall();
    return new HitTriageApp(c).multiView;
  }


  @grok.decorators.app({
    icon: 'images/icons/hit-design-icon.png',
    browsePath: 'Chem',
    name: 'Hit Design',
  })
  static async hitDesignApp(): Promise<DG.ViewBase> {
    const c = grok.functions.getCurrentCall();
    return new HitDesignApp(c).multiView;
  }


  @grok.decorators.app({
    icon: 'images/icons/pepti-hit-icon.png',
    browsePath: 'Peptides',
    name: 'PeptiHit',
  })
  static async peptiHitApp(): Promise<DG.ViewBase> {
    const c = grok.functions.getCurrentCall();
    await grok.functions.call('Bio:initBio', {});
    return new PeptiHitApp(c).multiView;
  }


  @grok.decorators.func({
    tags: ['HitTriageDataSource'],
    name: 'Demo Molecules 100',
  })
  static async demoFileIngest(): Promise<DG.DataFrame> {
    const df = grok.data.demo.molecules(100);
    df.name = '100 Molecules';
    return df;
  }


  @grok.decorators.func({
    tags: ['HitTriageDataSource'],
    name: 'Demo Molecules 5000',
  })
  static async demoFileIngest1(): Promise<DG.DataFrame> {
    const df = grok.data.demo.molecules(5000);
    df.name = '5000 Molecules';
    return df;
  }


  @grok.decorators.func({
    tags: ['HitTriageDataSource'],
    name: 'Demo Molecules variable',
  })
  static async demoFileIngest2(
    @grok.decorators.param({type: 'int', options: {description: 'Molecules counts'}}) numberOfMolecules: number,
  ): Promise<DG.DataFrame> {
    const df = grok.data.demo.molecules(numberOfMolecules);
    df.name = 'Variable Molecules number';
    return df;
  }


  @grok.decorators.func({
    tags: ['HitTriageSubmitFunction'],
    name: 'Demo File Submit',
  })
  static async demoFileSubmit(
    @grok.decorators.param({options: {description: 'Dataframe'}}) df: DG.DataFrame,
    @grok.decorators.param({options: {description: 'Molecules column name'}}) molecules: string): Promise<void> {
    grok.shell.info(df.rowCount);
    grok.shell.info(molecules);
  }


  @grok.decorators.func({
    meta: {
      cellType: 'customGasteigerPNG',
      columnTags: 'quality=customGasteigerPNG',
    },
    tags: ['cellRenderer'],
    name: 'gasteigerRenderer',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static gasteigerCellRenderer(): GasteigerPngRenderer {
    return new GasteigerPngRenderer();
  }


  @grok.decorators.func({
    tags: ['packageSettingsEditor'],
    name: 'Hit Triage package settings editor',
  })
  static async htPackageSettingEditor(
    @grok.decorators.param({'name': 'propList', 'type': 'object'}) properties: DG.Property[]) : Promise<DG.Widget> {
    return htPackageSettingsEditorWidget(properties);
  }
}
