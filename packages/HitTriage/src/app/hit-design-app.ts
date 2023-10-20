import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitDesignCampaign, HitDesignTemplate, HitTriageCampaignStatus} from './types';
import {HitDesignInfoView} from './hit-design-views/info-view';
import {CampaignIdKey, CampaignJsonName, CampaignTableName, EmptyStageCellValue, HDcampaignName, HitDesignCampaignIdKey,
  HitDesignMolColName, TileCategoriesColName, ViDColName, i18n} from './consts';
import {calculateSingleCellValues, getNewVid} from './utils/calculate-single-cell';
import '../../css/hit-triage.css';
import {_package} from '../package';
import {addBreadCrumbsToRibbons, checkRibbonsHaveSubmit, modifyUrl, toFormatedDateString} from './utils';
import {HitDesignSubmitView} from './hit-design-views/submit-view';
import {HitDesignTilesView} from './hit-design-views/tiles-view';
import {HitAppBase} from './hit-app-base';
import {HitBaseView} from './base-view';

export class HitDesignApp extends HitAppBase<HitDesignTemplate> {
  multiView: DG.MultiView;

  private _infoView: HitDesignInfoView;
  private _designView?: DG.TableView;
  public _submitView?: HitDesignSubmitView;
  private _tilesView?: HitDesignTilesView;
  private _designViewName = 'Design';
  private _campaignId?: string;
  private _dfName?: string;
  private _molColName: string = HitDesignMolColName;
  private _campaign?: HitDesignCampaign;
  public campaignProps: {[key: string]: any} = {};
  private processedValues: string[] = [];
  private _extraStageColsCount = 0;

  private currentDesignViewId?: string;
  private currentTilesViewId?: string;
  public mainView: DG.ViewBase;
  constructor(c: DG.FuncCall) {
    super(c);
    this._infoView = new HitDesignInfoView(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitBaseView)
        (this.multiView.currentView as HitBaseView<HitDesignTemplate, HitDesignApp>).onActivated();
    });
    this.multiView.parentCall = c;
    this.mainView = grok.shell.addView(this.multiView);
    grok.events.onCurrentViewChanged.subscribe(async () => {
      try {
        if (grok.shell.v?.name === this.currentDesignViewId || grok.shell.v?.name === this.currentTilesViewId) {
          grok.shell.windows.showHelp = false;

          this.setBaseUrl();
          modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');
          if (grok.shell.v?.name === this.currentTilesViewId)
            await this._tilesView?.render();
          const {sub} = addBreadCrumbsToRibbons(grok.shell.v, 'Hit Design', grok.shell.v?.name, () => {
            grok.shell.v = this.mainView;
            this._tilesView?.close();
            this._designView?.close();
            sub.unsubscribe();
          });
        }
      } catch (e) {
        console.error(e);
      }
    });
  }

  public async setTemplate(template: HitDesignTemplate, campaignId?: string) {
    if (!campaignId) {
      this._designView?.dataFrame && grok.shell.closeTable(this._designView.dataFrame);
      this._designView = undefined;
      campaignId = await this.getNewCampaignName('Hit Design/campaigns', template.key);
      modifyUrl(HitDesignCampaignIdKey, campaignId);
    } else {
      this.dataFrame = await _package.files.readCsv(`Hit Design/campaigns/${campaignId}/${CampaignTableName}`);
      if (this._campaign?.columnSemTypes) {
        Object.entries(this._campaign.columnSemTypes).forEach(([colName, semtype]) => {
          const col = this.dataFrame!.columns.byName(colName);
          if (col)
            col.semType = semtype;
        });
      }
      await this.dataFrame.meta.detectSemanticTypes();
    }

    if (!this.dataFrame) {
      console.error('DataFrame is empty');
      return;
    }
    this._molColName = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? HitDesignMolColName;
    this._dfName = this.dataFrame.name ??= 'Molecules';
    this._campaignId = campaignId;
    this.template = template;

    this._submitView ??= new HitDesignSubmitView(this);
    grok.shell.windows.showHelp = false;
    //add empty rows to define stages, used for tile categories;
    const stagesRow = this.dataFrame.getCol(TileCategoriesColName);
    if (stagesRow) {
      const categories = stagesRow.categories;
      if (categories && categories.length) {
        template.stages.forEach((s) => {
          if (!categories.includes(s)) {
            const newRow = this.dataFrame!.rows.addNew();
            const idx = newRow.idx;
            this.dataFrame!.set(TileCategoriesColName, idx, s);
            this.dataFrame!.set(ViDColName, idx, EmptyStageCellValue);
          }
        });
      }
    }
    this.dataFrame.rows.filter((r) => r[ViDColName] !== EmptyStageCellValue);
    this._extraStageColsCount = this.dataFrame!.rowCount - this.dataFrame.filter.trueCount;
    const designV = grok.shell.addView(this.designView);
    this.currentDesignViewId = designV.name;
    this._tilesView = new HitDesignTilesView(this);
    this._tilesView.parentCall = this.parentCall;
    const tilesV = grok.shell.addView(this._tilesView);
    grok.shell.v = designV;
    this.currentTilesViewId = tilesV.name;
    this._tilesView.onActivated();
    this.setBaseUrl();
    modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');
  }

  get campaignId(): string | undefined {return this._campaignId;}

  get designView(): DG.TableView {return this._designView = this.getDesignView();}

  get molColName() {
    return this._molColName ??= this.dataFrame?.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? HitDesignMolColName;
  }

  get campaign(): HitDesignCampaign | undefined {return this._campaign;}

  set campaign(campaign: HitDesignCampaign | undefined) {this._campaign = campaign;}

  private getDesignView(): DG.TableView {
    const isNew = this.dataFrame!.col(this.molColName)?.toList().every((m) => !m && m === '');
    const view = DG.TableView.create(this.dataFrame!, false);
    view.name = this._designViewName;
    this.processedValues = this.dataFrame!.getCol(this.molColName).toList();
    setTimeout(async () => {
      view._onAdded();
      await new Promise((r) => setTimeout(r, 1000)); // needed for substruct filter
      // apply layout.
      const layout = (await grok.dapi.layouts.filter(`friendlyName = "${this._designViewName}"`).list())
        .find((l) => l && l.getUserDataValue(HDcampaignName) === this._campaignId);
      if (layout)
        view.loadLayout(layout);


      if (isNew)
        grok.functions.call('Chem:editMoleculeCell', {cell: view.grid.cell(this._molColName, 0)});

      this.dataFrame!.onRowsAdded.subscribe(() => { // TODO, insertion of rows in the middle
        try {
          const newRowsNum = this.dataFrame!.rowCount - this.processedValues.length;
          this.processedValues.push(...new Array(newRowsNum).fill(''));
        } catch (e) {
          console.error(e);
        }
        //const lastCell = view.grid.cell(this.molColName, this.dataFrame!.rowCount - 1);
        //view.grid.onCellValueEdited
      });
      this.dataFrame?.onRowsRemoved.subscribe(() => {
        try {
          this.processedValues = this.dataFrame!.getCol(this.molColName).toList();
        } catch (e) {
          console.error(e);
        }
      });

      this.dataFrame!.onValuesChanged.subscribe(async () => {
        try {
          let newValueIdx: number | null = null;
          for (let i = 0; i < this.processedValues.length; i++) {
            if (this.processedValues[i] !== this.dataFrame!.get(this.molColName, i)) {
              newValueIdx = i;
              break;
            }
          }
          if (newValueIdx == null) return;
          this.processedValues[newValueIdx] = this.dataFrame!.get(this.molColName, newValueIdx);
          const newCellValue: string = this.processedValues[newValueIdx];
          this.dataFrame!.col(TileCategoriesColName)!.set(newValueIdx, this.template!.stages[0], false);
          if (!this.dataFrame!.col(ViDColName)?.get(newValueIdx))
          this.dataFrame!.col(ViDColName)!.set(newValueIdx, getNewVid(this.dataFrame!.col(ViDColName)!), false);

          const computeObj = this.template!.compute;
          if (!newCellValue || newCellValue === '')
            return;
          const calcDf =
            await calculateSingleCellValues(newCellValue, computeObj.descriptors.args, computeObj.functions);

          for (const col of calcDf.columns.toList()) {
            if (col.name === HitDesignMolColName) continue;
            if (!this.dataFrame!.columns.contains(col.name)) {
              const newCol = this.dataFrame!.columns.addNew(col.name, col.type);
              newCol.semType = col.semType;
            }
          this.dataFrame!.col(col.name)!.set(newValueIdx, col.get(0), false);
          }
        this.dataFrame!.fireValuesChanged();
        } catch (e) {
          console.error(e);
        }
      },
      );
    }, 300);
    const ribbons = view?.getRibbonPanels();
    if (ribbons) {
      const hasSubmit = checkRibbonsHaveSubmit(ribbons);
      if (!hasSubmit) {
        const submitButton = ui.bigButton('Submit', () => {
          const dialogContent = this._submitView?.render();
          if (dialogContent) {
            const dlg = ui.dialog('Submit');
            dlg.add(dialogContent);
            dlg.addButton('Save', ()=>{this.saveCampaign(); dlg.close();});
            dlg.addButton('Submit', ()=>{this._submitView?.submit(); dlg.close();});
            dlg.show();
          }
        });
        submitButton.classList.add('hit-design-submit-button');
        ribbons.push([submitButton]);
        view.setRibbonPanels(ribbons);
      }
    }
    view.parentCall = this.parentCall;
    return view;
  }

  getSummary(): {[_: string]: any} {
    const campaignProps = this.campaign?.campaignFields ?? this.campaignProps;
    return {
      'Template': this.template?.name ?? 'Molecules',
      'File path': ui.divH([ui.link(this._dfName ?? 'Molecules',
        () => this.download(this.dataFrame!, 'Molecules'), i18n.download)],
      {style: {alignItems: 'center'}}),
      ...campaignProps,
      'Number of molecules': (this.dataFrame!.rowCount - this._extraStageColsCount).toString(),
      'Enrichment methods': [this.template!.compute.descriptors.enabled ? 'descriptors' : '',
        ...this.template!.compute.functions.map((func) => func.name)].filter((f) => f && f.trim() !== '').join(', '),
    };
  }

  async saveCampaign(status?: HitTriageCampaignStatus, notify = true): Promise<any> {
    const campaignId = this.campaignId!;
    const templateName = this.template!.name;
    const enrichedDf = this.dataFrame!;
    const campaignName = campaignId;
    const columnSemTypes: {[_: string]: string} = {};
    enrichedDf.columns.toList().forEach((col) => columnSemTypes[col.name] = col.semType);
    const campaign: HitDesignCampaign = {
      name: campaignName,
      templateName,
      status: status ?? this.campaign?.status ?? 'In Progress',
      createDate: this.campaign?.createDate ?? toFormatedDateString(new Date()),
      campaignFields: this.campaign?.campaignFields ?? this.campaignProps,
      columnSemTypes,
      rowCount: enrichedDf.col(ViDColName)?.toList().filter((s) => s !== EmptyStageCellValue).length ?? 0,
      filteredRowCount: enrichedDf.filter.trueCount,
    };
    await _package.files.writeAsText(`Hit Design/campaigns/${campaignId}/${CampaignJsonName}`,
      JSON.stringify(campaign));

    const csvDf = DG.DataFrame.fromColumns(
      enrichedDf.columns.toList().filter((col) => !col.name.startsWith('~')),
    ).toCsv();
    await _package.files.writeAsText(`Hit Design/campaigns/${campaignId}/${CampaignTableName}`, csvDf);

    const newLayout = this._designView!.saveLayout();
    if (!newLayout) {
      grok.shell.warning('Layout cound not be saved');
      return;
    }

    const oldLayouts = (await grok.dapi.layouts.filter(`friendlyName = "${this._designViewName}"`).list())
      .filter((l) => l && l.getUserDataValue(HDcampaignName) === campaignId);
    for (const l of oldLayouts)
      await grok.dapi.layouts.delete(l);
    //save new layout
    newLayout.setUserDataValue(HDcampaignName, campaignId);
    const l = await grok.dapi.layouts.save(newLayout);
    const allGroup = await grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
    await grok.dapi.permissions.grant(l, allGroup, true);
    notify && grok.shell.info('Campaign saved successfully.');
  }
}
