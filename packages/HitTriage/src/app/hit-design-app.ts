import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AppName, HitDesignCampaign, HitDesignTemplate, HitTriageCampaignStatus, IFunctionArgs} from './types';
import {HitDesignInfoView} from './hit-design-views/info-view';
import {CampaignIdKey, CampaignJsonName, CampaignTableName,
  HTQueryPrefix, HTScriptPrefix, HitDesignCampaignIdKey,
  HitDesignMolColName, TileCategoriesColName, ViDColName, i18n} from './consts';
import {calculateColumns, calculateSingleCellValues, getNewVid} from './utils/calculate-single-cell';
import '../../css/hit-triage.css';
import {_package} from '../package';
import {addBreadCrumbsToRibbons, checkRibbonsHaveSubmit, modifyUrl, toFormatedDateString} from './utils';
import {HitDesignSubmitView} from './hit-design-views/submit-view';
import {getTilesViewDialog} from './hit-design-views/tiles-view';
import {HitAppBase} from './hit-app-base';
import {HitBaseView} from './base-view';
import {chemFunctionsDialog} from './dialogs/functions-dialog';
import {Subscription} from 'rxjs';
import {filter} from 'rxjs/operators';
import {defaultPermissions, PermissionsDialog} from './dialogs/permissions-dialog';

export class HitDesignApp<T extends HitDesignTemplate = HitDesignTemplate> extends HitAppBase<T> {
  multiView: DG.MultiView;
  protected _infoView: HitDesignInfoView;
  get infoView(): HitDesignInfoView {return this._infoView;}
  protected _designView?: DG.TableView;
  public _submitView?: HitDesignSubmitView;
  protected _designViewName = 'Design';
  protected _filePath = `System.AppData/HitTriage/${this.appName}/campaigns`;
  protected _campaignId?: string;
  protected _molColName: string = HitDesignMolColName;
  protected _campaign?: HitDesignCampaign;
  public campaignProps: {[key: string]: any} = {};

  protected currentDesignViewId?: string;
  public mainView: DG.ViewBase;
  protected get version() {return this._campaign?.version ?? 0;};
  constructor(c: DG.FuncCall, an: AppName = 'Hit Design',
    infoViewConstructor: (app: HitDesignApp) => HitDesignInfoView = (app) => new HitDesignInfoView(app)) {
    super(c, an);
    this._infoView = infoViewConstructor(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitBaseView)
        (this.multiView.currentView as HitBaseView<T, typeof this>).onActivated();
    });
    this.multiView.parentCall = c;

    this.mainView = this.multiView;
    this._initViewChangeSub();
  }

  _initViewChangeSub() {
    this.multiView.subs.push(grok.events.onCurrentViewChanged.subscribe(async () => {
      try {
        if (grok.shell.v?.name === this.currentDesignViewId) {
          grok.shell.windows.showHelp = false;

          this.setBaseUrl();
          modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');

          const {sub} = addBreadCrumbsToRibbons(grok.shell.v, this.appName, grok.shell.v?.name, () => {
            grok.shell.v = this.mainView;
            this._designView?.close();
            this._infoView.init();
            sub.unsubscribe();
          });
        }
      } catch (e) {
        console.error(e);
      }
    }));
  }

  public async setTemplate(template: T, campaignId?: string) {
    if (!campaignId) {
      this._designView?.dataFrame && grok.shell.closeTable(this._designView.dataFrame);
      this._designView = undefined;
      campaignId = await this.getNewCampaignName(`${this.appName}/campaigns`, template.key);
      modifyUrl(HitDesignCampaignIdKey, campaignId);
      this._filePath = `System.AppData/HitTriage/${this.appName}/campaigns/${campaignId}/${CampaignTableName}`;
    } else {
      const fileLoc = `System.AppData/HitTriage/${this.appName}/campaigns`;
      this._filePath = this.campaign?.savePath ?? `${fileLoc}/${campaignId}/${CampaignTableName}`;
      this.dataFrame = await grok.dapi.files.readCsv(this._filePath);
      if (this._campaign?.columnSemTypes) {
        Object.entries(this._campaign.columnSemTypes).forEach(([colName, semtype]) => {
          const col = this.dataFrame!.columns.byName(colName);
          if (col && semtype)
            col.semType = semtype;
        });
      }
      if (this._campaign?.columnTypes) {
        Object.entries(this._campaign.columnTypes).forEach(([colName, type]) => {
          const col = this.dataFrame!.columns.byName(colName);
          if (col && col.type !== type)
            try {this.dataFrame!.changeColumnType(colName, type as DG.COLUMN_TYPE);} catch (e) {console.error(e);}
        });
      }
      //await this.dataFrame.meta.detectSemanticTypes();
    }

    if (!this.dataFrame) {
      console.error('DataFrame is empty');
      return;
    }
    await this.dataFrame.meta.detectSemanticTypes();
    this._molColName = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? HitDesignMolColName;
    this._campaignId = campaignId;
    this.template = template;

    this._submitView ??= new HitDesignSubmitView(this);
    grok.shell.windows.showHelp = false;
    if (this.campaign)
      await this.setCanEdit(this.campaign);
    else
      this.hasEditPermission = true; // if the campaign is new, obviously the user can edit it
    const designV = this.designView;
    this.currentDesignViewId = designV.name;
    this.setBaseUrl();
    modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');
    if (this.campaign)
      this.campaign.template = template;
  }

  get campaignId(): string | undefined {return this._campaignId;}

  get designView(): DG.TableView {return this._designView = this.getDesignView();}

  get designViewName(): string {return this._designViewName;}

  get molColName() {
    return this._molColName ??= this.dataFrame?.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? HitDesignMolColName;
  }

  get campaign(): HitDesignCampaign | undefined {return this._campaign;}

  set campaign(campaign: HitDesignCampaign | undefined) {this._campaign = campaign;}

  private _duplicateVidCache?: {
    colVersion: number,
    valueCounts: Int32Array | Uint32Array,
    indexes: Int32Array | Uint32Array,
  };

  private cacheDuplicateVIDs() {
    const col = this.dataFrame?.col(ViDColName);
    if (!col) {
      this._duplicateVidCache = undefined;
      return;
    }
    const colVersion = col.version;
    const valueCounts = new Uint32Array(col.categories.length).fill(0);
    const indexes = (col.getRawData() as Uint32Array).subarray(0, col.length);
    indexes.forEach((v) => valueCounts[v]++);
    this._duplicateVidCache = {colVersion, valueCounts, indexes};
  }

  get duplicateVidCache() {
    const col = this.dataFrame?.col(ViDColName);
    if (!col)
      return undefined;
    if (!this._duplicateVidCache || col.version !== this._duplicateVidCache.colVersion)
      this.cacheDuplicateVIDs();
    return this._duplicateVidCache;
  }

  async performSingleCellCalculations(newValueIdx: number, newValue?: string) {
    const computeObj = this.template!.compute;
    if (!newValue || newValue === '')
      return;

    const calcDf =
              await calculateSingleCellValues(
                newValue, computeObj.descriptors.args, computeObj.functions, computeObj.scripts, computeObj.queries);

    for (const col of calcDf.columns.toList()) {
      if (col.name === HitDesignMolColName) continue;
      if (!this.dataFrame!.columns.contains(col.name)) {
        const newCol = this.dataFrame!.columns.addNew(col.name, col.type);
        newCol.semType = col.semType;
      }
            this.dataFrame!.col(col.name)!.set(newValueIdx, col.get(0), false);
    }
          this.dataFrame!.fireValuesChanged();
          this.saveCampaign(undefined, false);
  }

  protected initDesignViewRibbons(view: DG.TableView, subs: Subscription[]) {
    const onRemoveSub = grok.events.onViewRemoved.subscribe((v) => {
      if (v.id === view?.id) {
        subs.forEach((s) => s.unsubscribe());
        onRemoveSub.unsubscribe();
      }
    });
    const ribbons = view?.getRibbonPanels();
    if (ribbons) {
      const hasSubmit = checkRibbonsHaveSubmit(ribbons);
      if (!hasSubmit) {
        const getComputeDialog = async () => {
          chemFunctionsDialog(this, async (resultMap) => {
            const oldDescriptors = this.template!.compute.descriptors.args;
            const oldFunctions = this.template!.compute.functions;
            const oldScripts = this.template!.compute.scripts ?? [];
            const oldQueries = this.template!.compute.queries ?? [];
            const newDescriptors = resultMap.descriptors;
            const newComputeObj = {
              descriptors: {
                enabled: !!resultMap?.descriptors?.length,
                args: resultMap?.descriptors ?? [],
              },
              functions: Object.entries(resultMap?.externals ?? {}).map(([funcName, args]) => {
                const splitFunc = funcName.split(':');
                return ({
                  name: splitFunc[1],
                  package: splitFunc[0],
                  args: args,
                });
              }),
              scripts: Object.entries(resultMap?.scripts ?? {})
                .filter(([name, _]) => name.startsWith(HTScriptPrefix) && name.split(':').length === 3)
                .map(([scriptId, args]) => {
                  const scriptNameParts = scriptId.split(':');
                  return ({
                    name: scriptNameParts[1] ?? '',
                    id: scriptNameParts[2] ?? '',
                    args: args,
                  });
                }),
              queries: Object.entries(resultMap?.queries ?? {})
                .filter(([name, _]) => name.startsWith(HTQueryPrefix) && name.split(':').length === 3)
                .map(([queryName, args]) => {
                  const queryNameParts = queryName.split(':');
                  return ({
                    name: queryNameParts[1] ?? '',
                    id: queryNameParts[2] ?? '',
                    args: args,
                  });
                }),
            };
            this.template!.compute = newComputeObj;
            this.campaign!.template = this.template;
            const uncalculatedDescriptors = newDescriptors.filter((d) => !oldDescriptors.includes(d));

            const newFunctions: {[_: string]: IFunctionArgs} = {};
            Object.entries(resultMap?.externals ?? {})
              .filter(([funcName, args]) => {
                const oldFunc = oldFunctions.find((f) => `${f.package}:${f.name}` === funcName);
                if (!oldFunc)
                  return true;
                return !Object.entries(args).every(([key, value]) => oldFunc.args[key] === value);
              }).forEach(([funcName, args]) => {newFunctions[funcName] = args;});

            const newScripts: {[_: string]: IFunctionArgs} = {};
            Object.entries(resultMap?.scripts ?? {})
              .filter(([scriptName, args]) => {
                const oldScript = oldScripts.find((f) => `${HTScriptPrefix}:${f.name}:${f.id}` === scriptName);
                if (!oldScript)
                  return true;
                return !Object.entries(args).every(([key, value]) => oldScript.args[key] === value);
              }).forEach(([scriptName, args]) => {newScripts[scriptName] = args;});

            const newQueries: {[_: string]: IFunctionArgs} = {};
            Object.entries(resultMap?.queries ?? {})
              .filter(([queryName, args]) => {
                const oldQuery = oldQueries.find((f) => `${HTQueryPrefix}:${f.name}:${f.id}` === queryName);
                if (!oldQuery)
                  return true;
                return !Object.entries(args).every(([key, value]) => oldQuery.args[key] === value);
              }).forEach(([queryName, args]) => {newQueries[queryName] = args;});
            ui.setUpdateIndicator(view.grid.root, true);
            try {
              await calculateColumns({descriptors: uncalculatedDescriptors, externals: newFunctions,
                scripts: newScripts, queries: newQueries}, this.dataFrame!, this.molColName!);
              this.dataFrame!.fireValuesChanged();
            } finally {
              ui.setUpdateIndicator(view.grid.root, false);
              this.saveCampaign(undefined, false);
            }
          }, () => null, this.campaign?.template!, true);
        };

        const calculateRibbon = ui.iconFA('wrench', getComputeDialog, 'Calculate additional properties');
        const addNewRowButton = ui.icons.add(() => {this.dataFrame?.rows.addNew(null, true);}, 'Add new row');
        const permissionsButton = ui.iconFA('share', async () => {
          await (new PermissionsDialog(this.campaign?.permissions)).show((res) => {
            this.campaign!.permissions = res;
            this.saveCampaign(undefined, true);
          });
        }, 'Edit campaign permissions');
        const tilesButton = ui.bigButton('Progress tracker', () => {
          getTilesViewDialog(this, () => this._designView ?? null);
        });

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
        const ribbonButtons: HTMLElement[] = [submitButton];
        if (this.template?.stages?.length ?? 0 > 0)
          ribbonButtons.unshift(tilesButton);
        if (this.campaign && this.template && !this.campaign.template)
          this.campaign.template = this.template;

        if (this.hasEditPermission)
          ribbonButtons.unshift(permissionsButton);
        ribbonButtons.unshift(calculateRibbon);
        ribbonButtons.unshift(addNewRowButton);


        ribbons.push(ribbonButtons);
        // remove project save button from the ribbon
        view.setRibbonPanels(ribbons);
      }
    }
    return view;
  }

  protected getDesignView(): DG.TableView {
    const subs: Subscription[] = [];
    const isNew = this.dataFrame!.col(this.molColName)?.toList().every((m) => !m && m === '');
    const view = grok.shell.addTableView(this.dataFrame!);
    this._designViewName = this.campaign?.name ?? this._designViewName;
    view.name = this._designViewName;

    view._onAdded();
    const layoutViewState = this._campaign?.layout ?? this.template?.layoutViewState;
    if (layoutViewState) {
      try {
        const layout = DG.ViewLayout.fromViewState(layoutViewState);
        view.loadLayout(layout);
      } catch (e) {
        grok.shell.error('Failed to apply layout. Falling back to default layout.');
        console.error(e);
      }
    }

    if (isNew)
      grok.functions.call('Chem:editMoleculeCell', {cell: view.grid.cell(this._molColName, 0)});

    subs.push(this.dataFrame!.onRowsAdded.pipe(filter(() => !this.isJoining))
      .subscribe(() => { // TODO, insertion of rows in the middle
        try {
          if (this.template!.stages?.length > 0) {
            for (let i = 0; i < this.dataFrame!.rowCount; i++) {
              const colVal = this.dataFrame!.col(TileCategoriesColName)!.get(i);
              if (!colVal || colVal === '' || this.dataFrame!.col(TileCategoriesColName)?.isNone(i))
                this.dataFrame!.set(TileCategoriesColName, i, this.template!.stages[0]);
            }
          }
          let lastAddedCell: DG.GridCell | null = null;
          for (let i = 0; i < this.dataFrame!.rowCount; i++) {
            const cell = view.grid.cell(this.molColName, i);
            if (!cell)
              continue;
            if (cell.cell.value === '' || cell.cell.value === null)
              lastAddedCell = cell;
          }
          if (lastAddedCell)
            grok.functions.call('Chem:editMoleculeCell', {cell: lastAddedCell});
        } catch (e) {
          console.error(e);
        }
        //const lastCell = view.grid.cell(this.molColName, this.dataFrame!.rowCount - 1);
        //view.grid.onCellValueEdited
      }));
    subs.push(grok.events.onContextMenu.subscribe((args) => {
      try {
        const viewer: DG.Viewer = args?.args?.context;
        if (!viewer)
          return;
        if (viewer?.type !== DG.VIEWER.GRID)
          return;
        if (!viewer.tableView || viewer.tableView.id !== view.id)
          return;
        if (args?.args?.item?.tableColumn?.name !== this.molColName || !args?.args?.item?.isTableCell)
          return;
        const menu: DG.Menu = args?.args?.menu;
        if (!menu)
          return;
        menu.item('Add new row', () => {
            this.dataFrame!.rows.addNew(null, true);
        });
        menu.item('Duplicate molecule', () => {
          try {
            const row = this.dataFrame!.rows.addNew(null, true);
            const cell = row.get(this.molColName);
            if (cell != null && row.idx > -1)
              this.dataFrame!.cell(row.idx, this.molColName).value = args?.args?.item?.cell?.value ?? '';
          } catch (e) {
            console.error(e);
          }
        });

        const cellIndex = args?.args?.item?.tableRowIndex;
        const cellValue = args?.args?.item?.cell?.value;
        if (cellValue && (cellIndex ?? -1) > -1) {
          menu.item('Re-Run Calculations', async () => {
            try {
              await this.performSingleCellCalculations(cellIndex, cellValue);
            } catch (e) {
              console.error(e);
            }
          });
        }
      } catch (e: any) {
        grok.log.error(e);
      }
    }));

    if (!view?.grid) {
      grok.shell.error('Applied layout created view without grid. Resetting layout.');
      view.resetLayout();
    }
    view?.grid && subs.push(view.grid.onCellValueEdited.subscribe(async (gc) => {
      try {
        if (gc.tableColumn?.name === TileCategoriesColName) {
          await this.saveCampaign(undefined, false);
          return;
        }
        if (gc.tableColumn?.name !== this.molColName)
          return;
        const newValue = gc.cell.value;
        const newValueIdx = gc.tableRowIndex!;
        let newVid = this.dataFrame!.col(ViDColName)?.get(newValueIdx);
        let foundMatch = false;
        // try to find existing molecule
        if (newValue) {
          try {
            const canonicals = gc.tableColumn.toList().map((cv) => {
              try {
                return grok.chem.convert(cv, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
              } catch (e) {
                return '';
              }
            },
            );
            const canonicalNewValue =
                grok.chem.convert(newValue, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
            if (canonicals?.length === this.dataFrame!.rowCount) {
              for (let i = 0; i < canonicals.length; i++) {
                if (canonicals[i] === canonicalNewValue &&
                      i !== newValueIdx && this.dataFrame!.col(ViDColName)?.get(i)) {
                  newVid = this.dataFrame!.col(ViDColName)?.get(i);
                  foundMatch = true;
                  break;
                }
              }
            }
          } catch (e) {
            console.error(e);
          }
        }
        // if the vid was duplicated, generate a new one
        if (this.duplicateVidCache && !foundMatch &&
            this.duplicateVidCache.valueCounts[this.duplicateVidCache.indexes[newValueIdx]] > 1)
          newVid = null;

        if (!newVid || newVid === '')
          newVid = getNewVid(this.dataFrame!.col(ViDColName)!);

          this.dataFrame!.col(ViDColName)!.set(newValueIdx, newVid, false);

          this.performSingleCellCalculations(newValueIdx, newValue);
      } catch (e) {
        console.error(e);
      }
    }));

    view?.grid && subs.push(view.grid.onCellRender.subscribe((args) => {
      try {
        // color duplicate vid values
        const cell = args.cell;
        if (!cell || !cell.isTableCell || !cell.tableColumn || !this.duplicateVidCache ||
            cell.tableColumn.name !== ViDColName || (cell.tableRowIndex ?? -1) < 0)
          return;

        if (this.duplicateVidCache.valueCounts[this.duplicateVidCache.indexes[cell.tableRowIndex!]] > 1) {
          args.cell.style.backColor =
              DG.Color.setAlpha(DG.Color.getCategoricalColor(this.duplicateVidCache.indexes[cell.tableRowIndex!])
                , 150);
        }
      } catch (e) {}
    }));
    this.initDesignViewRibbons(view, subs);
    view.parentCall = this.parentCall;
    return view;
  }

  getSummary(): {[_: string]: any} {
    const getFolderPath = () => {
      const campaignIndex = this._filePath.indexOf((this.campaignId ?? this._campaign?.name)!);
      const folderPath = campaignIndex === -1 ? this._filePath : this._filePath.substring(0, campaignIndex - 1);
      return folderPath;
    };
    const campaignProps = {...(this.campaign?.campaignFields ?? this.campaignProps)};
    if (this.template && this.template.campaignFields) {
      Object.entries(campaignProps).forEach(([key, value]) => {
        const field = this.template!.campaignFields!.find((f) => f.name === key);
        if (field && field.type === 'Date' && value) {
          const dateString = (new Date(value)).toString();
          if (dateString !== 'Invalid Date')
            campaignProps[key] = (new Date(value)).toLocaleDateString();
          else {
            if (value.date)
              campaignProps[key] = new Date(value.date).toLocaleDateString();
            else
              campaignProps[key] = value;
          }
        } else {
          if (field && field.type === DG.SEMTYPE.MOLECULE && value)
            campaignProps[key] = grok.chem.drawMolecule(value);
        }
      });
    }

    const getPathEditor = () => {
      const editIcon = ui.icons.edit(() => {
        const folderPath = getFolderPath();
        const newPathInput = ui.input.string('Path', {value: folderPath});
        const labelElement = newPathInput.root.getElementsByTagName('label').item(0);
        if (labelElement)
          labelElement.remove();
        newPathInput.root.style.width = '100%';

        async function checkFolder() {
          const newPath = newPathInput.value;
          if (!newPath || newPath.trim() === '') {
            grok.shell.error('Path can not be empty');
            return false;
          }
          const exists = await grok.dapi.files.exists(newPath);
          if (!exists) {
            grok.shell.error('Given folder does not exist');
            return false;
          }
          return true;
        }

        const saveButton = ui.button('Save', async () => {
          const exists = await checkFolder();
          if (!exists)
            return;
          const newPath = newPathInput.value;
          if (newPath.endsWith('/'))
            this._filePath = `${newPath}${this.campaignId}/${CampaignTableName}`;
          else
            this._filePath = `${newPath}/${this.campaignId}/${CampaignTableName}`;

          if (this._campaign)
            this._campaign!.savePath = this._filePath;
          await this.saveCampaign(undefined, true);
          ui.empty(pathDiv);
          const folderPath = getFolderPath();
          link = ui.link(folderPath,
            () => this.download(this.dataFrame!, this.campaignId ?? this._campaign?.name ?? 'Molecules'),
            i18n.download);
          pathDiv.appendChild(link);
          pathDiv.appendChild(editIcon);
        });
        const cancelButton = ui.button('Cancel', () => {
          ui.empty(pathDiv);
          pathDiv.appendChild(link);
          pathDiv.appendChild(editIcon);
        });
        ui.empty(pathDiv);
        newPathInput.addOptions(cancelButton);
        newPathInput.addOptions(saveButton);
        pathDiv.appendChild(newPathInput.root);
      }, 'Edit file path');
      editIcon.style.marginLeft = '5px';
      const folderPath = getFolderPath();
      let link = ui.link(folderPath,
        () => this.download(this.dataFrame!, this.campaignId ?? this._campaign?.name ?? 'Molecules'), i18n.download);
      const pathDiv = ui.divH([link, editIcon], {style: {alignItems: 'center'}});
      return pathDiv;
    };

    return {
      'Template': this.template?.name ?? 'Molecules',
      'File path': getPathEditor(),
      ...campaignProps,
      'Number of molecules': (this.dataFrame!.rowCount).toString(),
      'Enrichment methods': [this.template!.compute.descriptors.enabled ? 'descriptors' : '',
        ...this.template!.compute.functions.map((func) => func.name)].filter((f) => f && f.trim() !== '').join(', '),
    };
  }

  async saveCampaign(status?: HitTriageCampaignStatus, notify = true): Promise<HitDesignCampaign> {
    const campaignId = this.campaignId!;
    const templateName = this.template!.name;
    const enrichedDf = this.dataFrame!;
    const campaignName = campaignId;
    const columnSemTypes: {[_: string]: string} = {};
    enrichedDf.columns.toList().forEach((col) => columnSemTypes[col.name] = col.semType);
    const colTypeMap: {[_: string]: string} = {};
    enrichedDf.columns.toList().forEach((col) => colTypeMap[col.name] = col.type);
    const sketchStateString = this.campaign?.tilesViewerFormSketch ?? undefined;

    // if its first time save author as current user, else keep the same
    const authorUserId = this.campaign?.authorUserId ?? grok.shell.user.id;
    const permissions = this.campaign?.permissions ?? defaultPermissions;

    const campaign: HitDesignCampaign = {
      name: campaignName,
      templateName,
      status: status ?? this.campaign?.status ?? 'In Progress',
      createDate: this.campaign?.createDate ?? toFormatedDateString(new Date()),
      campaignFields: this.campaign?.campaignFields ?? this.campaignProps,
      columnSemTypes,
      rowCount: enrichedDf.rowCount,
      filteredRowCount: enrichedDf.filter.trueCount,
      savePath: this._filePath,
      columnTypes: colTypeMap,
      tilesViewerFormSketch: sketchStateString,
      version: this.version + 1,
      authorUserId,
      permissions,
    };
    if (!this.hasEditPermission) {
      grok.shell.error('You do not have permission to modify this campaign');
      return campaign;
    }
    const campaignPath = `${this.appName}/campaigns/${campaignId}/${CampaignJsonName}`;
    // check if someone already saved the campaign
    let resDf = enrichedDf;


    if (await _package.files.exists(campaignPath)) {
      const prevCamp: HitDesignCampaign = JSON.parse(await _package.files.readAsText(campaignPath));
      if ((prevCamp.version ?? 0) > this.version) {
        campaign.version = Math.max(this.version, (prevCamp.version ?? 0)) + 1;
        if (await grok.dapi.files.exists(this._filePath)) {
          try {
            const prevCampDf = await grok.dapi.files.readCsv(this._filePath);
            const joined = this.unionDataframes(prevCampDf, enrichedDf, this.molColName);
            resDf = joined;
          } catch (e) {
            console.error(e);
          }
        }
      }
    }
    const csvDf = DG.DataFrame.fromColumns(
      resDf.columns.toList().filter((col) => !col.name.startsWith('~')),
    ).toCsv();
    await grok.dapi.files.writeAsText(this._filePath, csvDf);
    const newLayout = this._designView!.saveLayout();
    if (!newLayout)
      grok.shell.warning('Layout cound not be saved');
    else
      campaign.layout = newLayout.viewState;
    campaign.template = this.template;

    await _package.files.writeAsText(campaignPath,
      JSON.stringify(campaign));
    notify && grok.shell.info('Campaign saved successfully.');
    this.campaign = campaign;
    return campaign;
  }
}
