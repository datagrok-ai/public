/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {filter, map, take} from 'rxjs/operators';
import {getContextHelp} from '../../shared-utils/utils';

function addPopover(popover: HTMLElement) {
  stylePopover(popover);
  document.body.append(popover);
}

function displayPopover(icon: HTMLElement, popover: HTMLElement) {
  alignPopover(icon, popover);
  popover.showPopover();
}

function hidePopover(popover: HTMLElement) {
  popover.hidePopover();
}

function alignPopover(target: HTMLElement, popover: HTMLElement): void {
  const bounds = target.getBoundingClientRect().toJSON();
  popover.style.inset = 'unset';
  popover.style.top = bounds.y + 'px';
  popover.style.left = (bounds.x + 20) + 'px';
}

function stylePopover(popover: HTMLElement): void {
  popover.popover = 'auto';
  popover.style.cursor = 'default';
  popover.style.padding = '10px';
  popover.style.background = '#fdffe5';
  popover.style.border = '1px solid #E4E6CE';
  popover.style.borderRadius = '2px';
  popover.style.boxShadow = '0 0 5px #E4E6CE';
  popover.style.maxWidth = '500px';
  popover.style.color = 'var(--gray-6)';
}

function getBulletIcon() {
  const bulletIcon = ui.iconFA('times');
  bulletIcon.style.marginRight = '6px';
  bulletIcon.classList.remove('grok-icon');
  return bulletIcon;
}

async function requestMembership(groupName: string) {
  try {
    const groups = await grok.dapi.groups.filter(`friendlyName="${groupName}"`).list();
    if (groups.length === 0)
      throw new Error(`group with specified name is not found`);

    if (groups.length > 1)
      throw new Error(`found more than one group with specified name`);

    const group = groups[0];

    // Workaround till JS API is not ready: https://reddata.atlassian.net/browse/GROK-14160
    await fetch(`${window.location.origin}/api/groups/${group.id}/requests/${grok.shell.user.group.id}`, {method: 'POST'});

    grok.shell.info(`Request to join ${groupName} has been initiated. Please allow some time for approval.`);
  } catch (e: any) {
    grok.shell.error(e.toString());
  }
}

export class ModelHandler extends DG.ObjectHandler {
  override get type() {
    return 'Model';
  }

  override get helpUrl() {
    return null;
  }

  override async getById(id: string): Promise<DG.Func> {
    return await grok.dapi.functions.find(id);
  }

  override async refresh(x: DG.Script): Promise<DG.Func> {
    return await this.getById(x.id);
  }

  static async openHelp(func: DG.Func) {
    const readmeText = await getContextHelp(func);
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(ui.markdown(readmeText ?? ''));
  }

  static openModel(x: DG.Func) {
    const fc = x.prepare();
    fc.edit();
  }

  static openModelFromFunccall(funcCall: DG.FuncCall) {
    funcCall.edit();
  }

  static async getUserGroups() {
    // Workaround till JS API is not ready: https://reddata.atlassian.net/browse/GROK-14159
    const userGroups = (await(await fetch(`${window.location.origin}/api/groups/all_parents`)).json() as DG.Group[]);
    return userGroups;
  }

  static getMissingGroups(func: DG.Func, userGroups: DG.Group[]) {
    const mandatoryUserGroups = JSON.parse(
      func.options['mandatoryUserGroups'] ? `${func.options['mandatoryUserGroups']}` : '[]') as {name: string, help?: string}[];
    const missingMandatoryGroups = mandatoryUserGroups
      .filter((group) => !userGroups.map((userGroup) => userGroup.friendlyName).includes(group.name));
    return missingMandatoryGroups;
  }

  static openVersionSelectDialog(func: DG.Func) {
    const versions = JSON.parse(func.options['versions'] ? `${func.options['versions']}` : '[]') as string[];
    if ((versions?.length ?? 0) < 2) {
      grok.shell.info('No other versions avaliable');
      return;
    }
    const dlg = ui.dialog({title: 'Choose version'});
    const selInput = ui.input.choice('Select version', {items: versions, value: versions[0]});
    dlg.add(selInput).onOK(() => {
      const version = selInput.value;
      const fc = func.prepare({params: {version}});
      fc.edit();
    });
    dlg.show();
  }

  // Checks whether this is the handler for [x]
  override isApplicable(x: any) {
    const js = DG.toJs(x);
    return js instanceof DG.Func && js.hasTag('model');
  }

  private userGroups = new BehaviorSubject<DG.Group[] | undefined>(undefined);

  constructor() {
    super();
  }

  getLanguageIcon(language: string) {
    if (language == 'grok')
      return ui.iconSvg('project');
    return ui.iconImage('script', `/images/entities/${language}.png`);
  }

  override renderMarkup(x: DG.Func): HTMLElement {
    const markup = ui.divH([], {style: {justifyContent: 'space-between', width: '100%'}});

    setTimeout(async () => {
      const userGroups = await this.awaitUserGroups();
      const missingMandatoryGroups = ModelHandler.getMissingGroups(x, userGroups);
      const hasMissingMandatoryGroups = missingMandatoryGroups.length > 0;
      const mandatoryGroupsIcon = ui.iconFA('exclamation-triangle', null);
      mandatoryGroupsIcon.classList.remove('grok-icon');

      const mandatoryGroupsInfo = this.makeMandatoryGroupsInfo(missingMandatoryGroups);

      const label = ui.span([
        this.renderIcon(x),
        ui.label(x.friendlyName),
      ]);

      markup.onclick = () => {
        if (grok.shell.windows.help.visible)
          ModelHandler.openHelp(x);
      };

      if (!hasMissingMandatoryGroups)
        markup.ondblclick = () => {ModelHandler.openModel(x);};
      else
        label.classList.add('d4-disabled');

      if (hasMissingMandatoryGroups) {
        addPopover(mandatoryGroupsInfo);

        markup.addEventListener('click', () => {
          displayPopover(mandatoryGroupsIcon, mandatoryGroupsInfo);
        });
      }

      markup.append(label);
      if (hasMissingMandatoryGroups)
        markup.append(ui.span([mandatoryGroupsIcon], {style: {paddingLeft: '10px'}}));
    });

    return markup;
  }

  override renderIcon(func: DG.Func): HTMLElement {
    if (func.options['icon'] != null && (func.options['icon'].startsWith('http://') || func.options['icon'].startsWith('https://')))
      return ui.iconImage('model-icon', func.options['icon']);

    if (func instanceof DG.Script)
      return this.getLanguageIcon(func.language);

    func = DG.Func.find({package: func.package.name, name: func.name})[0];
    let iconUrl = func.package.getIconUrl();
    if (func.options['icon'] != null) {
      const packagePathSegments = iconUrl.split('/');
      packagePathSegments.pop();
      packagePathSegments.push(func.options['icon']);
      iconUrl = packagePathSegments.join('/');
    }
    return ui.iconImage(func.package.name, iconUrl);
  }

  override renderView(x: DG.Func) {
    return this.renderPreview(x).root;
  }

  override renderPreview(x: DG.Func): DG.View {
    const v = super.renderPreview(x);
    v.name = (x.friendlyName ?? x.name) + ' description';
    return DG.View.fromViewAsync(async () => {
      const help = await getContextHelp(x);
      const userGroups = await this.awaitUserGroups();
      const missingMandatoryGroups = ModelHandler.getMissingGroups(x, userGroups);
      const startBtnDiv = missingMandatoryGroups.length ?
        ui.div([this.makeMandatoryGroupsInfo(missingMandatoryGroups)], {style: {
          border: '2px solid var(--red-3)',
          padding: '10px',
          maxWidth: '400px',
        }}) :
        ui.div([ui.bigButton('Launch', () => x.prepare().edit())]);
      v.append(
        ui.div(
          [
            ui.h1(x.friendlyName ?? x.name),
            ui.divText(x.description ?? x.name, {style: {marginBottom: '10px'}}),
            startBtnDiv,
            ui.div([], {style: {padding: '20px '}}),
            ui.markdown(help ?? ''),
          ],
          {
            style: {
              display: 'flex',
              flexDirection: 'column',
              height: '100%',
              overflow: 'auto',
            },
          }),
      );
      return v;
    });
  }

  override renderProperties(func: DG.Func) {
    const a = ui.accordion('ComputeModel');
    a.context = func;
    const titleDiv = ui.div([
      ui.span([this.renderIcon(func), ui.label(func.friendlyName), ui.contextActions(func), ui.star(func.id)])]);
    a.addTitle(titleDiv);

    if (func.description != null)
      titleDiv.appendChild(ui.div([ui.markdown(func.description)], 'model-catalog-description'));
    titleDiv.appendChild(ui.tags(func));
    if ( func.options['dev.status'] != null)
      titleDiv.appendChild(ui.tableFromMap({Status: func.options['dev.status']}));

    a.end();
    return a.root;
  }

  override renderTooltip(x: DG.Func) {
    const h = ui.span([this.renderIcon(x), ui.label(x.friendlyName)]);
    h.classList.add('d4-link-label');
    const card = ui.bind(x, ui.divV([
      ui.h2(h),
      ui.divV([ui.divText(x.description, 'ui-description')], {style: {lineHeight: '150%'}}),
    ]));
    return card;
  }

  override renderCard(x: DG.Func): HTMLElement {
    const h = this.renderMarkup(x);
    h.classList.add('d4-link-label');
    const card = ui.bind(x, ui.h2(h));
    return card;
  }

  override init(): void {
    this.registerParamFunc('Choose version', (func: DG.Func) => {
      ModelHandler.openVersionSelectDialog(DG.toJs(func));
    });

    this.registerParamFunc('Help', (func: DG.Func) => {
      ModelHandler.openHelp(DG.toJs(func));
    });

    setTimeout(async () => {
      const userGroups = await ModelHandler.getUserGroups();
      this.userGroups.next(userGroups);
    });
  }

  private makeMandatoryGroupsInfo(missingMandatoryGroups: {name: string, help?: string}[]) {
    const mandatoryGroupsInfo = ui.div(ui.divV([
      ui.label('You should be a member of the following group(s):', {style: {marginLeft: '0px'}}),
      ...missingMandatoryGroups.map((group) => ui.divV([
        ui.span([getBulletIcon(), group.name], {style: {fontWeight: '600'}}),
        ...group.help ? [ui.span([group.help], {style: {marginLeft: '16px'}})]: [],
        ui.link(`Request group membership`, async () => {
          await requestMembership(group.name);
          hidePopover(mandatoryGroupsInfo);
        }, undefined, {style: {marginLeft: '16px'}}),
      ])),
    ], {style: {gap: '10px'}})) as HTMLElement;
    return mandatoryGroupsInfo;
  }

  private async awaitUserGroups() {
    return this.userGroups.pipe(
      filter((g) => !!g),
      map((g) => g!),
      take(1),
    ).toPromise();
  }
}
