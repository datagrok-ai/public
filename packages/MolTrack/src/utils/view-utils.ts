import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package, getBatchByCorporateId, getCompoundByCorporateId } from '../package';
import { excludedScopes, MOLTRACK_APP_PATH, Scope } from './constants';
import { EntityBaseView } from '../views/registration-entity-base';
import { RegistrationView } from '../views/registration-tab';
import { u2 } from '@datagrok-libraries/utils/src/u2';

export function createPath(viewName: string) {
  let path = `${MOLTRACK_APP_PATH}/`;
  path += encodeURIComponent(viewName);
  return path;
}

export function createPathFromArr(view: DG.ViewBase, pathComponents: string[]) {
  let path = view.path.toLowerCase().includes(MOLTRACK_APP_PATH.toLowerCase()) ? `` : `${MOLTRACK_APP_PATH}`;
  path += `/${pathComponents.map((it) => encodeURIComponent(it)).join('/')}`;
  return path;
}

export async function buildRegistrationView({
  title,
  smiles,
  pathQueryParam,
  propsList,
  batchSection,
  path,
}: {
  title: string;
  smiles: string;
  pathQueryParam: string;
  propsList: any[];
  batchSection: boolean,
  path: string,
}) {
  const registrationView = new EntityBaseView(false);
  registrationView.initialSmiles = smiles;
  registrationView.singleRetrieved = true;
  registrationView.title = title;
  registrationView.isBatchSectionExpanded = batchSection;
  await registrationView.buildUIMethod();

  registrationView.formBackingObject = {};
  for (const prop of propsList) {
    const value = prop.value_string ?? prop.value_num ?? prop.value_datetime ?? prop.value_uuid;
    if (value !== null && value !== undefined)
      registrationView.formBackingObject[prop.name] = value;
  }

  for (const input of registrationView.inputs) {
    const propName = input.property.name;
    if (registrationView.formBackingObject.hasOwnProperty(propName))
      input.value = registrationView.formBackingObject[propName];
  }

  registrationView.show();
  registrationView.view.path = `${createPath(path)}?${pathQueryParam}`;

  return registrationView.view;
}

export async function compoundView(corporateCompoundId: string) {
  const compound = await getCompoundByCorporateId(corporateCompoundId);
  const { canonical_smiles: smiles, properties = [] } = compound;
  return buildRegistrationView({
    title: `Compound: ${corporateCompoundId}`,
    smiles,
    pathQueryParam: `corporate_compound_id=${encodeURIComponent(corporateCompoundId)}`,
    propsList: properties,
    batchSection: false,
    path: 'Compound',
  });
}

export async function batchView(corporateBatchId: string) {
  const batch = await getBatchByCorporateId(corporateBatchId);
  const { canonical_smiles: smiles, properties: compoundProps } = batch.compound;
  const combinedPropsMap = new Map(
    [...compoundProps, ...batch.properties].map((prop: any) => [prop.name, prop]),
  );


  return buildRegistrationView({
    title: `Batch: ${corporateBatchId}`,
    smiles,
    pathQueryParam: `corporate_batch_id=${encodeURIComponent(corporateBatchId)}`,
    propsList: Array.from(combinedPropsMap.values()),
    batchSection: true,
    path: 'Batch',
  });
}

export function initRegisterView(entity: 'Compound' | 'Batch', setPath: boolean = true) {
  const isBatch = entity === 'Batch';
  const view = new EntityBaseView(!isBatch);

  if (isBatch) {
    view.title = `Register a new ${entity.toLowerCase()}`;
    view.isBatchSectionExpanded = true;
    view.path = entity;
    view.buildUIMethod();
  }

  view.view.name = `Register a ${entity.toLowerCase()}`;
  if (setPath) view.view.path = createPath(entity);

  view.show();
  return view.view;
}

export function initBulkRegisterView(setPath: boolean = true) {
  const registrationView = new RegistrationView();
  registrationView.view.path = createPath('Bulk');
  registrationView.show();
  return registrationView.view;
}
export function getAppHeader(): HTMLElement {
  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/moltrack.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/MolTrack/README.md',
    description: '- Chemical compound registration system\n' +
        '- Analyze assay data\n' +
        '- Find contextual information on molecules.\n',
    appTitle: 'MolTrack',
    appSubTitle: 'Track, analyze, and manage chemical data',
    bottomLine: true,
  });
  return appHeader;
}

export async function getStatisticsWidget(onCountClick: (...args: any[]) => void, isSearch?: boolean):
Promise<HTMLTableElement> {
  let scopes = Object.values(Scope);
  if (isSearch)
    scopes = scopes.filter((scope) => !excludedScopes.includes(scope));
  const rows: any[][] = await Promise.all(
    scopes.map(async (entity) => {
      try {
        const df = await grok.functions.call('MolTrack:retrieveEntity', { scope: entity, flatten: true });
        const count = df?.rowCount ?? 0;

        return [
          entity,
          ui.link(count.toString(), () => {
            const funcParams = isSearch ? [entity.charAt(0).toUpperCase() + entity.slice(1), entity] : [df];
            onCountClick(...funcParams);
          }),
        ];
      } catch (e) {
        grok.shell.error(`Failed to retrieve ${entity}: ${e}`);
        return [entity, 'Error'];
      }
    }),
  );

  return ui.table(rows, (row) => row);
};
