import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package, getBatchByCorporateId, getCompoundByCorporateId } from '../package';
import { excludedScopes, MOLTRACK_APP_PATH, Scope, SEARCH_NODE } from './constants';
import { EntityBaseView } from '../views/registration-entity-base';
import { RegistrationView } from '../views/registration-tab';
import { u2 } from '@datagrok-libraries/utils/src/u2';
import { createSearchExpandableNode, createSearchView } from '../views/search';
import { PropertySchemaView } from '../views/schema-view';

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

export function initBulkRegisterView() {
  const registrationView = new RegistrationView();
  registrationView.view.path = createPath('Bulk');
  registrationView.show();
  return registrationView.view;
}

export async function initSchemaView() {
  const schemaView = new PropertySchemaView();
  await schemaView.init();
  schemaView.view.path = createPath('Schema');
  schemaView.show();
  return schemaView.view;
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

const snakeCaseToTitleCase = (str: string) =>
  str.replace(/(^|_)([a-z])/g, (_, __, c) => ' ' + c.toUpperCase()).trim();

export async function getStatisticsWidget(
  onCountClick: (...args: any[]) => void,
): Promise<HTMLTableElement> {
  const scopes = Object.values(Scope);
  const rows: (string | HTMLElement)[][] = await Promise.all(
    scopes.map(async (entity) => {
      try {
        const df = await grok.functions.call('MolTrack:retrieveEntity', { scope: entity, flatten: true });
        const count = df?.rowCount ?? 0;
        const entityTitle = snakeCaseToTitleCase(entity);

        const link = ui.link(count.toString(), () => {
          const args = [entityTitle, entity];
          onCountClick(...args);
        });

        return [entityTitle, link];
      } catch (e: any) {
        grok.shell.error(`Failed to retrieve ${entity}: ${e.message ?? e}`);
        return [snakeCaseToTitleCase(entity), 'Error'];
      }
    }),
  );

  return ui.table(rows, (row) => row);
}

export function getQuickActionsWidget(): HTMLElement {
  const boltIcon = ui.iconFA('bolt');
  boltIcon.classList.add('moltrack-bolt-icon');

  const quickActionsLabel = ui.label('Quick Actions');
  quickActionsLabel.classList.add('moltrack-quick-label');

  const quickActionsTitle = ui.divH([boltIcon, quickActionsLabel]);
  quickActionsTitle.classList.add('moltrack-quick-actions');

  const linksConfig: { label: string; createView: () => Promise<DG.ViewBase> }[] = [
    {
      label: 'Register compound',
      createView: async () => initRegisterView('Compound', true),
    },
    {
      label: 'Register batch',
      createView: async () => initRegisterView('Batch', true),
    },
    {
      label: 'Register bulk',
      createView: async () => initBulkRegisterView(),
    },
    {
      label: 'Search',
      createView: async () => {
        return await createSearchExpandableNode([SEARCH_NODE], () => getStatisticsWidget(createSearchView));
      },
    },
  ];

  const quickLinks = ui.divV(
    linksConfig.map((cfg) =>
      ui.link(cfg.label, async () => {
        grok.shell.v.close();
        await cfg.createView();
      }),
    ),
  );
  quickLinks.classList.add('moltrack-quick-links');

  const container = ui.divV([quickActionsTitle, quickLinks]);
  container.classList.add('moltrack-quick-container');

  return container;
}
