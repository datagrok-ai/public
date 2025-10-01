import { getBatchByCorporateId, getCompoundByCorporateId } from '../package';
import { MOLTRACK_APP_PATH } from './constants';
import { EntityBaseView } from '../views/registration-entity-base';

export function createPath(viewName: string) {
  let path = `${MOLTRACK_APP_PATH}/`;
  path += encodeURIComponent(viewName);
  return path;
}

export function createPathFromArr(pathComponents: string[]) {
  let path = `${MOLTRACK_APP_PATH}`;
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
