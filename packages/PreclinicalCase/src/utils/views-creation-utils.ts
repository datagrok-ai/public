import * as DG from 'datagrok-api/dg';
import {MATRIX_TABLE_VIEW_NAME, MEASUREMENT_PROFILE_TABLE_VIEW_NAME,
  MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME,
  VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import {hideValidationColumns} from './utils';
import {TableView} from '../types/types';
import {createValidationView} from '../views/validation-table-view';
import {createMatrixTableView} from '../views/matrix-table-view';
import {createMeasurementProfileTableView} from '../views/measurement-profile-table-view';
import {createMICrossDomainView} from '../views/mi-cross-domain-analysis';
import {awaitCheck} from '@datagrok-libraries/test/src/test';

export type StudyTableViewParams = {
  df: DG.DataFrame;
  onTableViewAddedFunc?: (tv: DG.TableView) => void;
}

export function createTableView(studyId: string, viewName: string, 
  createTableViewFunc: (studyId: string) => StudyTableViewParams): DG.TableView {
  const {df, onTableViewAddedFunc} = createTableViewFunc(studyId);
  const tableView = DG.TableView.create(df, false);
  if (onTableViewAddedFunc)
    onTableViewAddedFunc(tableView as DG.TableView);
  awaitCheck(() => (tableView as DG.TableView).grid !== null, `${viewName} hasn't been added`, 10000)
    .then(() => hideValidationColumns(tableView as DG.TableView))
    .catch(() => {});
  tableView.name = viewName;

  return tableView;
}
