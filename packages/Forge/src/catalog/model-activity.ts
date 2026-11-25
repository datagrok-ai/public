import type {Dayjs} from 'dayjs';
import {ApplicationRow, ApplicationStatus, forgeDb} from '../generated/db';

const MAX_ACTIVITY_ROWS = 100;

export interface ModelActivity {
  /** The newest applications, at most MAX_ACTIVITY_ROWS, newest first. */
  applications: ApplicationRow[];
  /** The number of all the model's applications. */
  count: number;
  /** The newest application of any status; undefined when the model was never applied. */
  lastRun?: {when: Dayjs; status: ApplicationStatus};
}

export async function modelActivity(id: string): Promise<ModelActivity> {
  const applications = await forgeDb.applications.query().where('model_id', '=', id).orderBy('created_on', true)
    .top(MAX_ACTIVITY_ROWS);
  // Below the cap, the list is all of them.
  const count = applications.length < MAX_ACTIVITY_ROWS ? applications.length :
    await forgeDb.applications.query().where('model_id', '=', id).count();
  if (applications.length === 0)
    return {applications, count};
  const {created_on: when, status} = applications[0];
  return {applications, count, lastRun: {when, status}};
}
