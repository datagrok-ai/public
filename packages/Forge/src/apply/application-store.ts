import * as rxjs from 'rxjs';
import {ApplicationInsert, forgeDb} from '../generated/db';

/** The id of the model an application was just recorded for, whatever its status. */
export const applicationRecorded = new rxjs.Subject<string>();

export async function recordApplication(record: ApplicationInsert): Promise<void> {
  await forgeDb.applications.insert(record);
  applicationRecorded.next(record.model_id);
}
