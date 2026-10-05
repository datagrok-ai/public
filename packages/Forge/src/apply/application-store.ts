import {ApplicationInsert, forgeDb} from '../generated/db';

export async function recordApplication(record: ApplicationInsert): Promise<void> {
  await forgeDb.applications.insert(record);
}
