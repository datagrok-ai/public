import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {openDB, DBSchema, StoreNames} from 'idb';
import {computedAsync} from '@vueuse/core';

export function useLayoutDb<S extends DBSchema>(
  dbName: string,
  storeName: StoreNames<S>,
) {
  const layoutDatabase = computedAsync(async () => {
    const db = await openDB<S>(dbName, 2, {
      blocked: () => {
        grok.shell.error(`Layout database requires update. Please close all webpages with models opened.`);
      },
      upgrade: (db, oldVersion) => {
        if (oldVersion === 0)
          db.createObjectStore(storeName);

        const hasStore = db.objectStoreNames.contains(storeName);
        if (!hasStore)
          db.createObjectStore(storeName);
      },
    });

    return db;
  }, null);
  return {layoutDatabase};
}
