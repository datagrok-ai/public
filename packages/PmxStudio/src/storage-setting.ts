/* The domainStorage package setting (on by default): plans and studies saved as rows of PMX Studio's domain tables
   (store.ts); off, as files in each user's Home. Apart from store.ts so that any module may read it. */
import {_package} from './package';

/** Whether package settings say rows: anything but an explicit false (a setting never saved is the default, on). */
export const storageOnFrom = (settings: {[key: string]: unknown} | null | undefined): boolean =>
  String(settings?.domainStorage) !== 'false';
export const domainStorageOn = (): boolean => storageOnFrom(_package.settings);
