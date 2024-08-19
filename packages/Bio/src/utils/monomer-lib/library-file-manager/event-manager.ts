import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BehaviorSubject, Observable, Subject} from 'rxjs';
import {debounceTime, tap, skip} from 'rxjs/operators';

import {IMonomerLibFileEventManager} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

export class MonomerLibFileEventManager implements IMonomerLibFileEventManager {
  // WARNING: this must be a singleton because it manages the unique state
  private constructor() {}

  private static _instance: MonomerLibFileEventManager;

  static getInstance(): MonomerLibFileEventManager {
    if (!MonomerLibFileEventManager._instance)
      MonomerLibFileEventManager._instance = new MonomerLibFileEventManager();

    return MonomerLibFileEventManager._instance;
  }

  private _libFilesUpdateSubject$ = new BehaviorSubject<string[]>([]);
  private _setFilesUpdateSubject$ = new BehaviorSubject<string[]>([]);

  private _addLibraryFilesSubject$ = new Subject<void>();
  private _librarySelectionSubject$ = new Subject<[string, boolean]>();

  getValidLibPathList(): string[] {
    return this._libFilesUpdateSubject$.getValue();
  }

  getValidSetPathList(): string[] {
    return this._setFilesUpdateSubject$.getValue();
  }

  // TODO: remove after adding init from user data storage
  // WARNING: a temporary solution
  async getValidLibraryPathsAsynchronously(): Promise<string[]> {
    return new Promise((resolve) => {
      const sub = this._libFilesUpdateSubject$.pipe<string[]>(
        skip(1)
      ).subscribe((fileNames) => {
        resolve(fileNames);
        sub.unsubscribe();
      });
    });
  }

  changeValidLibPathList(newLibPathList: string[]): void {
    this._libFilesUpdateSubject$.next(newLibPathList);
  }

  changeValidSetPathList(newSetPathList: string[]): void {
    this._setFilesUpdateSubject$.next(newSetPathList);
  }

  get updateUIControlsRequested$(): Observable<string[]> {
    return this._libFilesUpdateSubject$.pipe(
      // debounceTime(1000)
    );
  }

  get updateValidLibraryFileListRequested$(): Observable<string[]> {
    return this._libFilesUpdateSubject$.pipe(
      // debounceTime(3000)
    );
  }

  get updateValidSetFileListRequested$(): Observable<string[]> {
    return this._setFilesUpdateSubject$.pipe(
      // debounceTime(1000)
    );
  }

  get addLibraryFileRequested$(): Observable<void> {
    return this._addLibraryFilesSubject$.pipe(
      // debounceTime(1000)
    );
  }

  addLibraryFile(): void {
    this._addLibraryFilesSubject$.next();
  }

  get librarySelectionRequested$(): Observable<[string, boolean]> {
    return this._librarySelectionSubject$;
  }

  updateLibrarySelectionStatus(fileName: string, isSelected: boolean): void {
    this._librarySelectionSubject$.next([fileName, isSelected]);
  }
}
