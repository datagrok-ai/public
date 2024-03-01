import * as rxjs from 'rxjs';
import {debounceTime, tap, skip} from 'rxjs/operators';

export class MonomerLibFileEventManager {
  // WARNING: this must be a singleton because it manages the unique state
  private constructor() {}

  private static _instance: MonomerLibFileEventManager;

  static getInstance(): MonomerLibFileEventManager {
    if (!MonomerLibFileEventManager._instance)
      MonomerLibFileEventManager._instance = new MonomerLibFileEventManager();

    return MonomerLibFileEventManager._instance;
  }

  private _libraryFilesUpdate$ = new rxjs.BehaviorSubject<string[]>([]);
  private _addLibraryFiles$ = new rxjs.Subject<void>();
  private _librarySelection$ = new rxjs.Subject<[string, boolean]>();
  private _saveLibrarySettings$ = new rxjs.Subject<void>();
  private _resetLibrarySettings$ = new rxjs.Subject<void>();

  getValidFilesPathList(): string[] {
    return this._libraryFilesUpdate$.getValue();
  }

  // TODO: remove after adding init from user data storage
  // WARNING: a temporary solution
  async getValidLibraryPathsAsynchronously(): Promise<string[]> {
    return new Promise((resolve) => {
      this._libraryFilesUpdate$.pipe(
        skip(1)
      ).subscribe((fileNames) => {
        resolve(fileNames);
      });
    });
  }

  changeValidFilesPathList(newList: string[]): void {
    this._libraryFilesUpdate$.next(newList);
  }

  get updateUIControlsRequested$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdate$.pipe(
      // debounceTime(1000)
    );
  }

  get updateValidLibraryFileListRequested$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdate$.pipe(
      // debounceTime(3000)
    );
  }

  get addLibraryFileRequested$(): rxjs.Observable<void> {
    return this._addLibraryFiles$.pipe(
      // debounceTime(1000)
    );
  }

  addLibraryFile(): void {
    this._addLibraryFiles$.next();
  }

  get librarySelectionRequested$(): rxjs.Observable<[string, boolean]> {
    return this._librarySelection$;
  }

  updateLibrarySelectionStatus(fileName: string, isSelected: boolean): void {
    this._librarySelection$.next([fileName, isSelected]);
  }

  saveLibrarySettings(): void {
    this._saveLibrarySettings$.next();
  }

  get saveLibrarySettingsRequested$(): rxjs.Observable<void> {
    return this._saveLibrarySettings$;
  }

  resetLibrarySettings(): void {
    this._resetLibrarySettings$.next();
  }

  get resetLibrarySettingsRequested$(): rxjs.Observable<void> {
    return this._resetLibrarySettings$;
  }
}
