import * as rxjs from 'rxjs';
import {debounceTime, tap} from 'rxjs/operators';

export class MonomerLibFileEventManager {
  // WARNING: this must be a singleton because it manages the unique state
  private constructor() {}

  private static _instance: MonomerLibFileEventManager;

  static getInstance(): MonomerLibFileEventManager {
    if (!MonomerLibFileEventManager._instance)
      MonomerLibFileEventManager._instance = new MonomerLibFileEventManager();

    return MonomerLibFileEventManager._instance;
  }

  private _libraryFilesUpdateSubject$ = new rxjs.BehaviorSubject<string[]>([]);
  private _addLibraryFilesSubject$ = new rxjs.Subject<void>();
  private _librarySelectionSubject$ = new rxjs.Subject<[string, boolean]>();

  getValidFilesPathList(): string[] {
    return this._libraryFilesUpdateSubject$.getValue();
  }

  changeValidFilesPathList(newList: string[]): void {
    this._libraryFilesUpdateSubject$.next(newList);
  }

  get updateUIControlsRequested$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdateSubject$.pipe(
      // debounceTime(1000)
    );
  }

  get updateValidLibraryFileListRequested$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdateSubject$.pipe(
      // debounceTime(3000)
    );
  }

  get addLibraryFileRequested$(): rxjs.Observable<void> {
    return this._addLibraryFilesSubject$.pipe(
      // debounceTime(1000)
    );
  }

  addLibraryFile(): void {
    this._addLibraryFilesSubject$.next();
  }

  get librarySelectionRequested$(): rxjs.Observable<[string, boolean]> {
    return this._librarySelectionSubject$;
  }

  updateLibrarySelectionStatus(fileName: string, isSelected: boolean): void {
    this._librarySelectionSubject$.next([fileName, isSelected]);
  }
}
