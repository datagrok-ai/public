import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';

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

  getValidFilesPathList(): string[] {
    return this._libraryFilesUpdateSubject$.getValue();
  }

  changeValidFilesPathList(newList: string[]): void {
    this._libraryFilesUpdateSubject$.next(newList);
  }

  get updateUIControlsRequested$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdateSubject$.pipe(
      debounceTime(1000)
    );
  }

  get updateValidLibraryFileListRequested$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdateSubject$.pipe(
      debounceTime(3000)
    );
  }
}
