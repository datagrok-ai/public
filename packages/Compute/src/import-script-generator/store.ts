import * as rxjs from "rxjs";

export class ImportGeneratorStore {
  scriptName: string = '';
  description: string = '';
  dfSheetsToSeek: string[] = [];
  scalarsToSeek: string[] = []

  scriptText = new rxjs.BehaviorSubject('');
  infoText = new rxjs.BehaviorSubject('');
  lastSavedScriptId = new rxjs.BehaviorSubject('');
}