/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {STRANDS, STRAND} from './const';
import {PatternConfiguration, StrandType, TerminalType} from './types';
import {NucleotideSequences, PhosphorothioateLinkageFlags, StrandTerminusModifications} from './types';
import {GRAPH_SETTINGS_KEYS as G, LEGEND_SETTINGS_KEYS as L} from './const';

import * as rxjs from 'rxjs';
import {map, debounceTime} from 'rxjs/operators';
import {PatternDefaultsProvider} from './defaults-provider';

export class EventBus {
  private _patternName$: rxjs.BehaviorSubject<string>;

  // todo: redundant, remove
  private _isAntisenseStrandActive$: rxjs.BehaviorSubject<boolean>;

  private _nucleotideSequences$: rxjs.BehaviorSubject<NucleotideSequences>;
  private _phosphorothioateLinkageFlags: rxjs.BehaviorSubject<PhosphorothioateLinkageFlags>;
  private _terminalModifications: rxjs.BehaviorSubject<StrandTerminusModifications>;
  private _comment$: rxjs.BehaviorSubject<string>;
  private _modificationsWithNumericLabels$: rxjs.BehaviorSubject<string[]>;
  private _sequenceBase$: rxjs.BehaviorSubject<string>;

  private _patternListUpdated$ = new rxjs.Subject<void>();
  private _patternLoadRequested$ = new rxjs.Subject<string>();
  private _patternSaveRequested$ = new rxjs.Subject<void>();
  private _sequenceBaseReplacementRequested$ = new rxjs.Subject<string>();
  private _uniqueNucleotideBases$ = new rxjs.BehaviorSubject<string[]>([]);

  private _patternDeletionRequested$ = new rxjs.Subject<string>();
  private _tableSelection$ = new rxjs.BehaviorSubject<DG.DataFrame | null>(null);

  private _svgSaveRequested$ = new rxjs.Subject<void>();

  constructor(defaults: PatternDefaultsProvider) {
    this.initializeDefaultState(defaults);

    this._sequenceBaseReplacementRequested$.subscribe((newBase: string) => {
      const oldNucleotideSequences = this._nucleotideSequences$.getValue();
      const newNucleotideSequences = {} as NucleotideSequences;
      STRANDS.forEach((strand) => {
        newNucleotideSequences[strand] = oldNucleotideSequences[strand].map((_: string) => newBase);
      });
    });

    this._nucleotideSequences$.subscribe(() => {
      this.updateUniqueNucleotideBases();
    });
  }

  private updateUniqueNucleotideBases(): void {
    const nucleotideSequences = this._nucleotideSequences$.getValue();
    const uniqueNucleotideBases = new Set<string>();
    STRANDS.forEach((strand) => {
      nucleotideSequences[strand].forEach((nucleotide: string) => {
        uniqueNucleotideBases.add(nucleotide);
      });
    });
    this._uniqueNucleotideBases$.next(Array.from(uniqueNucleotideBases).sort());
  }

  get nucleotideSequencesChanged$(): rxjs.Observable<NucleotideSequences> {
    return this._nucleotideSequences$.asObservable();
  }

  private initializeDefaultState(defaults: PatternDefaultsProvider) {
    this._patternName$ = new rxjs.BehaviorSubject(defaults.getPatternName());
    this._isAntisenseStrandActive$ = new rxjs.BehaviorSubject(defaults.getAntiSenseStrandVisibilityFlag());
    this._nucleotideSequences$ = new rxjs.BehaviorSubject(defaults.getNucleotideSequences());
    this._phosphorothioateLinkageFlags = new rxjs.BehaviorSubject(defaults.getPhosphorothioateLinkageFlags());
    this._terminalModifications = new rxjs.BehaviorSubject(defaults.getTerminusModifications());
    this._comment$ = new rxjs.BehaviorSubject(defaults.getComment());
    this._modificationsWithNumericLabels$ = new rxjs.BehaviorSubject(defaults.getModificationsWithNumericLabels());
    this._sequenceBase$ = new rxjs.BehaviorSubject(defaults.fetchDefaultNucleobase());
  }

  getPatternName(): string {
    return this._patternName$.getValue();
  }

  updatePatternName(patternName: string) {
    this._patternName$.next(patternName);
  }

  get antisenseStrandToggled$(): rxjs.Observable<boolean> {
    return this._isAntisenseStrandActive$.asObservable();
  }

  toggleAntisenseStrand(isActive: boolean) {
    if (!isActive)
      this.updateStrandLength(STRAND.ANTISENSE, 0)
    else
      this.updateStrandLength(STRAND.ANTISENSE, this.getNucleotideSequences()[STRAND.SENSE].length);

    this._isAntisenseStrandActive$.next(isActive);
  }

  isAntiSenseStrandVisible(): boolean {
    return this._isAntisenseStrandActive$.getValue();
  }

  getNucleotideSequences(): NucleotideSequences {
    return this._nucleotideSequences$.getValue();
  }

  updateNucleotideSequences(nucleotideSequences: NucleotideSequences) {
    this._nucleotideSequences$.next(nucleotideSequences);
  }

  updateStrandLength(strand: StrandType, newStrandLength: number): void {
    const sequence = this.getNucleotideSequences()[strand];
    if (sequence.length === newStrandLength) return;

    const phosphorothioateLinkageFlags = this.getPhosphorothioateLinkageFlags()[strand];

    if (newStrandLength === 0) {
      this.updateNucleotideSequences({...this.getNucleotideSequences(), [strand]: []});
      this.updatePhosphorothioateLinkageFlags({
        ...this.getPhosphorothioateLinkageFlags(),
        [strand]: []
      });
      return;
    }

    if (sequence.length > newStrandLength) {
      const newSequence = sequence.slice(0, newStrandLength);
      const newFlags = phosphorothioateLinkageFlags.slice(0, newStrandLength + 1);
      this.updateNucleotideSequences({...this.getNucleotideSequences(), [strand]: newSequence});
      this.updatePhosphorothioateLinkageFlags({
        ...this.getPhosphorothioateLinkageFlags(),
        [strand]: newFlags
      });
      return;
    }

    const appendedNucleotidesLength = newStrandLength - sequence.length;
    const newSequence = sequence.concat(new Array(newStrandLength - sequence.length).fill(this._sequenceBase$.getValue()));
    const appendedFlagsLength = (sequence.length === 0) ? newStrandLength + 1 : appendedNucleotidesLength;
    const newFlags = phosphorothioateLinkageFlags.concat(
      new Array(appendedFlagsLength).fill(true)
    );
    this.updateNucleotideSequences({
      ...this.getNucleotideSequences(),
      [strand]: newSequence
    });
    this.updatePhosphorothioateLinkageFlags({
      ...this.getPhosphorothioateLinkageFlags(),
      [strand]: newFlags
    });
  }

  getPhosphorothioateLinkageFlags(): PhosphorothioateLinkageFlags {
    return this._phosphorothioateLinkageFlags.getValue();
  }

  updatePhosphorothioateLinkageFlags(phosphorothioateLinkageFlags: PhosphorothioateLinkageFlags) {
    this._phosphorothioateLinkageFlags.next(phosphorothioateLinkageFlags);
  }

  get phosphorothioateLingeFlagsChanged$(): rxjs.Observable<PhosphorothioateLinkageFlags> {
    return this._phosphorothioateLinkageFlags.asObservable();
  }

  getTerminalModifications(): StrandTerminusModifications {
    return this._terminalModifications.getValue();
  }

  updateTerminalModifications(terminalModifications: StrandTerminusModifications) {
    this._terminalModifications.next(terminalModifications);
  }

  getComment(): string {
    return this._comment$.getValue();
  }

  updateComment(comment: string) {
    this._comment$.next(comment);
  }

  getModificationsWithNumericLabels(): string[] {
    return this._modificationsWithNumericLabels$.getValue();
  }

  updateModificationsWithNumericLabels(modificationsWithNumericLabels: string[]) {
    const newValue = Array.from(new Set(modificationsWithNumericLabels)).sort();
    this._modificationsWithNumericLabels$.next(newValue);
  }

  changeSequenceBase(base: string) {
    this._sequenceBase$.next(base);
  }

  get patternLoadRequested$(): rxjs.Observable<string> {
    return this._patternLoadRequested$.asObservable();
  }

  requestPatternLoad(patternName: string) {
    this._patternLoadRequested$.next(patternName);
  }

  get patternListUpdated$(): rxjs.Observable<void> {
    return this._patternListUpdated$.asObservable();
  }

  updatePatternList() {
    this._patternListUpdated$.next();
  }

  get tableSelectionChanged$(): rxjs.Observable<DG.DataFrame | null> {
    return this._tableSelection$.asObservable();
  }

  selectTable(table: DG.DataFrame | null) {
    this._tableSelection$.next(table);
  }

  getTableSelection(): DG.DataFrame | null {
    return this._tableSelection$.getValue();
  }

  deletePattern(patternName: string) {
    this._patternDeletionRequested$.next(patternName);
  }

  requestPatternSave() {
    this._patternSaveRequested$.next();
  }

  patternSaveRequested$(): rxjs.Observable<void> {
    return this._patternSaveRequested$.asObservable();
  }

  replaceSequenceBase(newNucleobase: string) {
    this._sequenceBaseReplacementRequested$.next(newNucleobase);
  }

  get patternStateChanged$(): rxjs.Observable<void> {
    const observable = rxjs.merge(
      this._patternName$,
      this._isAntisenseStrandActive$,
      this._nucleotideSequences$,
      this._phosphorothioateLinkageFlags,
      this._terminalModifications,
      this._comment$,
      this._modificationsWithNumericLabels$,
    ) as rxjs.Observable<void>;

    return observable.pipe(debounceTime(50));
  }

  getSequenceBase(): string {
    return this._sequenceBase$.getValue();
  }

  uniqueNucleotideBasesChanged$(): rxjs.Observable<string[]> {
    return this._uniqueNucleotideBases$.asObservable();
  }

  getUniqueNucleotideBases(): string[] {
    return this._uniqueNucleotideBases$.getValue();
  }

  get svgSaveRequested$(): rxjs.Observable<void> {
    return this._svgSaveRequested$.asObservable();
  }

  requestSvgSave() {
    this._svgSaveRequested$.next();
  }
  
  setAllPTOLinkages(value: boolean) {
    const flags = this.getPhosphorothioateLinkageFlags();
    STRANDS.forEach((strand) => {
      flags[strand] = flags[strand].map(() => value);
    });
    this.updatePhosphorothioateLinkageFlags(flags);
  }

  setPatternConfig(patternConfiguration: PatternConfiguration) {
    this._patternName$.next(patternConfiguration[L.PATTERN_NAME]);
    this._isAntisenseStrandActive$.next(patternConfiguration[G.IS_ANTISENSE_STRAND_INCLUDED]);
    this._nucleotideSequences$.next(patternConfiguration[G.NUCLEOTIDE_SEQUENCES]);
    this._phosphorothioateLinkageFlags.next(patternConfiguration[G.PHOSPHOROTHIOATE_LINKAGE_FLAGS]);
    this._terminalModifications.next(patternConfiguration[G.STRAND_TERMINUS_MODIFICATIONS]);
    this._comment$.next(patternConfiguration[L.PATTERN_COMMENT]);
    this._modificationsWithNumericLabels$.next(patternConfiguration[L.NUCLEOTIDES_WITH_NUMERIC_LABELS]);
  }

  getPatternConfig(): PatternConfiguration {
    return {
      [L.PATTERN_NAME]: this.getPatternName(),
      [G.IS_ANTISENSE_STRAND_INCLUDED]: this.isAntiSenseStrandVisible(),
      [G.NUCLEOTIDE_SEQUENCES]: this.getNucleotideSequences(),
      [G.PHOSPHOROTHIOATE_LINKAGE_FLAGS]: this.getPhosphorothioateLinkageFlags(),
      [G.STRAND_TERMINUS_MODIFICATIONS]: this.getTerminalModifications(),
      [L.PATTERN_COMMENT]: this.getComment(),
      [L.NUCLEOTIDES_WITH_NUMERIC_LABELS]: this.getModificationsWithNumericLabels(),
    };
  }

  setPhosphorothioateLinkageFlag(strand: StrandType, index: number, newValue: boolean) {
    const flags = this.getPhosphorothioateLinkageFlags();
    flags[strand][index] = newValue;
    this.updatePhosphorothioateLinkageFlags(flags);
  }

  setNucleotideBase(strand: StrandType, index: number, value: string) {
    const sequences = this.getNucleotideSequences();
    sequences[strand][index] = value;
    const labelledModifications = this.getModificationsWithNumericLabels();
    this.updateModificationsWithNumericLabels(labelledModifications.concat(value));
    this.updateNucleotideSequences(sequences);
  }

  get updatePatternEditor$(): rxjs.Observable<void> {
    return rxjs.merge(
      this._isAntisenseStrandActive$.asObservable().pipe(map(() => {})),
      this._nucleotideSequences$.asObservable().pipe(map(() => {})),
    ).pipe(debounceTime(50)) as rxjs.Observable<void>;
  }
}
