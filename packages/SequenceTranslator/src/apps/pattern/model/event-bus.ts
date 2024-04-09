/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime, map, skip, switchMap} from 'rxjs/operators';

import {GRAPH_SETTINGS_KEYS as G, LEGEND_SETTINGS_KEYS as L, PATTERN_RECORD_KEYS as R, STRAND, STRANDS} from './const';
// import {PatternDefaultsProvider} from './defaults-provider';
import {
  NucleotideSequences, PatternConfigRecord, PatternConfiguration,
  PhosphorothioateLinkageFlags, StrandTerminusModifications, StrandType
} from './types';
import {
  getMostFrequentNucleotide, getUniqueNucleotides, getUniqueNucleotidesWithNumericLabels, StrandEditingUtils
} from './utils';
import {DataManager} from './data-manager';

/** Manager of all events in the application, *the* central state manager.
 * Use for communication between app's components to avoid tight coupling. */
export class EventBus {
  private _patternAuthorSelection$: rxjs.BehaviorSubject<string>;

  private _patternName$: rxjs.BehaviorSubject<string>;

  private _isAntisenseStrandActive$: rxjs.BehaviorSubject<boolean>;

  private _nucleotideSequences$: rxjs.BehaviorSubject<NucleotideSequences>;
  private _phosphorothioateLinkageFlags: rxjs.BehaviorSubject<PhosphorothioateLinkageFlags>;
  private _terminalModifications: rxjs.BehaviorSubject<StrandTerminusModifications>;
  private _comment$: rxjs.BehaviorSubject<string>;
  private _modificationsWithNumericLabels$: rxjs.BehaviorSubject<string[]>;
  private _sequenceBase$: rxjs.BehaviorSubject<string>;

  private _patternListUpdated$ = new rxjs.Subject<void>();

  private _patternLoadRequested$ = new rxjs.Subject<string>();
  private _patternLoaded$ = new rxjs.Subject<string>();
  private _uniqueNucleotides$ = new rxjs.BehaviorSubject<string[]>([]);

  private _patternDeletionRequested$ = new rxjs.Subject<string>();
  private _tableSelection$ = new rxjs.BehaviorSubject<DG.DataFrame | null>(null);

  private _svgSaveRequested$ = new rxjs.Subject<void>();

  constructor(
    private dataManager: DataManager,
    initialPaternConfigRecord: PatternConfigRecord
  ) {
    this.initializeAuthorSelection(initialPaternConfigRecord);
    this.initializePatternState(initialPaternConfigRecord);

    this._nucleotideSequences$.subscribe(() => {
      this.updateUniqueNucleotides();
      this.updateSequenceBase();
    });
  }

  private updateUniqueNucleotides(): void {
    const sequences = this._nucleotideSequences$.getValue();
    const uniqueNucleotides = getUniqueNucleotides(sequences);
    this._uniqueNucleotides$.next(uniqueNucleotides);
  }

  private updateSequenceBase(): void {
    const nucleotideSequences = this._nucleotideSequences$.getValue();
    const mostFrequentNucleotide = getMostFrequentNucleotide(nucleotideSequences);
    this._sequenceBase$.next(mostFrequentNucleotide);
  }

  get nucleotideSequencesChanged$(): rxjs.Observable<NucleotideSequences> {
    return this._nucleotideSequences$.asObservable();
  }

  private initializeAuthorSelection(
    initialConfigRecord: PatternConfigRecord
  ) {
    const patternAuthorId = initialConfigRecord[R.AUTHOR_ID];
    if (this.dataManager.isCurrentUserId(patternAuthorId))
      this._patternAuthorSelection$ = new rxjs.BehaviorSubject(this.dataManager.getCurrentUserAuthorshipCategory());
    else
      this._patternAuthorSelection$ = new rxjs.BehaviorSubject(this.dataManager.getOtherUsersAuthorshipCategory());
  }

  private initializePatternState(
    initialConfigRecord: PatternConfigRecord
  ) {
    const initialPattern = initialConfigRecord[R.PATTERN_CONFIG];
    this._patternName$ = new rxjs.BehaviorSubject(initialPattern[L.PATTERN_NAME]);
    this._isAntisenseStrandActive$ = new rxjs.BehaviorSubject(
      initialPattern[G.IS_ANTISENSE_STRAND_INCLUDED]
    );
    this._nucleotideSequences$ = new rxjs.BehaviorSubject(
      initialPattern[G.NUCLEOTIDE_SEQUENCES]
    );
    this._phosphorothioateLinkageFlags = new rxjs.BehaviorSubject(
      initialPattern[G.PHOSPHOROTHIOATE_LINKAGE_FLAGS]
    );
    this._terminalModifications = new rxjs.BehaviorSubject(
      initialPattern[G.STRAND_TERMINUS_MODIFICATIONS]
    );
    this._comment$ = new rxjs.BehaviorSubject(
      initialPattern[L.PATTERN_COMMENT]
    );
    this._modificationsWithNumericLabels$ = new rxjs.BehaviorSubject(
      initialPattern[L.NUCLEOTIDES_WITH_NUMERIC_LABELS]
    );
    this._sequenceBase$ = new rxjs.BehaviorSubject(
      getMostFrequentNucleotide(initialPattern[G.NUCLEOTIDE_SEQUENCES])
    );
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

  isAntisenseStrandActive(): boolean {
    return this._isAntisenseStrandActive$.getValue();
  }

  toggleAntisenseStrand(isActive: boolean) {
    if (!isActive)
      this.updateStrandLength(STRAND.ANTISENSE, 0);
    else
      this.updateStrandLength(STRAND.ANTISENSE, this.getNucleotideSequences()[STRAND.SENSE].length);

    this._isAntisenseStrandActive$.next(isActive);
  }

  isAntiSenseStrandActive(): boolean {
    return this._isAntisenseStrandActive$.getValue();
  }

  getNucleotideSequences(): NucleotideSequences {
    return this._nucleotideSequences$.getValue();
  }

  updateNucleotideSequences(nucleotideSequences: NucleotideSequences) {
    this._nucleotideSequences$.next(nucleotideSequences);
  }

  updateStrandLength(strand: StrandType, newStrandLength: number): void {
    const originalNucleotides = this.getNucleotideSequences()[strand];
    if (originalNucleotides.length === newStrandLength) return;

    const originalPTOFlags = this.getPhosphorothioateLinkageFlags()[strand];

    if (newStrandLength === 0) {
      this.setNewStrandData([], [], strand);
      return;
    }

    if (originalNucleotides.length > newStrandLength) {
      const {nucleotides, ptoFlags} = StrandEditingUtils.getTruncatedStrandData(
        originalNucleotides, originalPTOFlags, newStrandLength
      );
      this.setNewStrandData(nucleotides, ptoFlags, strand);
      return;
    }

    const sequenceBase = this.getSequenceBase();
    const {nucleotides, ptoFlags} = StrandEditingUtils.getExtendedStrandData(
      originalNucleotides, originalPTOFlags, newStrandLength, sequenceBase
    );
    this.setNewStrandData(nucleotides, ptoFlags, strand);
  }

  private setNewStrandData(
    newNucleotides: string[],
    newPTOFlags: boolean[],
    strand: StrandType
  ): void {
    this.updateNucleotideSequences({
      ...this.getNucleotideSequences(),
      [strand]: newNucleotides
    });
    this.updatePhosphorothioateLinkageFlags({
      ...this.getPhosphorothioateLinkageFlags(),
      [strand]: newPTOFlags
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
    const newValue = getUniqueNucleotidesWithNumericLabels(modificationsWithNumericLabels);
    this._modificationsWithNumericLabels$.next(newValue);
  }

  get patternLoadRequested$(): rxjs.Observable<string> {
    return this._patternLoadRequested$.asObservable();
  }

  requestPatternLoad(patternHash: string) {
    this._patternLoadRequested$.next(patternHash);
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

  requestPatternDeletion(patternName: string) {
    this._patternDeletionRequested$.next(patternName);
  }

  get patternDeletionRequested$(): rxjs.Observable<string> {
    return this._patternDeletionRequested$.asObservable();
  }

  replaceSequenceBase(newNucleobase: string) {
    const oldNucleotideSequences = this._nucleotideSequences$.getValue();
    const newNucleotideSequences = {} as NucleotideSequences;
    STRANDS.forEach((strand) => {
      newNucleotideSequences[strand] = oldNucleotideSequences[strand].map(() => newNucleobase);
    });
    this._nucleotideSequences$.next(newNucleotideSequences);

    const labelledNucleotides = this._modificationsWithNumericLabels$.getValue();
    if (!labelledNucleotides.includes(newNucleobase))
      this.updateModificationsWithNumericLabels(labelledNucleotides.concat(newNucleobase));
  }

  get patternStateChanged$(): rxjs.Observable<void> {
    const observable = rxjs.merge(
      this._patternName$.pipe(debounceTime(300), map(() => {})),
      this._isAntisenseStrandActive$,
      this._nucleotideSequences$,
      this._phosphorothioateLinkageFlags,
      this._terminalModifications,
      this._comment$.pipe(debounceTime(300)),
      this._modificationsWithNumericLabels$
    ) as rxjs.Observable<void>;

    return observable;
  }

  getSequenceBase(): string {
    return this._sequenceBase$.getValue();
  }

  uniqueNucleotidesChanged$(): rxjs.Observable<string[]> {
    // WARNING: switchMap is necessary to preserve order of events
    const observable = this.patternStateChanged$.pipe(
      switchMap(() => this._uniqueNucleotides$)
    );

    return observable;
  }

  getUniqueNucleotides(): string[] {
    return this._uniqueNucleotides$.getValue();
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
      [G.IS_ANTISENSE_STRAND_INCLUDED]: this.isAntiSenseStrandActive(),
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

  setNucleotide(strand: StrandType, index: number, value: string) {
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
      this._patternLoaded$.asObservable().pipe(map(() => {}))
    ) as rxjs.Observable<void>;
  }

  updateControlsUponPatternLoaded(patternHash: string) {
    this._patternLoaded$.next(patternHash);
  }

  get patternLoaded$(): rxjs.Observable<string> {
    return this._patternLoaded$.asObservable();
  }

  get userSelection$(): rxjs.Observable<string> {
    return this._patternAuthorSelection$.asObservable().pipe(skip(1));
  }

  selectUser(username: string) {
    this._patternAuthorSelection$.next(username);
  }

  getSelectedUser(): string {
    return this._patternAuthorSelection$.getValue();
  }
}

