/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsStyleMap} from '../data-loading-utils/json-loader';
import {isPatternCreatedByCurrentUser} from '../../model/pattern-app/oligo-pattern';
import {USER_STORAGE_KEY} from './const';
import {EventBus} from './event-bus';

type PatternsRecord = {[patternName: string]: string};

export class AppDataManager {
  private static instance: AppDataManager;
  private constructor(eventBus: EventBus) {
    this.patternListManager = new PatternListManager(eventBus);
    this.patternListManager.init();
  }
  private patternListManager: PatternListManager;

  static getInstance(eventBus: EventBus): AppDataManager {
    if (!AppDataManager.instance) {
      AppDataManager.instance = new AppDataManager(eventBus);
    }
    return AppDataManager.instance;
  }

  fetchNucleotideBases(): string[] {
    const nucleotideBases: string[] = Object.keys(axolabsStyleMap);
    return nucleotideBases;
  }

  fetchCurrentUserPatterns(): string[] {
    return this.patternListManager.getCurrentUserPatterns();
  }

  fetchOtherUsersPatterns(): string[] {
    return this.patternListManager.getOtherUsersPatterns();
  }
}

class PatternListManager {
  // private initialized = false;
  private currentUserPatterns: string[] = [];
  private otherUsersPatterns: string[] = [];
  private eventBus: EventBus

  constructor(eventBus: EventBus) {
    this.eventBus = eventBus;
  }

  async init(): Promise<void> {
    const patternsRecords = await this.fetchAllPatterns();
    const categorizedPatterns = await this.categorizePatternsByUserOwnership(patternsRecords);

    this.currentUserPatterns = categorizedPatterns.currentUserPatterns;
    this.otherUsersPatterns = categorizedPatterns.otherUsersPatterns;

    this.eventBus.updatePatternList();
  }

  getCurrentUserPatterns(): string[] {
    return this.currentUserPatterns;
  }

  getOtherUsersPatterns(): string[] {
    return this.otherUsersPatterns;
  }

  async savePatternToUserStorage(patternName: string, pattern: string): Promise<void> {
    // todo: implement
  }

  private async fetchAllPatterns() {
    const patternsData = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false) as PatternsRecord;
    return patternsData;
  }

  private async categorizePatternsByUserOwnership(
    patternsRecords: PatternsRecord
  ): Promise<{ currentUserPatterns: string[], otherUsersPatterns: string[] }> {
    const currentUserPatterns: string[] = [];
    const otherUsersPatterns: string[] = [];

    for (const patternName of Object.keys(patternsRecords)) {
      if (await isPatternCreatedByCurrentUser(patternName))
        otherUsersPatterns.push(patternName);
      else
        currentUserPatterns.push(patternName);
    }

    return { currentUserPatterns, otherUsersPatterns };
  }
}
