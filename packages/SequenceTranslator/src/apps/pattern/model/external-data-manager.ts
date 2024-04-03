/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';

import {
  STORAGE_NAME, GRAPH_SETTINGS_KEY_LIST as GKL, LEGEND_SETTINGS_KEYS as L, PATTERN_RECORD_KEYS as R, OTHER_USERS
} from './const';
import {PatternConfiguration, PatternExistsError, PatternNameExistsError} from './types';
import {EventBus} from './event-bus';

import objectHash from 'object-hash';
import {DataInitializer} from './data-initializer';

export class PatternAppDataManager {
  // WARNING: not a singleton, encapsulated async init logic in a static method
  private constructor(
    private eventBus: EventBus,
    private patternConfigManager: PatternConfigManager
  ) { }

  static async getInstance(
    eventBus: EventBus,
    dataInitializer: DataInitializer
  ): Promise<PatternAppDataManager> {
    const patternConfigManager = new PatternConfigManager(eventBus, dataInitializer);
    return new PatternAppDataManager(eventBus, patternConfigManager);
  }

  getCurrentUserPatternNames(): string[] {
    return this.patternConfigManager.getCurrentUserPatternNames() || [];
  }

  getOtherUsersPatternNames(): string[] {
    return this.patternConfigManager.getOtherUsersPatternNames() || [];
  }

  getCurrentUserName(): string {
    return this.patternConfigManager.getCurrentUserName();
  }

  async savePatternToUserStorage(): Promise<void> {
    const patternConfig = this.eventBus.getPatternConfig();
    try {
      await this.patternConfigManager.savePatternToUserStorage(patternConfig);
    } catch (e) {
      if (e instanceof PatternNameExistsError || e instanceof PatternExistsError)
        throw e;
      else
        console.error('Error while saving pattern to user storage', e);
    }
  }

  getPatternHash(patternName: string, isCurrentUserPattern: boolean): string {
    return this.patternConfigManager.getPatternHash(patternName, isCurrentUserPattern);
  }

  async getPatternConfig(hash: string): Promise<PatternConfiguration> {
    return await this.patternConfigManager.getPatternConfig(hash);
  }

  getPatternNameByHash(hash: string): string {
    return this.patternConfigManager.getPatternNameByHash(hash);
  }
}

class PatternConfigManager {
  private currentUserName: string;
  private otherUsersPatternNameToHash = new Map<string, string>();
  private currentUserPatternNameToHash = new Map<string, string>();

  constructor(
    private eventBus: EventBus,
    private dataInitializer: DataInitializer
  ) {
    this.currentUserName = this.dataInitializer.getCurrentUserName();
    this.currentUserPatternNameToHash = this.dataInitializer.getCurrentUserPatternNameToHash();
    this.otherUsersPatternNameToHash = this.dataInitializer.getOtherUsersPatternNameToHash();

    this.eventBus.patternDeletionRequested$.subscribe(async (patternName: string) => {
      await this.deletePattern(patternName);
    });
  }

  getCurrentUserPatternNames(): string[] {
    return Array.from(this.currentUserPatternNameToHash.keys())
      .sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
  }

  getOtherUsersPatternNames(): string[] {
    return Array.from(this.otherUsersPatternNameToHash.keys())
      .sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
  }

  getCurrentUserName(): string {
    return this.currentUserName;
  }

  async savePatternToUserStorage(patternConfig: PatternConfiguration): Promise<void> {
    const hash = this.getHash(patternConfig);
    this.validatePatternUniqueness(hash);

    const patternName = patternConfig[L.PATTERN_NAME];
    this.validatePatternNameUniqueness(patternName);

    const record = await this.getRecordFromPattern(patternConfig);
    await grok.dapi.userDataStorage.postValue(STORAGE_NAME, hash, record, false);
    this.currentUserPatternNameToHash.set(patternName, hash);
    this.eventBus.updatePatternList();
  }

  private validatePatternNameUniqueness(patternName: string): void {
    if (this.currentUserPatternNameToHash.has(patternName))
      throw new PatternNameExistsError(`Pattern with name ${patternName} already exists`);
  }

  private validatePatternUniqueness(hash: string): void {
    const existingHashes = Array.from(this.currentUserPatternNameToHash.values())
      .concat(Array.from(this.otherUsersPatternNameToHash.values()));

    if (existingHashes.includes(hash))
      throw new PatternExistsError(`Pattern with hash ${hash} already exists`);
  }

  getPatternHash(patternName: string, isCurrentUserPattern: boolean): string {
    const patternHash = isCurrentUserPattern ?
      this.currentUserPatternNameToHash.get(patternName) :
      this.otherUsersPatternNameToHash.get(patternName);

    console.warn('pattern data', patternName, isCurrentUserPattern, patternHash);

    if (patternHash === undefined) {
      console.warn('pattern data', patternName, isCurrentUserPattern, patternHash);
      throw new Error(`Pattern with name ${patternName} not found`);
    }

    return patternHash;
  }

  async getPatternConfig(hash: string): Promise<PatternConfiguration> {
    const patternConfig = await grok.dapi.userDataStorage.getValue(STORAGE_NAME, hash, false);
    const config = JSON.parse(patternConfig)[R.PATTERN_CONFIG] as PatternConfiguration;
    return config;
  }

  private getAuthorCategoryByHash(hash: string): string {
    if (this.isCurrentUserPattern(hash))
      return this.currentUserName;
    else if (this.isOtherUserPattern(hash))
      return OTHER_USERS;
    else
      throw new Error(`Pattern with hash ${hash} not found`);
  }

  private isCurrentUserPattern(hash: string): boolean {
    return Array.from(this.currentUserPatternNameToHash.values()).includes(hash);
  }

  private isOtherUserPattern(hash: string): boolean {
    return Array.from(this.otherUsersPatternNameToHash.values()).includes(hash);
  }

  getPatternNameByHash(hash: string): string {
    const maps = [this.currentUserPatternNameToHash, this.otherUsersPatternNameToHash];
    for (const map of maps) {
      for (const [patternName, patternHash] of map.entries()) {
        if (patternHash === hash)
          return patternName;
      }
    }
    throw new Error(`Pattern with hash ${hash} not found`);
  }

  private async deletePattern(patternName: string): Promise<void> {
    const hash = this.currentUserPatternNameToHash.get(patternName);
    if (hash === undefined)
      throw new Error(`Pattern with name ${patternName} not found`);
    await grok.dapi.userDataStorage.remove(STORAGE_NAME, hash, false);
    this.currentUserPatternNameToHash.delete(patternName);
    this.eventBus.updatePatternList();

    const patternToBeLoaded = this.getCurrentUserPatternNames()[0];
    const hashOfPatternToBeLoaded = this.currentUserPatternNameToHash.get(patternToBeLoaded);

    if (hashOfPatternToBeLoaded && patternToBeLoaded)
      this.eventBus.requestPatternLoad(hashOfPatternToBeLoaded);

    // todo: load default pattern in case there are no patterns left
  }

  private async getRecordFromPattern(patternConfig: PatternConfiguration): Promise<string> {
    const record = {
      [R.PATTERN_CONFIG]: patternConfig,
      [R.AUTHOR_ID]: await grok.dapi.users.current().then((u) => u.id),
    };
    const stringifiedRecord = JSON.stringify(record);
    return stringifiedRecord;
  }

  private getHash(patternConfig: PatternConfiguration): string {
    // WARNING: hash is computed over the graph settings only, as defining the
    // pattern's identity
    const graphSettings = GKL.reduce((acc, key) => {
      acc[key] = patternConfig[key as keyof PatternConfiguration];
      return acc;
    }, {} as any);
    return objectHash.sha1(graphSettings);
  }
}
