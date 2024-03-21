/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  STORAGE_NAME, GRAPH_SETTINGS_KEY_LIST as GKL, LEGEND_SETTINGS_KEYS as L, PATTERN_RECORD_KEYS as R, OTHER_USERS
} from './const';
import {PatternConfiguration, PatternExistsError, PatternNameExistsError} from './types';
import {EventBus} from './event-bus';

import objectHash from 'object-hash';

type RawPatternRecords = {[patternName: string]: string};

export class PatternAppDataManager {
  // WARNING: not a singleton, encapsulated async init logic in a static method
  private constructor(
    private eventBus: EventBus,
    private patternConfigManager: PatternConfigManager
  ) { }

  static async getInstance(eventBus: EventBus): Promise<PatternAppDataManager> {
    const patternConfigManager = await PatternConfigManager.getInstance(eventBus);
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

  getPatternConfig(patternName: string, isCurrentUserPattern: boolean): Promise<PatternConfiguration> {
    return this.patternConfigManager.getPatternConfig(patternName, isCurrentUserPattern);
  }
}

class PatternConfigManager {
  // todo: remove
  private currentUserName: string;
  private currentUserId: string;

  private otherUsersPatternNameToHash = new Map<string, string>();
  private currentUserPatternNameToHash = new Map<string, string>();

  private constructor(private eventBus: EventBus) { }

  private async init(): Promise<void> {
    try {
      this.currentUserName = await this.fetchCurrentUserName();
      this.currentUserId = await this.fetchCurrentUserId();
      const patternRecords = await PatternConfigLoader.fetchPatterns();
      await this.initializePatterns(patternRecords);

      this.eventBus.updatePatternList();
    } catch (e) {
      console.error('Error while initializing pattern list manager', e);
    }
  }

  static async getInstance(eventBus: EventBus): Promise<PatternConfigManager> {
    const instance = new PatternConfigManager(eventBus);
    await instance.init();

    eventBus.patternDeletionRequested$.subscribe(async (patternName: string) => {
      await instance.deletePattern(patternName);
    });

    return instance;
  }

  private async initializePatterns(patternRecords: RawPatternRecords) {
    if (!this.currentUserId)
      throw new Error('Current user ID is not set');

    const userIdsToUserNames = new Map<string, string>();

    for (const [patternHash, stringifiedRecord] of Object.entries(patternRecords))
      await this.extractDataFromRecordToMaps(patternHash, stringifiedRecord, userIdsToUserNames);
  }

  private async extractDataFromRecordToMaps(
    patternHash: string,
    stringifiedRecord: string,
    userIdsToUserNames: Map<string, string>
  ) {
    const record = JSON.parse(stringifiedRecord);
    const patternConfig = record[R.PATTERN_CONFIG] as PatternConfiguration;
    const patternName = patternConfig.patternName;
    const authorID = record[R.AUTHOR_ID];
    if (authorID === this.currentUserId) {
      this.currentUserPatternNameToHash.set(patternName, patternHash);
    } else {
      if (!userIdsToUserNames.has(authorID)) {
        const userFriendlyName = (await grok.dapi.users.find(authorID)).friendlyName;
        userIdsToUserNames.set(authorID, userFriendlyName);
      }
      const fullPatternName = patternName + ` (created by ${userIdsToUserNames.get(authorID)})`;
      this.otherUsersPatternNameToHash.set(fullPatternName, patternHash);
    }
  }

  private async fetchCurrentUserId(): Promise<string> {
    return (await grok.dapi.users.current()).id;
  }

  getCurrentUserPatternNames(): string[] {
    return Array.from(this.currentUserPatternNameToHash.keys())
      .sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
  }

  getOtherUsersPatternNames(): string[] {
    return Array.from(this.otherUsersPatternNameToHash.keys())
      .sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
  }

  private async fetchCurrentUserName(): Promise<string> {
    const friendlyName = (await grok.dapi.users.current()).friendlyName;
    return friendlyName;
  }

  getCurrentUserName(): string {
    return this.currentUserName;
  }

  async savePatternToUserStorage(patternConfig: PatternConfiguration): Promise<void> {
    const hash = PatternConfigLoader.getHash(patternConfig);
    this.validatePatternUniqueness(hash);

    const patternName = patternConfig[L.PATTERN_NAME];
    this.validatePatternNameUniqueness(patternName);

    const record = await PatternConfigLoader.getRecordFromPattern(patternConfig);
    await grok.dapi.userDataStorage.postValue(STORAGE_NAME, hash, record, false);
    this.currentUserPatternNameToHash.set(patternName, hash);
    this.eventBus.updatePatternList();
  }

  private validatePatternNameUniqueness(patternName: string): void {
    if (this.currentUserPatternNameToHash.has(patternName))
      throw new PatternNameExistsError(`Pattern with name ${patternName} already exists`);
  }

  private validatePatternUniqueness(hash: string): void {
    const existingHashes = Array.from(this.currentUserPatternNameToHash.values());

    if (existingHashes.includes(hash))
      throw new PatternExistsError(`Pattern with hash ${hash} already exists`);
  }

  async getPatternConfig(patternName: string, isCurrentUserPattern: boolean): Promise<PatternConfiguration> {
    const patternHash = isCurrentUserPattern ?
      this.currentUserPatternNameToHash.get(patternName) :
      this.otherUsersPatternNameToHash.get(patternName);

    if (patternHash === undefined)
      throw new Error(`Pattern with name ${patternName} not found`);

    return await this.getPatternConfigByHash(patternHash);
  }

  private async getPatternConfigByHash(hash: string): Promise<PatternConfiguration> {
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

  private getPatternNameByHash(hash: string): string {
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
    if (patternToBeLoaded)
      this.eventBus.requestPatternLoad(patternToBeLoaded);

    // todo: load default pattern in case there are no patterns left
  }
}

namespace PatternConfigLoader {
  export async function fetchPatterns() {
    const patternsRecord = await grok.dapi.userDataStorage.get(STORAGE_NAME, false) as RawPatternRecords;
    return patternsRecord;
  }

  export async function getRecordFromPattern(patternConfig: PatternConfiguration): Promise<string> {
    const record = {
      [R.PATTERN_CONFIG]: patternConfig,
      [R.AUTHOR_ID]: await grok.dapi.users.current().then((u) => u.id),
    };
    const stringifiedRecord = JSON.stringify(record);
    return stringifiedRecord;
  }

  export function getHash(patternConfig: PatternConfiguration): string {
    // WARNING: hash is computed over the graph settings only, as defining the
    // pattern's identity
    const graphSettings = GKL.reduce((acc, key) => {
      acc[key] = patternConfig[key as keyof PatternConfiguration];
      return acc;
    }, {} as any);
    return objectHash.sha1(graphSettings);
  }
}
