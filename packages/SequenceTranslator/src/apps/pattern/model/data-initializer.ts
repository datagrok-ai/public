/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';

import {PATTERN_RECORD_KEYS as R, STORAGE_NAME} from './const';
import {PatternConfiguration, RawPatternRecords} from './types';
import {PatternDefaultsProvider} from './defaults-provider';

export class DataInitializer {
  private currentUserName: string;
  private currentUserId: string;
  private otherUsersPatternNameToHash = new Map<string, string>();
  private currentUserPatternNameToHash = new Map<string, string>();

  // WARNING: not a singleton, encapsulates async init logic in getInstance
  // method
  private constructor() { }

  static async getInstance(
    defaultsProvider: PatternDefaultsProvider
  ): Promise<DataInitializer> {
    const instance = new DataInitializer();

    instance.currentUserName = await instance.fetchCurrentUserName();
    instance.currentUserId = await instance.fetchCurrentUserId();

    const patternRecords = await instance.fetchPatterns();
    await instance.initializePatternMaps(patternRecords);

    return instance;
  }

  getCurrentUserName(): string {
    return this.currentUserName;
  }

  getCurrentUserId(): string {
    return this.currentUserId;
  }

  getOtherUsersPatternNameToHash(): Map<string, string> {
    return this.otherUsersPatternNameToHash;
  }

  getCurrentUserPatternNameToHash(): Map<string, string> {
    return this.currentUserPatternNameToHash;
  }

  private async fetchCurrentUserName(): Promise<string> {
    return (await grok.dapi.users.current()).friendlyName;
  }

  private async fetchCurrentUserId(): Promise<string> {
    return (await grok.dapi.users.current()).id;
  }

  private async fetchPatterns() {
    const patternsRecord = await grok.dapi.userDataStorage.get(STORAGE_NAME, false) as RawPatternRecords;
    return patternsRecord;
  }

  private async initializePatternMaps(patternRecords: RawPatternRecords) {
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
}
