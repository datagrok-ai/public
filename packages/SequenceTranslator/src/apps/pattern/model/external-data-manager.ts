/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {STORAGE_NAME, GRAPH_SETTINGS_KEY_LIST as GKL, PATTERN_RECORD_KEYS as R} from './const';
import {PatternGraphSettings, PatternConfigRecord, PatternConfiguration, PatternLegendSettings} from './types';
import {EventBus} from './event-bus';
import {PatternConfigurationManager} from './pattern-config-manager';

import objectHash from 'object-hash';

type RawPatternRecords = {[patternName: string]: string};

export class PatternAppDataManager {
  constructor(private eventBus: EventBus) {
    this.patternConfigManager = new PatternConfigManager(eventBus);
    this.patternConfigManager.init();
  }
  private patternConfigManager: PatternConfigManager;

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
    const patternConfig = PatternConfigurationManager.getConfig(this.eventBus);
    await PatternConfigLoader.savePattern(patternConfig);
  }

}

class PatternConfigManager {
  // todo: remove
  private currentUserName: string;
  private currentUserId: string;

  private otherUsersPatternNameToHash = new Map<string, string>();
  private currentUserPatternNameToHash = new Map<string, string>();

  constructor(private eventBus: EventBus) { }

  async init(): Promise<void> {
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

  private async initializePatterns(patternRecords: RawPatternRecords) {
    if (!this.currentUserId) {
      throw new Error('Current user ID is not set');
    }

    const userIdsToUserNames = new Map<string, string>();

    for (const [patternHash, stringifiedRecord] of Object.entries(patternRecords)) {
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

  private async fetchCurrentUserId(): Promise<string> {
    return (await grok.dapi.users.current()).id;
  }

  getCurrentUserPatternNames(): string[] {
    return Array.from(this.currentUserPatternNameToHash.keys());
  }

  getOtherUsersPatternNames(): string[] {
    return Array.from(this.otherUsersPatternNameToHash.keys());
  }

  private async fetchCurrentUserName(): Promise<string> {
    const friendlyName = (await grok.dapi.users.current()).friendlyName;
    return friendlyName;
  }

  getCurrentUserName(): string {
    return this.currentUserName;
  }
}

namespace PatternConfigLoader {
  export async function fetchPatterns() {
    const patternsRecord = await grok.dapi.userDataStorage.get(STORAGE_NAME, false) as RawPatternRecords;
    console.log(`patterns record:`, patternsRecord);
    return patternsRecord;
  }

  export async function savePattern(patternConfig: PatternConfiguration) {
    const hash = getHash(patternConfig);
    const record = {
      [R.PATTERN_CONFIG]: patternConfig,
      [R.AUTHOR_ID]: await grok.dapi.users.current().then(u => u.id),
    };
    const stringifiedRecord = JSON.stringify(record);
    // console.log(`hash:`, hash);
    // console.log(`record:`, record);
    await grok.dapi.userDataStorage.postValue(STORAGE_NAME, hash, stringifiedRecord, false);
  }

  function getHash(patternConfig: PatternConfiguration): string {
    // WARNING: hash is computed over the graph settings only, as they define pattern identity
    const graphSettings = GKL.reduce((acc, key) => {
      acc[key] = patternConfig[key as keyof PatternConfiguration];
      return acc;
    }, {} as any);
    return objectHash.sha1(graphSettings);
  }
}
