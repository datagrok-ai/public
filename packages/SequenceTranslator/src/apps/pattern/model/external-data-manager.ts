/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {isPatternCreatedByCurrentUser} from './utils';
import {STORAGE_NAME, GRAPH_SETTINGS_KEY_LIST as GKL} from './const';
import {PatternGraphSettings, PatternConfigRecord, PatternConfiguration, PatternLegendSettings} from './types';
import {EventBus} from './event-bus';
import {PatternConfigurationManager} from './pattern-config-manager';

import * as getHash from 'object-hash';

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

  async savePatternToUserStorage(patternName: string): Promise<void> {
    const patternConfig = PatternConfigurationManager.getConfig(this.eventBus);
    const userName = this.patternConfigManager.getCurrentUserName();
    await PatternConfigLoader.savePattern(patternConfig, patternName, userName);
  }

}

class PatternConfigManager {

  private otherUsersPatterns = new Map<string, string>();
  private currentUserPatterns = new Map<string, string>();

  constructor(private eventBus: EventBus) { }

  async init(): Promise<void> {
    try {
      const patternRecords = await PatternConfigLoader.fetchPatterns();

      this.eventBus.updatePatternList();
    } catch (e) {
      console.error('Error while initializing pattern list manager', e);
    }
  }

  getCurrentUserPatternNames(): string[] {
    return Array.from(this.currentUserPatterns.keys());
  }

  getOtherUsersPatternNames(): string[] {
    return Array.from(this.otherUsersPatterns.keys());
  }

  private async fetchCurrentUserName(): Promise<string> {
    const friendlyName = (await grok.dapi.users.current()).friendlyName;
    return friendlyName;
  }
}

namespace PatternConfigLoader {
  export async function fetchPatterns() {
    const patternsRecord = await grok.dapi.userDataStorage.get(STORAGE_NAME, false) as RawPatternRecords;
    return patternsRecord;
  }

  export async function savePattern(patternConfig: PatternConfiguration, patternName: string) {
    // create graphSettings from patternConfig
    const graphSettings = GKL.reduce((acc, key) => {
      acc[key] = patternConfig[key];
      return acc;
    }, {} as PatternGraphSettings);
    // await grok.dapi.userDataStorage.postValue(STORAGE_NAME, hash, stringifiedRecord, false);
  }

  function createPatternNamePostfix(friendlyUserName: string): string {
    return ` (created by ${friendlyUserName})`;
  }
}
