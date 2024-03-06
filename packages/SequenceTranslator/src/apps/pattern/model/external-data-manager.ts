/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {isPatternCreatedByCurrentUser} from './utils';
import {USER_STORAGE_KEY} from './const';
import {PatternConfiguration} from './types';
import {EventBus} from './event-bus';
import {PatternConfigurationManager} from './pattern-config-manager';

type PatternsRecord = {[patternName: string]: string};

export class PatternAppDataManager {
  constructor(private eventBus: EventBus) {
    this.patternListManager = new PatternListManager(eventBus);
    this.patternListManager.init();
  }
  private patternListManager: PatternListManager;

  getCurrentUserPatternNames(): string[] {
    return this.patternListManager.getCurrentUserPatternNames() || [];
  }

  getOtherUsersPatternNames(): string[] {
    return this.patternListManager.getOtherUsersPatternNames() || [];
  }

  getCurrentUserName(): string {
    return this.patternListManager.getCurrentUserName();
  }

  async savePatternToUserStorage(patternName: string): Promise<void> {
    const patternConfig = PatternConfigurationManager.getConfig(this.eventBus);
    const userName = this.patternListManager.getCurrentUserName();
    await PatternConfigLoader.savePatternToUserDataStorage(patternConfig, patternName, userName);
  }

}

class PatternListManager {
  private currentUserPatterns: string[] = [];
  private otherUsersPatterns: string[] = [];
  private currentUserFriendlyName = '';

  private eventBus: EventBus

  constructor(eventBus: EventBus) {
    this.eventBus = eventBus;
  }

  async init(): Promise<void> {
    try {
      const patternRecords = await PatternConfigLoader.fetchPatternsFromUserDataStorage();
      const categorizedPatterns = await this.categorizePatternsByUserOwnership(patternRecords);

      this.currentUserPatterns = categorizedPatterns.currentUserPatterns;
      this.otherUsersPatterns = categorizedPatterns.otherUsersPatterns;
      this.currentUserFriendlyName = await this.fetchCurrentUserName();

      this.eventBus.updatePatternList();
    } catch (e) {
      console.error('Error while initializing pattern list manager', e);
    }
  }

  getCurrentUserPatternNames(): string[] {
    return this.currentUserPatterns;
  }

  getOtherUsersPatternNames(): string[] {
    return this.otherUsersPatterns;
  }

  getCurrentUserName(): string {
    return this.currentUserFriendlyName;
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

  private async fetchCurrentUserName(): Promise<string> {
    const friendlyName = (await grok.dapi.users.current()).friendlyName;
    return friendlyName;
  }
}

namespace PatternConfigLoader {
  export async function fetchPatternsFromUserDataStorage() {
    const patternsData = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false) as PatternsRecord;
    return patternsData;
  }

  export async function savePatternToUserDataStorage(patternConfig: PatternConfiguration, patternName: string, userName: string) {
    const stringifiedPatternConfig = JSON.stringify(patternConfig);
    const fullPatternName = patternName + createPatternNamePostfix(userName);
    await grok.dapi.userDataStorage.postValue(USER_STORAGE_KEY, fullPatternName, stringifiedPatternConfig, false);
  }

  function createPatternNamePostfix(userName: string): string {
    return ` (created by ${userName})`;
  }
}
