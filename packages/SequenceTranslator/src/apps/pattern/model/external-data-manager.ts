import {PatternConfiguration, PatternExistsError, PatternNameExistsError} from './types';
import {EventBus} from './event-bus';
import {PatternConfigManager} from './config-manager';

export class PatternAppDataManager {
  constructor(
    private eventBus: EventBus,
    private patternConfigManager: PatternConfigManager
  ) {
    this.eventBus.patternDeletionRequested$.subscribe(async (patternName: string) => {
      await this.deletePattern(patternName);
    });
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
      this.eventBus.updatePatternList();
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

  async deletePattern(patternName: string): Promise<void> {
    await this.patternConfigManager.deletePattern(patternName);
    this.eventBus.updatePatternList();

    const hashOfPatternToBeLoaded = this.patternConfigManager.getHashOfPatternToBeLoadedAfterDeletion();
    this.eventBus.requestPatternLoad(hashOfPatternToBeLoaded);
  }
}

