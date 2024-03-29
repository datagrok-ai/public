import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';

export class Router {
  private searchParams: URLSearchParams;

  constructor(private eventBus: EventBus, private patternAppDataManager: PatternAppDataManager) {
    this.searchParams = new URLSearchParams(window.location.search);
  }

  async navigate(): Promise<void> {
    if (!this.searchParams.has('pattern'))
      return;

    const patternHash = this.searchParams.get('pattern');
    if (!patternHash)
      return;

    console.log(`pattern hash:`, patternHash);
    try {
      const patternConfig = await this.patternAppDataManager.getPatternConfigByHash(patternHash);
      this.eventBus.setPatternConfig(patternConfig);
      console.log(`pattern config:`, patternConfig);
    } catch (e) {
      this.searchParams.delete('pattern');
      window.history.pushState({}, '', window.location.pathname);
      console.warn('Error while getting pattern config by hash', e);
    }
  }
}
