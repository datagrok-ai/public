import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {DataInitializer} from './data-initializer';

// export namespace URLRouter {
//   export function navigate(
//     searchParams: URLSearchParams,
//     dataInitializer: DataInitializer
//   ): void {
//     if (!searchParams.has('pattern'))
//       return;

//     const patternHash = searchParams.get('pattern');
//     if (!patternHash)
//       return;

//     try {
//       const patternConfig = await patternAppDataManager.getPatternConfig(patternHash);
//       eventBus.setPatternConfig(patternConfig);
//     } catch (e) {
//       searchParams.delete('pattern');
//       window.history.pushState({}, '', window.location.pathname);
//       console.warn('Error while getting pattern config by hash', e);
//     }
//   }
// }

export class URLRouter {
  private searchParams: URLSearchParams;

  constructor(private eventBus: EventBus, private patternAppDataManager: PatternAppDataManager) {
    this.searchParams = new URLSearchParams(window.location.search);

    this.eventBus.patternLoaded$.subscribe((hash) => this.setURL(hash));
  }

  async navigate(): Promise<void> {
    if (!this.searchParams.has('pattern'))
      return;

    const patternHash = this.searchParams.get('pattern');
    if (!patternHash)
      return;

    try {
      const patternConfig = await this.patternAppDataManager.getPatternConfig(patternHash);
      this.eventBus.setPatternConfig(patternConfig);
    } catch (e) {
      this.searchParams.delete('pattern');
      window.history.pushState({}, '', window.location.pathname);
      console.warn('Error while getting pattern config by hash', e);
    }
  }

  private setURL(patternHash: string): void {
    this.searchParams.set('pattern', patternHash);
    window.history.pushState({}, '', `${window.location.pathname}?${this.searchParams}`);
  }
}
