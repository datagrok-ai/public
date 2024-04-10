import {EventBus} from '../model/event-bus';

const PATTERN_KEY = 'pattern';

export class URLRouter {
  private urlSearchParams: URLSearchParams;
  constructor() {
    this.urlSearchParams = new URLSearchParams(window.location.search);
  }

  subscribeToObservables(eventBus: EventBus): void {
    eventBus.urlStateUpdated$.subscribe((hash) => this.setPatternURL(hash));
    eventBus.loadPatternInNewTabRequested$.subscribe((hash) => {
      const url = `${window.location.origin}${window.location.pathname}?${PATTERN_KEY}=${hash}`;
      window.open(url, '_blank');
    });

    window.addEventListener('popstate', () => {
      this.urlSearchParams = new URLSearchParams(window.location.search);
      const patternHash = this.getPatternHash();
      if (patternHash === null)
        return;

      eventBus.requestPatternLoad(patternHash);
    });
  }

  getPatternHash(): string | null {
    return this.urlSearchParams.get(PATTERN_KEY);
  }

  setPatternURL(patternHash: string): void {
    if (patternHash === null || patternHash === '') {
      this.clearPatternURL();
      return;
    }

    this.urlSearchParams.set(PATTERN_KEY, patternHash);
    window.history.pushState({}, '', `${window.location.pathname}?${this.urlSearchParams}`);
  }

  clearPatternURL(): void {
    this.urlSearchParams.delete(PATTERN_KEY);
    window.history.pushState({}, '', `${window.location.pathname}`);
  }
}
