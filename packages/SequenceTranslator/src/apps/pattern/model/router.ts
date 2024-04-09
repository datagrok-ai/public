import {EventBus} from '../model/event-bus';

const PATTERN_KEY = 'pattern';

export class URLRouter {
  private urlSearchParams: URLSearchParams;
  constructor() {
    this.urlSearchParams = new URLSearchParams(window.location.search);
  }

  updateURLOnPatternLoaded(eventBus: EventBus): void {
    eventBus.patternLoaded$.subscribe((hash) => this.setPatternURL(hash));
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
