/** EXPERIMENTAL - USE AT YOUR OWN RISK - API might change */
import * as ui from "../../ui";
import {Dialog} from "../widgets";
let api = <any>window;

/** EXPERIMENTAL - USE AT YOUR OWN RISK - API might change */
interface WizardPage {

  root: HTMLDivElement;

  /** Caption to be displayed on top of the panel */
  caption?: HTMLElement | string;

  /** Called when the page is activated */
  onActivated?: () => void;

  /** Returns error message (and stops wizard from proceeding to the next page),
   * or null if validated */
  validate?: () => string | null;
}

/** A set of pages that user navigates using the "<<" and ">>" buttons */
/** EXPERIMENTAL - USE AT YOUR OWN RISK - API might change */
export class Wizard extends Dialog {
  captionHost = ui.div([], 'ui-wizard-page-title');
  contentHost = ui.div([], 'ui-wizard-page-content');
  wizardRoot: HTMLDivElement = ui.divV([this.captionHost, this.contentHost]);
  pages: WizardPage[] = [];
  _currentPage: WizardPage = { root: ui.div() };
  okButton = this.getButton('OK');
  prevButton = this.addButton('<<', () => this.prev());
  nextButton = this.addButton('>>', () => this.next());

  constructor(options?: {title?: string, helpUrl?: string}) {
    super(api.grok_Dialog(options?.title, options?.helpUrl, true, true));
    this.add(this.wizardRoot);
  }

  page(p: WizardPage): Wizard {
    this.pages.push(p);
    if (this.pageIndex == -1)
      this.currentPage = p;
    this._updateButtonStates();
    return this;
  }

  get pageIndex(): number { return this.pages.indexOf(this._currentPage); }

  get completable(): boolean {
    return this.pageIndex == this.pages.length - 1 &&
      (this._currentPage.validate == undefined || this._currentPage.validate() == null);
  }

  get currentPage(): WizardPage { return this._currentPage; }
  set currentPage(w) {
    ui.empty(this.captionHost);
    ui.empty(this.contentHost);
    this.contentHost.appendChild(w.root);
    this.captionHost.appendChild(ui.render(w.caption));
    this._currentPage = w;
    if (w.onActivated)
      w.onActivated();
    this._updateButtonStates();
  }

  _updateButtonStates() {
    ui.setClass(this.getButton('OK'), 'disabled', !this.completable);
    ui.setClass(this.getButton('<<'), 'disabled', this.pageIndex == 0);
    ui.setClass(this.getButton('>>'), 'disabled', this.pageIndex == this.pages.length - 1);
  }

  /** Activates the previous page */
  prev() {
    this.currentPage = this.pages[this.pageIndex - 1];
    this._updateButtonStates();
  }

  /** Activates the previous page */
  next() {
    this.currentPage = this.pages[this.pageIndex + 1];
    this._updateButtonStates();
  }
}