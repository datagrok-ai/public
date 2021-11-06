/** EXPERIMENTAL - API might change */
import * as ui from "../../ui";
import {Dialog} from "../widgets";
import $ from "cash-dom";

interface WizardPage {

  root: HTMLDivElement;

  /** Caption to be displayed on top of the panel */
  caption?: HTMLElement | string;

  /** Called when the page is activated */
  onActivated?: () => void;

  /** Returns error message (and stops wizard from proceeding to the next page),
   * or null if validated */
  validate?: () => string;
}

/** A set of pages that user navigates using the "<<" and ">>" buttons */
export class Wizard {
  captionHost = ui.div([], 'ui-wizard-page-title');
  contentHost = ui.div([], 'ui-wizard-page-content');
  root: HTMLDivElement = ui.divV([this.captionHost, this.contentHost]);
  pages: WizardPage[] = [];
  _currentPage: WizardPage = { root: ui.div() };
  dialog: Dialog = ui.dialog().add(this.root);
  okButton = this.dialog.getButton('OK');
  prevButton = this.dialog.addButton('<<', () => this.prev());
  nextButton = this.dialog.addButton('>>', () => this.next());

  Wizard(options?: {title: string}) {
    this.dialog.title = options?.title ?? '';
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
    ui.setClass(this.dialog.getButton('OK'), 'disabled', !this.completable);
    ui.setClass(this.dialog.getButton('<<'), 'disabled', this.pageIndex == 0);
    ui.setClass(this.dialog.getButton('>>'), 'disabled', this.pageIndex == this.pages.length - 1);
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