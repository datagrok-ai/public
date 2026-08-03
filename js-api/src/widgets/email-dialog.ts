/**
 * EmailDialog — generic email dialog with markdown body and attachments.
 * @module widgets/email-dialog
 */

import {toDart, toJs} from '../wrappers';
import {IDartApi} from '../api/grok_api.g';
import {DartWrapper} from './base';
import {FileInfo} from '../entities';

const api: IDartApi = (typeof window !== 'undefined' ? window : (global as any).window) as any;


/**
 * Final values produced by [EmailDialog] when SEND is pressed.
 * Map directly into `grok.dapi.admin.sendEmail` (`html`/`text`/`attachments`).
 */
export interface EmailDialogResult {
  subject: string;
  to: string[];
  bcc: string[];
  /** Rendered HTML body. Base64 `<img>` tags are hoisted to `cid:` inline
   * parts by the server (Mailgun) — the dialog ships the body as-is. */
  html: string;
  /** Plain-text fallback derived from [html]. */
  text: string;
  attachments: FileInfo[];
}


export interface EmailDialogOptions {
  /** Dialog title. */
  title?: string;
  subject?: string;
  to?: string[];
  bcc?: string[];
  /** Initial markdown body. */
  body?: string;
  attachments?: FileInfo[];
  toEditable?: boolean;
  bccEditable?: boolean;
  helpUrl?: string;
  /** Max total attachment size in bytes. Defaults to 20 MB (server limit). */
  maxAttachmentBytes?: number;
  /**
   * Handler invoked with the built result when SEND is pressed. The dialog
   * disables the SEND button, awaits the handler, closes on success, and
   * shows a balloon + keeps itself open on error.
   */
  onSend?: (result: EmailDialogResult) => void | Promise<void>;
}


/**
 * Generic markdown-body email dialog with attachments.
 */
export class EmailDialog extends DartWrapper {

  static create(opts?: EmailDialogOptions): EmailDialog {
    const d = new EmailDialog(api.grok_EmailDialog());
    if (opts != null) {
      if (opts.title != null) d.title = opts.title;
      if (opts.subject != null) d.subject = opts.subject;
      if (opts.to != null) d.to = opts.to;
      if (opts.bcc != null) d.bcc = opts.bcc;
      if (opts.body != null) d.body = opts.body;
      if (opts.attachments != null) d.attachments = opts.attachments;
      if (opts.toEditable != null) d.toEditable = opts.toEditable;
      if (opts.bccEditable != null) d.bccEditable = opts.bccEditable;
      if (opts.helpUrl != null) d.helpUrl = opts.helpUrl;
      if (opts.maxAttachmentBytes != null) d.maxAttachmentBytes = opts.maxAttachmentBytes;
      if (opts.onSend != null) d.onSend(opts.onSend);
    }
    return d;
  }

  set title(v: string) { api.grok_EmailDialog_Set_Title(this.dart, v); }
  set subject(v: string) { api.grok_EmailDialog_Set_Subject(this.dart, v); }
  set to(v: string[]) { api.grok_EmailDialog_Set_To(this.dart, v); }
  set bcc(v: string[]) { api.grok_EmailDialog_Set_Bcc(this.dart, v); }
  set body(v: string) { api.grok_EmailDialog_Set_Body(this.dart, v); }
  set attachments(v: FileInfo[]) { api.grok_EmailDialog_Set_Attachments(this.dart, toDart(v)); }
  set toEditable(v: boolean) { api.grok_EmailDialog_Set_ToEditable(this.dart, v); }
  set bccEditable(v: boolean) { api.grok_EmailDialog_Set_BccEditable(this.dart, v); }
  set helpUrl(v: string) { api.grok_EmailDialog_Set_HelpUrl(this.dart, v); }
  set maxAttachmentBytes(v: number) { api.grok_EmailDialog_Set_MaxAttachmentBytes(this.dart, v); }

  /** Shows the dialog. Optionally centers near [centerAt]. */
  show(opts?: {centerAt?: HTMLElement}): this {
    api.grok_EmailDialog_Show(this.dart, opts?.centerAt ?? null);
    return this;
  }

  close(): void { api.grok_EmailDialog_Close(this.dart); }

  onSend(handler: (result: EmailDialogResult) => void | Promise<void>): this {
    api.grok_EmailDialog_SetOnSend(this.dart, (resultDart: any) => handler(_readResult(resultDart)));
    return this;
  }
}


function _readResult(dart: any): EmailDialogResult {
  return {
    subject: api.grok_EmailDialogResult_Get_Subject(dart),
    to: toJs(api.grok_EmailDialogResult_Get_To(dart)) ?? [],
    bcc: toJs(api.grok_EmailDialogResult_Get_Bcc(dart)) ?? [],
    html: api.grok_EmailDialogResult_Get_Html(dart),
    text: api.grok_EmailDialogResult_Get_Text(dart),
    attachments: toJs(api.grok_EmailDialogResult_Get_Attachments(dart)) ?? [],
  };
}
