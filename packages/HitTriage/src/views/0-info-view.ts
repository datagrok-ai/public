import {HitTriageBaseView} from "./hit-triage-base-view";
import {_package} from "../package";
import * as ui from "datagrok-api/ui";
import {HitTriageApp} from "../hit-triage-app";

export class InfoView extends HitTriageBaseView {

  constructor(app: HitTriageApp) {
    super(app);
    _package.files.readAsText('README.md').then((md) => {
      this.root.appendChild(ui.markdown(md));
    });
  }
}