import {RichFunctionView} from '../../../function-views';


export class RFVPopup extends RichFunctionView {
  override detach() { }

  customClose() {
    super.detach();
    super.close();
  }
}
