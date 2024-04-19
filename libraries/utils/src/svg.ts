import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {testEvent} from './test';
import {Subject} from 'rxjs';

export function svgToImage(svgEl: SVGSVGElement, dpr: number, timeout: number = 1000): Promise<HTMLImageElement> {
  const w: number = svgEl.width.baseVal.value;
  const h: number = svgEl.height.baseVal.value;
  const svgImg = new Image(dpr * w, dpr * h);
  const svgXml = new XMLSerializer().serializeToString(svgEl);

  const svg64Prefix: string = 'data:image/svg+xml;base64,';
  const svg64: string = window.btoa(svgXml);

  const svgImgOnLoad = new Subject<void>();
  const svgImgOnLoadListener = () => {
    svgImgOnLoad.next();
  };
  svgImg.addEventListener('load', svgImgOnLoadListener);

  return new Promise<HTMLImageElement>((resolve, reject) => {
    testEvent(svgImgOnLoad,
      () => { resolve(svgImg); },
      () => { svgImg.src = svg64Prefix + svg64; },
      timeout)
      .catch((err) => { reject(err); })
      .finally(() => {
        svgImg.removeEventListener('load', svgImgOnLoadListener);
      });
  });
}
