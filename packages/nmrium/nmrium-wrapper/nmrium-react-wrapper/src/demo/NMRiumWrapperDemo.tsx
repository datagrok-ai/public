// /** @jsxImportSource @emotion/react */
// import { css } from '@emotion/react';
// import { NMRiumData } from 'nmrium';
// import Button from 'nmrium/lib/component/elements/Button';

// import NMRiumWrapper from '../NMRiumWrapper';
// import events from '../events';
// import { loadFilesFromURLs } from '../utilities/loadFilesFromURLs';

// import jsonData from './data/test.json';

// const styles = {
//   container: css`
//     display: flex;
//     flex-direction: column;
//     height: 100%;
//   `,
//   header: css`
//     height: 40px;
//     width: 100%;
//     padding: 5px;
//     display: flex;
//   `,
// };

// export default function NMRiumWrapperDemo() {
//   return (
//     <div css={styles.container}>
//       <div id="header" css={styles.header}>
//         <Button.Done
//           style={{ marginRight: '10px' }}
//           onClick={() => {
//             events.trigger('load', {
//               data: jsonData as NMRiumData,
//               type: 'nmrium',
//             });
//           }}
//         >
//           Test load from json
//         </Button.Done>

//         <Button.Done
//           style={{ marginRight: '10px' }}
//           onClick={() => {
//             events.trigger('load', {
//               data: [
//                 'https://cheminfo.github.io/nmr-dataset-demo/cytisine/13c.jdx',
//                 'https://cheminfo.github.io/nmr-dataset-demo/cytisine/1h.jdx',
//                 'https://cheminfo.github.io/bruker-data-test/data/zipped/aspirin-1h.zip',
//                 // '../data/13c.zip',
//                 // 'https://cloud.uni-jena.de/s/y72GbCX8bJbmpJT/download/10.zip',
//                 // 'https://cloud.uni-jena.de/s/jsMed9fmqWZzo6r/download/53.zip',
//               ],
//               type: 'url',
//             });
//           }}
//         >
//           Test Load from URLS
//         </Button.Done>

//         <Button.Done
//           style={{ marginRight: '10px' }}
//           onClick={() => {
//             events.trigger('load', {
//               data: [
//                 'https://cheminfo.github.io/bruker-data-test/data/zipped/aspirin-1h',
//               ],
//               type: 'url',
//             });
//           }}
//         >
//           Test Load URL without extension
//         </Button.Done>
//         <Button.Done
//           style={{ marginRight: '10px' }}
//           onClick={() => {
//             void loadFilesFromURLs([
//               '../data/COSY-12.zip',
//               '../data/HMBC-13.zip',
//               '../data/13c.zip',
//             ]).then((files) => {
//               events.trigger('load', {
//                 data: files,
//                 type: 'file',
//                 activeTab: '13C',
//               });
//             });
//           }}
//         >
//           Test Load Files
//         </Button.Done>
//         <Button.Done
//           className="logger-btn"
//           onClick={() => {
//             void loadFilesFromURLs(['../data/sample-with-error.zip']).then(
//               (files) => {
//                 events.trigger('load', {
//                   data: files,
//                   type: 'file',
//                   activeTab: '13C',
//                 });
//               },
//             );
//           }}
//         >
//           Test Logger
//         </Button.Done>
//       </div>

//       <NMRiumWrapper />
//     </div>
//   );
// }
