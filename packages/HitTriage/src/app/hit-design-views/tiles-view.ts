import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {TileCategoriesColName} from '../consts';
import './utils.css';

export function getTilesViewDialog(app: HitDesignApp, getTableView: () => DG.TableView | null) {
  const tilesViewerSketchStateString = app.campaign?.tilesViewerFormSketch;
  let sketchState: any | null = null;
  if (tilesViewerSketchStateString && tilesViewerSketchStateString.length > 0) {
    try {
      sketchState = JSON.parse(tilesViewerSketchStateString);
    } catch (e) {
      console.error('Failed to parse sketch state string', e);
    }
  }
  if (sketchState) {
    // after saving, some columns might be removed in the meantime, so we need to update the
    //sketch state and remove them so that the viewer can be initialized correctly.
    if (sketchState.elementStates && (sketchState.elementStates.length ?? 0) > 0) {
      sketchState.elementStates = sketchState.elementStates.filter((elementState: any) =>
        elementState?.viewerSettings?.column && app.dataFrame!.col(elementState.viewerSettings.column),
      );
      // similarly, because the table name can be duplicated or changed,
      // we need to make sure that table name in all viewer settings corresponds to actual table name
      sketchState.elementStates.forEach((elState: any) => {
        if (elState?.viewerSettings?.table)
          elState.viewerSettings.table = app.dataFrame!.name;
      });
    }
  }

  const tileOpts = {lanesColumnName: TileCategoriesColName, lanes: app.template?.stages ?? [],
    ...((sketchState?.elementStates?.length ?? 0) > 0 ? {sketchState} : {})};

  const tv = getTableView();

  if (!tv) {
    grok.shell.error('Failed to create tiles viewer. table view is not available.');
    return;
  }
  let v: DG.Viewer | null = null;
  const modal = ui.dialog('Progress Tracker');
  try {
    v = tv.addViewer(DG.VIEWER.TILE_VIEWER, tileOpts);

    const opts = v.getOptions();
    opts.look.sketchState = (sketchState?.elementStates?.length ?? 0) > 0 ? sketchState : null;

    v.setOptions(opts);
    v.copyViewersLook(v); // hacky way to apply sketch state
  } catch (e) {
    grok.shell.error('Failed to create tiles viewer with preset layout.\n Falling back to default layout.');
    console.error('Failed to create tiles viewer', e);
    if (v) {
      v.close();
      v.detach();
    }
    v = tv.addViewer(DG.VIEWER.TILE_VIEWER,
      {lanesColumnName: TileCategoriesColName, lanes: app.template?.stages ?? []});
  }

  if (!v) {
    grok.shell.error('Failed to create tiles viewer. check the console for more details.');
    return;
  }
  modal.add(v.root);

  // from this modal, only way to go to another view is through edit form.
  // when this happens, we should not destroy the viewer,
  //but just hide it so that when we come back, we can show it again.
  const hideModal = () => modal.root.style.setProperty('display', 'none', 'important');
  const showModal = () => modal.root.style.display = 'flex';

  const viewChangeSub = grok.events.onCurrentViewChanged.subscribe(() => {
    if (grok.shell.v === tv)
      showModal();
    else
      hideModal();
  });


  const closeSub = modal.onClose.subscribe(() => {
    closeSub.unsubscribe();
    viewChangeSub.unsubscribe();
    // save the sketch state
    try {
      const sketchState = v.props.sketchState;
      if (sketchState && typeof sketchState == 'object') {
        const stringState = JSON.stringify(sketchState);
        if (app.campaign)
          app.campaign.tilesViewerFormSketch = stringState;
      }
    } catch (e) {
      grok.shell.error('Failed to save sketch state. check the console for more details.');
      console.error('Failed to save sketch state', e);
    }
    v.detach();
    v.close();
  });

  modal.showModal(true);


  // remove the modal background
  document.querySelector('.d4-modal-background')?.remove();
}
