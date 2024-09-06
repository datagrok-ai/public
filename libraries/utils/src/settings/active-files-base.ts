import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type FilesSettings = {
  excluded: string[],
  explicit: string[]
}

export class ActiveFiles {
  private settings?: FilesSettings;

  protected constructor(
    public readonly path: string,
    public readonly userStorageName: string,
    /**file extension to filter files in file open dialog*/public readonly ext: string = '.csv'
  ) {
    if ('/' !== this.path[this.path.length - 1]) this.path += '/';
    if ('.' !== this.ext[0]) this.ext = '.' + this.ext;
  }

  async getUserSettings(): Promise<FilesSettings> {
    if (this.settings)
      return this.settings;

    const resStr: string | undefined = grok.userSettings.getValue(this.userStorageName, 'Settings');
    const res = resStr ? JSON.parse(resStr) : {excluded: [], explicit: []};

    res.explicit = res.explicit instanceof Array ? res.explicit : [];
    res.excluded = res.excluded instanceof Array ? res.excluded : [];

    return this.settings = res;
  }

  setUserSettings(value: FilesSettings) {
    this.settings = value;
    grok.userSettings.add(this.userStorageName, 'Settings', JSON.stringify(value));
  }

  protected async getAllAvailable(): Promise<string[]> {
    const list = await grok.dapi.files.list(this.path);
    const paths = list.map((fileInfo) => {
      return fileInfo.fullPath.replace(`${this.path}`, '');
    });

    return paths;
  }

  public async getActive(): Promise<string[]> {
    const settings = await this.getUserSettings();
    if (settings.explicit.length > 0)
      return settings.explicit;

    return (await this.getAllAvailable()).filter((available) => !settings.excluded.includes(available));
  }

  protected createInput(available: string, isChecked: boolean): DG.InputBase<boolean> {
    const cb = ui.input.bool(available, {
      value: isChecked,
      onValueChanged: (inp, value) => this.updateSelectionStatus(available, value)
    }) as DG.InputBase<boolean>;
    cb.addOptions(ui.button(ui.iconFA('trash'), () => {
      const remDialog = ui.dialog({title: 'Warning'}).add(ui.divText(`Delete file '${available}'?`)).onOK(() => {
        cb.root.remove();
        this.availableRemove(available); //async
      });

      remDialog.show();
    }, `Delete ${available}`));

    return cb;
  }

  protected async getInputs(): Promise<DG.InputBase<boolean>[]> {
    const all = await this.getAllAvailable();
    const active = await this.getActive();
    const inputs: DG.InputBase<boolean>[] = new Array<DG.InputBase<boolean>>(all.length);

    for (let i = 0; i < all.length; i++) {
      const isChecked = active.includes(all[i]);
      inputs[i] = this.createInput(all[i], isChecked);
    }

    return inputs;
  }

  public async getForm(): Promise<HTMLDivElement> {
    const inputs = await this.getInputs();
    const inputsDiv = ui.divV(inputs);
    const formDiv = ui.divV( [
      inputsDiv,
      ui.button('ADD', async () => {
        let newAvailable: string = '';
        newAvailable = await this.getNewAvailable();
        if (newAvailable !== '')
          inputsDiv.append(this.createInput(newAvailable, true).root);
      })
    ]);

    return formDiv;
  }

  private async updateSelectionStatus(available: string, newisChecked: boolean): Promise<void> {
    const settings = await this.getUserSettings();
    const oldisChecked = !settings.excluded.includes(available);
    if (oldisChecked === newisChecked)
      return;

    if (newisChecked) {
      const index = settings.excluded.indexOf(available);
      if (index > -1)
        settings.excluded.splice(index, 1);
    } else
      settings.excluded.push(available);

    this.setUserSettings(settings);
  }

  private async availableRemove(available: string): Promise<void> {
    const settings = await this.getUserSettings();
    const index = settings.excluded.indexOf(available);
    if (index > -1)
      settings.excluded.splice(index, 1);
    this.setUserSettings(settings);
    await grok.dapi.files.delete(this.path + available);
    grok.shell.info(`File ${available} successfully deleted`);
  }

  private async getNewAvailable(): Promise<string> {
    return new Promise<string>((resolve, reject) => {
      DG.Utils.openFile({
        accept: this.ext,
        open: async (selectedFile) => {
          const path = selectedFile.name;
          const buffer = await selectedFile.arrayBuffer();
          await grok.dapi.files.write(this.path + `${selectedFile.name}`, (new Uint8Array(buffer)) as any);
          resolve(path);
        }
      });
    });
  }
}
