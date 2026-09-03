import * as grok from 'datagrok-api/grok';

/** A sample: row-level security, a natural key, replay-safe creates. */
@grok.decorators.entity({schema: 'lab', securityMode: 'row', businessKey: ['name'], idempotency: true})
export class Sample {
  @grok.decorators.column({required: true, unique: true, isName: true})
  name!: string;

  @grok.decorators.column({min: 0})
  count!: number;

  @grok.decorators.column()
  active!: boolean;

  @grok.decorators.column()
  measured_on!: Date;

  @grok.decorators.column()
  tags!: string[];

  @grok.decorators.column({type: 'user'})
  owner!: string;

  @grok.decorators.column({choices: ['new', 'done'], default: 'new'})
  status!: string;

  notes?: string;

  describe(): string {
    return `${this.name} (${this.count})`;
  }
}

@grok.decorators.entity({schema: 'lab', friendlyName: 'Plates', singularName: 'Plate'})
export class Plate {
  @grok.decorators.column({required: true, unique: true})
  barcode!: string;

  @grok.decorators.column({type: 'int', min: 1, max: 32})
  rows!: number;
}

@grok.decorators.entity({schema: 'lab', securityMode: 'master', delegate: 'plate_id', audit: false,
  businessKey: ['plate_id', 'well_position']})
export class PlateWell {
  @grok.decorators.column({ref: () => Plate, onDelete: 'cascade', required: true})
  plate_id!: string;

  @grok.decorators.column({ref: 'sample', onDelete: 'setnull'})
  sample_id?: string;

  @grok.decorators.column({name: 'well_position', required: true, friendlyName: 'Position'})
  position!: string;
}

@grok.decorators.entity({schema: 'lab', table: 'reader', description: 'Plate readers'})
export class Reader {
  @grok.decorators.column({required: true})
  name!: string;

  @grok.decorators.column({type: 'file'})
  config!: string;
}
