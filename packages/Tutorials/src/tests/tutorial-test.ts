import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { before, category, test, assure, expect } from '@datagrok-libraries/test/src/test';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import { Track } from '@datagrok-libraries/tutorials/src/track';


class TestTutorial extends Tutorial {
  name = 'Test tutorial';
  description = 'Test description';
  icon = '';
  steps = 1;

  protected _run(): Promise<void> {
    throw new Error('Method "_run" not implemented.');
  }
}

category('Tutorials', () => {
  const tutorial = new TestTutorial();
  const track = new Track('Test track', [tutorial], 'test-track-url');

  test('Create a track', async () => {
    expect(track.name, 'Test track');
    expect(track.helpUrl, 'test-track-url');
    expect(track.tutorials.length, 1);
    expect(track.tutorials[0].name, tutorial.name);
  });

  test('Run a tutorial', async () => {
    try {
      await tutorial.run();
    } catch (e) {
      expect(e instanceof Error ? e.message : e, 'Method "_run" not implemented.');
    }
  });
});
