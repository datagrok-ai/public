import sharp from 'sharp';
import { readdirSync, mkdirSync, copyFileSync, renameSync, unlinkSync } from 'fs';
import { join, extname, basename } from 'path';

const INPUT_DIR  = './static/docusaurus_img/team';
const BACKUP_DIR = './static/docusaurus_img/team/_originals';
const SIZE = 400;

mkdirSync(BACKUP_DIR, { recursive: true });

const files = readdirSync(INPUT_DIR).filter(f => {
  const ext = extname(f).toLowerCase();
  return (ext === '.jpg' || ext === '.jpeg' || ext === '.png') && f !== 'no-avatar.png';
});

console.log(`Processing ${files.length} images...`);

for (const file of files) {
  const src  = join(INPUT_DIR, file);
  const bak  = join(BACKUP_DIR, file);
  const dest = join(INPUT_DIR, basename(file, extname(file)) + '.jpg');
  const tmp  = dest + '.tmp.jpg';

  // Back up original (skip if already backed up)
  try { copyFileSync(src, bak); } catch {}

  try {
    await sharp(src)
      // flatten handles transparent PNGs — fills with light neutral grey
      .flatten({ background: { r: 242, g: 242, b: 242 } })
      // smart crop: 'attention' finds the visually dominant region (faces)
      .resize(SIZE, SIZE, { fit: 'cover', position: 'attention' })
      // convert to greyscale
      .grayscale()
      // modest sharpening to counter softness from resizing
      .sharpen({ sigma: 0.6 })
      // save as high-quality JPEG
      .jpeg({ quality: 90, mozjpeg: true })
      .toFile(tmp);

    // replace original with processed version
    try { unlinkSync(dest); } catch {}
    renameSync(tmp, dest);

    console.log(`  ✓ ${file}`);
  } catch (err) {
    console.error(`  ✗ ${file}: ${err.message}`);
  }
}

console.log(`Done. Originals saved to ${BACKUP_DIR}`);
