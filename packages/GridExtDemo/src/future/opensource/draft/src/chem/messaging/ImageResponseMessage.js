export class ImageResponseMessage
{
  constructor(bitmap)
  {
      this.m_bitmap = bitmap;
  }

  getBitmapImage() {return this.m_bitmap;}
}