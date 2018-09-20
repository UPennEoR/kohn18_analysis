import glob
import os
import shutil

abscal_file = 'gc.2457549.uvcRP.abscal.uvfits'

msfile = abscal_file.strip('uvfits') + 'MS'
if os.path.exists(msfile):
    shutil.rmtree(msfile)
importuvfits(vis=msfile, fitsfile=abscal_file)

image_fn = msfile[:-3] + '_abscal'
if os.path.exists(image_fn + '.image'):
    image_files = glob.glob('{}*'.format(image_fn))
    for f in image_files:
        shutil.rmtree(f)
clean(msfile, msfile[:-3] + '_abscal', niter=0, weighting='briggs',
      robust=0, imsize=[512, 512], pbcor=False, cell=['500 arcsec'],
      mode='mfs', nterms=1, spw='0:150~900', stokes='IQUV',
      interactive=True, npercycle=5, threshold='0.1mJy/beam')
