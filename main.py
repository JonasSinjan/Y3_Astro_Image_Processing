from astropy.io import fits

with fits.open("A1_mosaic.fits") as hdulist:
    for key,val in hdulist[0].header.items():
        print(f"{key},{val}")