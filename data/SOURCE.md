# Data sources

## TESS light curve

- File: `tess2019058134432-s0009-0000000035516889-0139-s_lc.fits`
- Archive: Mikulski Archive for Space Telescopes (MAST), TESS SPOC light-curve product
- TESS sector: 9
- TIC target ID: 35516889
- MAST observation ID: 62892306
- MAST data URI: `mast:TESS/product/tess2019058134432-s0009-0000000035516889-0139-s_lc.fits`
- Exact download URL: <https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS%2Fproduct%2Ftess2019058134432-s0009-0000000035516889-0139-s_lc.fits>
- Collection DOI: [10.17909/t9-nmc8-f686](https://doi.org/10.17909/t9-nmc8-f686) (TESS 2-minute light curves, all sectors; sector 9 used here)
- Retrieved: 2026-08-15
- SHA-256: `70a40881930e7b374632a7afabb34ba2b18ff81841f6a622f03cf1d0268a06c0`

The FITS file is stored unmodified. The analysis reads `TIME`, `PDCSAP_FLUX`,
`PDCSAP_FLUX_ERR`, and `QUALITY`. PDCSAP flux is the SPOC light curve with common
instrumental trends removed and aperture/crowding corrections applied; this does
not make it free of residual stellar or instrumental systematics.

## System parameters

- File: `system_parameters.csv`
- Service: NASA Exoplanet Archive TAP, `pscomppars` table
- Exact query: <https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name%2Chostname%2Cra%2Cdec%2Cpl_orbper%2Cpl_tranmid%2Cpl_trandur%2Cpl_rade%2Cpl_bmasse%2Cpl_eqt%2Cpl_orbsmax%2Csy_dist%2Csy_tmag%2Cst_teff%2Cst_rad%2Cst_mass%2Cdisc_year%2Cdiscoverymethod%2Cdisc_refname%2Cdisc_pubdate%2Cdisc_facility+from+pscomppars+where+pl_name%3D%27WASP-19+b%27&format=csv>
- Retrieved: 2026-08-15

The saved row is the input actually used by `scripts/analyze_transit.py`; the
analysis does not query a changing live service at run time.


## Additional TESS sectors for robustness analysis

All are unmodified standard-cadence SPOC light curves from the same [MAST TESS collection](https://doi.org/10.17909/t9-nmc8-f686).

- Sector 9: `tess2019058134432-s0009-0000000035516889-0139-s_lc.fits` (1,848,960 bytes)
  - MAST URI: `mast:TESS/product/tess2019058134432-s0009-0000000035516889-0139-s_lc.fits`
  - SHA-256: `70a40881930e7b374632a7afabb34ba2b18ff81841f6a622f03cf1d0268a06c0`
- Sector 62: `tess2023043185947-s0062-0000000035516889-0254-s_lc.fits` (1,880,640 bytes)
  - MAST URI: `mast:TESS/product/tess2023043185947-s0062-0000000035516889-0254-s_lc.fits`
  - SHA-256: `c96bb25dadab40209390793713ad8ca6a579a0b224b06b9889f689cce555599a`
- Sector 63: `tess2023069172124-s0063-0000000035516889-0255-s_lc.fits` (1,941,120 bytes)
  - MAST URI: `mast:TESS/product/tess2023069172124-s0063-0000000035516889-0255-s_lc.fits`
  - SHA-256: `8d61f6d1c3bc90f061ae4ce37ce9509ed52f4836022576ef6581510061a5601f`
