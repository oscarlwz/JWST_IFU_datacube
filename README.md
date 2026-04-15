This is a collection of script used for the paper Liu et al 2026.

To run the scripts, please download all files within this folder.

Jupyter Notebooks:
nrs_ifu_template_naturepaper.ipynb is an example to reduce the IFU data cube in this project.
qso46_pyqsofit_submission.ipynb is used for all analyses and figure production.
qso46_pyqsofit_2o3_radialbin_submission.ipynb is used to plot the radial profiles of outflows.

Data:
The JWST 1-D spectra of all 27 quasars are in the folder all1dtxt_27sources.
For each file, the three columns are the wavelength (mircon), flux density (cgs units) and errors (cgs units).
All required data for running qso46_pyqsofit_submission.ipynb are located at as a zip file:
https://zenodo.org/records/19027791?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjU0NTczODI0LTZjYjUtNDExZS05OWNlLWZmMDlmZWU1M2I1ZiIsImRhdGEiOnt9LCJyYW5kb20iOiI5YzI4MDFiMDc1OWMxNTJkYTI5YzdkODhhMTQ1YmY0YSJ9.yYdo672FNHbT7UCoHKBTmkt0-urmdNkeu1cBHRkcjtVMPMyHDr71FS1hmgpB0RDXa0ykKj_aAeJXxn-g5KHm8w

Specifically, the Shen2016_table2.fits and spectra in the folder Shen_oirspec are obtained from https://iopscience.iop.org/article/10.3847/0004-637X/817/1/55 and should obey the copyright therein.

dr16q_prop_May01_2024.fits.gz is obtained from https://iopscience.iop.org/article/10.3847/1538-4365/ac9ead 

ERQ_Perrotta_table3.txt is compiled from https://academic.oup.com/mnras/article/488/3/4126/5538862?login=false

Software dependence:
PyQSOFit can be found at: https://github.com/legolason/PyQSOFit/tree/master
