## Tuned exchange imaging: Can the filter exchange imaging pulse sequence be adapted for applications with thin slices and restricted diffusion?

Samo Lasič<sup>1,2*</sup>, Arthur Chakwizira<sup>3</sup>, Henrik Lundell<sup>2,4</sup>, Carl-Fredrik Westin<sup>5</sup>, and Markus Nilsson<sup>6</sup>

1. Department of Diagnostic Radiology, Lund University, Lund, Sweden
2. Danish Research Centre for Magnetic Resonance, Centre for Functional and Diagnostic Imaging and Research, Copenhagen University Hospital - Amager and Hvidovre, Copenhagen, Denmark
3. Department of Medical Radiation Physics, Lund University, Lund, Sweden
4. Magnetic Resonance Section, DTU Health Technologies, Technical University of Denmark, Kgs. Lyngby, Denmark
5. Department of Radiology, Brigham and Women's Hospital, Harvard Medical School, Boston, MA, United States
6. Department of Clinical Sciences Lund, Radiology, Lund University, Lund, Sweden


*Corresponding Author:\
Samo Lasič, PhD\
Department of Diagnostic Radiology, Lund, Sweden
SE-22185 Lund, Sweden
Email address: samo.lasic@med.lu.se


### Reference
If you use these resources, please cite:\
[Lasič et al. Tuned exchange imaging: Can the filter exchange imaging pulse sequence be adapted for applications with thin slices and restricted diffusion? NMR in Biomedicine. 2024;e5208. doi:10.1002/nbm.5208](http://doi.org/10.1002/nbm.5208)

### Additinal resources
[Multidimensional diffusion MRI repository](https://github.com/markus-nilsson/md-dmri)

### Overview
Filter exchange imaging (FEXI) with thin imaging slices requires strong crushers, which bias exchange rate estimates and confound exchange with restricted diffusion. The proposed modified tuned exchange imaging (TEXI) protocol provides accurate exchange rates regardless of slice thickness and restriction size, even with strong crushers.

This repository contains FEXI and TEXI pulse sequence parameters, results from Monte Carlo simulations and Matlab code used in the publication [Lasič et al., NMR Biomed. 2024;e5208. doi:10.1002/nbm.5208](http://doi.org/10.1002/nbm.5208).

Matlab scripts:

- `TEXI_make_seq.m`\
Make parameters for FEXI and TEXI sequences. TEXI waveforms have constant $bV_\omega$. Only the sequence parameters and info to generate effective gradient waveforms are saved in the "sequence" folder as `VT_*.mat` (for $V_\omega$ tuned TEXI) and `FEXI_*.mat` (for the original FEXI). The waveforms can be generated with `def_FEXI.m` and `def_sequence.m`.\
This script can save ADCs from frequency-domain GPA calculations and tune sequences based on the entire restricted diffusion spectrum for a given geometry and restriction size (not used in the manuscript). For a more general tuning, the sequence parameters and the corresponding ADC values are saved for a wide range of encodings with constant b-values (`b_*.mat` and `*_DvsR.mat` files). The tuned sequence parameters are then selected based on equal ADC values and the "size"-tuned sequence is saved as `ST_*.mat` file.

- `TEXI_make_data_and_seq_reduce_b_VT.m`\
Make new data and sequence files with a reduced set of b values for the VT sequence and save them as `*_reduceBrange.mat` files. 

- `TEXI_show_seq.m`\
Display sequence parameters:
  - crusher gradient magnitude
  - b-value
  - gradient magnitude
  - maximum gradient magnitude
  - pulse duration, $\delta$
  - mixing time, $t_m$
  - exchange-weighting time, $\Gamma$
  - spectral variance, $V_\omega$


- `TEXI_fit.m`\
Fit signals from simulations using different models. \
Fitting models used for the published paper:
  - the AXR model for the original FEXI, also included in the [Multidimensional diffusion MRI repository](https://github.com/markus-nilsson/md-dmri). See `fexi11_fit2data.m`.
  - the Ning model (Ning et al., J. Chem. Phys. 148, 2018) using the exchange-weighting function $h(k)$ and the fitting parameters: MD, V and k, where S0 = 1. See `Ning_h_norm_fit2data.m`.

  Additional fitting models:

  - the Ning model (Ning et al., J. Chem. Phys. 148, 2018) using the exchange-weighting function $h(k)$ and the fitting parameters: MD, V, k and S0. See `Ning_h_fit2data.m`.
  - the Ning model using the exchange-weighting time $\Gamma$ and the fitting parameters: MD, V, k and S0. See `Ning_Gamma_fit2data.m`.
  - the modified Ning model using the exchange-weighting function $h(k)$ and Gamma distributed diffusivity, with the fitting parameters: MD, V, k, S0. See `Ning_g_fit2data`.
  - the modified Ning model using the exchange-weighting function $h(k)$ and 3 cumulants, with the fitting parameters: MD, V, k, c3. See `Ning_h3_norm_fit2data`.
  - the modified Ning model using the exchange-weighting function $h(k)$ and Gamma distributed diffusivity, with the fitting parameters: MD, V, k. See `Ning_g_norm_fit2data`.
  - the modified Ning model using the exchange-weighting function $h(k)$ scaled with b-value, with the fitting parameters: MD, V, k, c3. See `Ning_h3b_norm_fit2data`.
  - the modified Ning model using the exchange-weighting function $h(k)$ and an offset due to the intra-compartmental variance of substrate mean diffusivities in directional averaging, with the fitting parameters: MD, V, k, $V_{intra}$. See `Ning_h_intra_norm_fit2data`.

- `TEXI_plot_res_list.m`\
  Combine plots with different substrates, crushers and sizes. Assumes each result is for a single substrate and size.

- `TEXI_show_q4.m`\
Display the fourth-order autocorrelation function of dephasing $q_4(t)$ for different sequence parameters.

- `TEXI_show_power_spectrum.m`\
Display the dephasing power spectra and cumulative attenuations for different sequence parameters (encoding waveforms).

- `TEXI_fit_precision.m`\
Estimate fit precision for selected data.  

- `TEXI_precision_plot.m`\
Plot results from the precision estimation.
  
- `TEXI_Ning_FEXI_Gauss.m`\
Simulate blood-brain barrier exchange with Gaussian IVIM diffusion contrast in FEXI with crushers. Exchange rates are estimated with the original AXR model and with the Ning model.

- `TEXI_res_list_Gauss.m`\
Plot results for the blood-brain barrier exchange estimated with the original AXR fit and the Ning model.


Additional scripts used for various tests in this project.

- `TEXI_plot_res_single.m`\
  Inspect raw signals, fitted signals and fit parameters for selected substrate/size, crusher and exchange rate.  

- `TEXI_make_data_and_seq_with_sel_tm.m`\
A sub-sampling script. Make new data and sequence files with selected subsets of mixing times and save them as `*_sub_tm.mat` files.

- `TEXI_make_data_and_seq_without_SDE.m`\
A sub-sampling script. Make new data and sequence files without the SDE in VT and ST sequences and save them as `*_noSDE.mat` files.

- `TEXI_add_b0_data_and_seq.m`\
Make new data and sequence files by additing S = 1 @ b = 0 and save them as `*_b0.mat` files.
