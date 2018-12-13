<h1>Triple-decomposition-analysis-for-instationary-air-water-flows</h1>

<p>The software is able to perform a triple decomposition of raw data collected in instationary air-water flows providing a wide range of air-water flow properties for the raw data as well as the decomposed signals (high-pass filtered signal, band-pass filtered signal and low-pass filtered signal). The software is suitable for two simultaneously sampled signals of a phase-detection intrusive probes including one double-tip probe or two single tip probes.</p>
<p>The software was developed by Dr Stefan Felder (UNSW Sydney) during his PhD project (Felder 2013) and has been used in several scientific publications. If you decide to use the software for post-processing of your data, please credit the creator of the work as:</p>
<p>Felder, S. and Chanson, H. (2014). "Triple Decomposition Technique in Air–Water Flows: Application to Instationary Flows on a Stepped Spillway." International Journal of Multiphase Flow, Vol. 58, pp. 139-153 & 3 videos.</p>
<p>The software has been already made public as a digital appendix in the PhD thesis of Felder (2013) and is still available for download here: https://espace.library.uq.edu.au/view/UQ:301329. Note that the readme file of the digital appendix is also provided as a PDF document for further documentation on the software (Note that the digital appendix refers to two further data processing tools including for a double-tip conductivity probe and two side-by-side single tip probes which are available as a separate repository).</p> 
<p>The software is able to analyse two simultaneously sampled raw Voltage signals of double-tip phase-detection intrusive probes (or two single tip probes). The input data may be in ASCII or binary format:</p>

- For ASCII raw data use “Triple_decomposition_analysis_txt_input_file”
- For binary raw data use “Triple_decomposition_analysis_binary_input_file”

The software is able to provide result files for the raw data and the decomposed data (including band-pass, high-pass and low-pass filtered data):
- Void fraction C (-) (raw data only)
- Bubble count rate F (Hz) (raw data only)
- Interfacial velocity V (m/s) 
- Turbulence intensity Tu (-) 
- Bubble/droplet chord time tch (s) (raw data only)
- Bubble/droplet chord length ch (m) (raw data only)
- Auto-correlation function Rxx (-)
- Cross-correlation function Rxy (-)
- Maximum cross-correlation coefficient (Rxy)max (-)
- Auto-correlation time scale Txx (s) 
- Cross-correlation time-scale Txy (s) 

<b>1 Content:</b>
<p>The software is written in Fortran Compaq. It has been successfully compiled with any Intel Fortran compiler including the latest version of Intel Fortran using Parallel Studio XE 2018 updates 2 or 3. The repository contains:</p>

- Full source code for triple decomposition analysis using raw data with ASCII format 
- Full source code for triple decomposition analysis using raw data with binary format
- Digital Appendix of Felder (2013) with additional documentation
- “Supplementary input data files for software”, which are needed to run the software comprising three input text files:
	- parameter.txt	: This file contains important input parameters (see Table DA-4 in the Digital Appendix of Felder (2013))
	- binary.txt :  The recorded input files are in binary format and must be of 8 character lengths. (Details regarding the required binary file format are available in the source code: lines 359 to 384.)
	- positions.txt :  The measurement locations must be of 4 character lengths and will appear in the summary result files.
The raw data were acquired with a LabVIEW data acquisition software, which was also developed by Felder (2013). (Further details on this acquisition software can be found in Felder (2013) and the Digital Appendix file.)

<b>2 How to run the code:</b>
- Compile the source code with a Fortran compiler (it works with any form of Intel or Compaq Fortran compilers) to create an executable file.
- Copy this file into the same folder as your ASCII or binary raw data together with the three supplementary input text files “parameter”, “binary” and “positions”.
- Run the executable file and follow the prompts in the software. The software will analyse all raw data files and will produce a range of result text files containing the wide range of air-water flow properties for the raw data and the decomposed signal components (see Digital Appendix of Felder (2013) for further details and documentation)

<b>3 Contact:</b>
<p>For feedback, questions and recommendations, please use the issue-section or contact the author via Email:
Stefan Felder, Senior Lecturer, Water Research Laboratory, School of Civil and Environmental Engineering, UNSW Sydney, Australia; Email: s.felder@unsw.edu.au </p>

<b>4. Selected References:</b>
<p>The software has been used in a number of studies and research projects. Some selected references are:</p>

Felder, S. (2013). “Air-Water Flow Properties on Stepped Spillways for Embankment Dams: Aeration, Energy Dissipation and Turbulence on Uniform, Non-Uniform and Pooled Stepped Chutes.” PhD thesis, The University of Queensland.

Felder, S. and Chanson, H. (2012). "Air-Water Flow Measurements in Instationary Free-Surface Flows: a Triple Decomposition Technique." Hydraulic Model Report No. CH85/12, School of Civil Engineering, The University of Queensland, Brisbane, Australia, 161 pages (ISBN 9781742720494).

Felder, S. and Chanson, H. (2014). "Triple Decomposition Technique in Air–Water Flows: Application to Instationary Flows on a Stepped Spillway." International Journal of Multiphase Flow, Vol. 58, pp. 139-153 & 3 videos.

Wang, H., Felder, S. and Chanson, H. (2014). "An Experimental Study of Turbulent Two-Phase Flow in Hydraulic Jumps and Application of a Triple Decomposition Technique." Experiments in Fluids, Vol. 55, No. 7, Paper 1775, 18 pages & 2 video movies.
