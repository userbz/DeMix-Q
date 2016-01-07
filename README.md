### DeMix-Q: Quantification-centered LC-MS/MS data processing
------
__Dependencies:__  
* __Python 2.7__
* __OpenMS 2.0__
* numpy  
* pandas
* lxml  
* pymzml http://pymzml.github.io/  
* pyteomics http://pythonhosted.org//pyteomics/  

Optional:
* msconvert (ProteoWizard) http://proteowizard.sourceforge.net/index.shtml  
* MS-GF+ http://omics.pnl.gov/software/ms-gf


---
__Procedures:__   

- Copy the consensus2edta.ttd config file to the OpenMS installation path (e.g. C:\Program Files\OpenMS-2.0\share\OpenMS\TOOLS\EXTERNAL) in order to enable external the external tools for XIC extraction.

- Load the __DeMix.toppas__ pipeline in the TOPPAS program.

- Prepare centroid MS1 spectra in mzML format and load in __"Node 1"__ of the DeMix.toppas pipeline.

- Convert corresponding MS/MS identification lists into OpenMS-compatible __".idXML "__ format. MS-GF+ resultant mzid files can be converted by the python script: __mzid2idxml_theomass.py__

- Modify __"Node 3"__ in the pipeline by selecting the converted idXML files.

- Modify __"Node 4"__ in the pipeline by choosing one of the idXML file (the best run) as the reference run. RT scales of other runs will be re-calibrated to the scale of the reference run after alignment.

- Modify __"Node 17"__ in the pipeline, set the __converter_script__ option to point to the python script: __consensus2edta.py__

- Execute the pipeline in TOPPAS.

- Collect output files from processes of TextExporter (output 16) and EICExtractor (output 20), and do the post processing as shown in the example __iPython notebook__ in the script folder.


__Note: __  
The example workflow is tested with high resolution Orbitrap (FTMS 70,000) datasets. You may need to adjust the parameters to fit with your dataset.

Map RT alignment can be done without MS/MS identifications, or done with other tools then convert into the OpenMS-compatible format (i.e. trafoXML).

For a large dataset (e.g. >10 samples), feature detection and linking may require considerable computing resources: CPU, RAM and disk space. 

---

__Download the example data set__   
https://ki.box.com/s/15kv73pu2cvi4gwakzd6kdxcwch7dbl0


---
##### Reference
<ol> <li> Zhang, B., Käll, L. & Zubarev, R. A., (2016). DeMix-Q: Quantification-centered Data Processing Workflow. Molecular & Cellular Proteomics, mcp.O115.055475  
http://www.ncbi.nlm.nih.gov/pubmed/26729709 </li>

<li>Zhang, B., Pirmoradian, M., Chernobrovkin, A., & Zubarev, R. A. (2014). DeMix Workflow for Efficient Identification of Co-fragmented Peptides in High Resolution Data-dependent Tandem Mass Spectrometry. Molecular & Cellular Proteomics : MCP. doi:10.1074/mcp.O114.038877
http://www.ncbi.nlm.nih.gov/pubmed/25100859</li></ol>
