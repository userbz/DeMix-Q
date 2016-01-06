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

- Load the DeMix.toppas pipeline in the TOPPAS program.

- Prepare centroid MS1 spectra in mzML format and load in __"Node 1"__ of the DeMix.toppas pipeline.

- Convert corresponding MS/MS identification lists into OpenMS-compatible __".idXML "__ format. MS-GF+ resultant mzid files can be converted by the python script: __mzid2idxml_theomass.py__

- Modify __"Node 3"__ in the pipeline by selecting the converted idXML files.

- Modify __"Node 4"__ in the pipeline by choosing one of the idXML file (the best run) as the reference run. RT scales of other runs will be re-calibrated to the scale of the reference run after alignment.

- Execute the pipeline in TOPPAS.

- Collect output files from processes of TextExporter (output 16) and EICExtractor (output 20), and do the post processing as shown in the example __iPython notebook (FIXME: not uploaded yet, sorry)__.

---
__To Do:__
* Example data
* iPython notebook.

---
##### Reference
<ol> <li> Zhang, B., Kall, L. & Zubarev, R. A., (2016). DeMix-Q: Quantification-centered Data Processing Workflow. Molecular & Cellular Proteomics, mcp.O115.055475  
http://www.ncbi.nlm.nih.gov/pubmed/26729709 </li>

<li>Zhang, B., Pirmoradian, M., Chernobrovkin, A., & Zubarev, R. A. (2014). DeMix Workflow for Efficient Identification of Co-fragmented Peptides in High Resolution Data-dependent Tandem Mass Spectrometry. Molecular & Cellular Proteomics : MCP. doi:10.1074/mcp.O114.038877
http://www.ncbi.nlm.nih.gov/pubmed/25100859</li></ol>
