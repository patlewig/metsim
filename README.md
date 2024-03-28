# Metabolism Simulation Framework (MetSim)
==============================

Code repository supporting manuscript: Groff L, Williams A, Shah I, Patlewicz G*. 2024. MetSim: Integrated Programmatic Access and Pathway Management for Xenobiotic Metabolism Simulators. Chem. Res. Toxicol *in press*

Associated data files for the manuscript are accessible at <https://doi.org/10.23645/epacomptox.25463926> as a gz zipped tarball. When cloning the repository, a data folder structure as shown below needs to be created. Use  tar xvzf metsim_data.tar.gz.

MetSim comprises three main components: a graph based schema to allow metabolism information to be harmonised for subsequent analysis; an application programming interface (API) for 4 metabolic simulators: BioTransformer, the OECD Toolbox, EPA's Chemical Transformation Simulator (CTS) and TIssue Metabolism Simulator (TIMES); and functions to help evaluate simulator performance.

## Prequisites

Preprocessing and standardisation of chemical information rely on the API calls from the EPA CompTox Chemicals Dashboard and the EPA Cheminformatics Modules. Further details on both these web applications can be found at <www.comptox.epa.gov>. Functions to retrieve chemical information on the basis of SMILES are provided in the utilities folder of the metsim package.

**CTS** predictions can be made using the web-based API at <https://qed.epa.gov/cts/rest>. For the other 3 simulators, local installations are required.

**BioTransformer 3.0** is an open access software tool available from <https://biotransformer.ca/> 
Further information about BioTransformer have been published here: 
Wishart DS, Tian S, Allen D, Oler E, Peters H, Lui VW, Gautam V, Djoumbou Feunang Y, Greiner R, Metz TO; BioTransformer 3.0 – A Web Server for Accurately Predicting Metabolic Transformation Products [Submitted in Nucleic Acids Research, Webserver Issue.Apr.2022] 

Djoumbou Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R, and Wishart DS; BioTransformer: A Comprehensive Computational Tool for Small Molecule Metabolism Prediction and Metabolite Identification; Journal of Cheminformatics; 2019; Journal of Cheminformatics 11:2; DOI: 10.1186/s13321-018-0324-5 

**The OECD QSAR Toolbox** is a freely available tool downloadable from <https://qsartoolbox.org/> upon registration. Once installed, the local web server instance needs to be started. 

**TIMES** is a commercial tool, licenced from the Laboratory of Mathematical Chemistry, University As. Zlatarov, Burgas, Bulgaria. Further information can be found here <http://oasis-lmc.org/products/software/times.aspx>. Whilst there are means of accessing TIMES programmatically, the functions in metsim are limited to processing the report output produced from TIMES after running the rat metabolism simulators.
    
## Easy start

Running the notebooks in this repository requires Python, Anaconda (or Miniforge), Jupyter and some additional configuration. 
1. Install Python 3, anaconda/conda and Jupyter Lab
2. Clone the repository
3. Download the data files from the figshare address and create the relevant data files directory
4. Create a conda or virtual environment to install relevant libraries - pandas, matplotlib, seaborn, rdkit, numpy are needed as minimum.
5. Add the environment as a kernel to jupyter-lab

Further details are provided in the relevant notebooks.


Project Organization
------------
   
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external         <- MetSim Tool details and Performance table
    │   ├── interim           <- Intermediate data that has been transformed.
    │   ├── processed      <- Processed MetSim json files for the 112 chemical database
    │   └── raw                <- ClassyFire outputs for all 112 chemicals used in the analysis
    │
    │
    ├── notebooks          <- Jupyter notebooks to walk through each tool and the analysis supporting the manuscript
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── metsim          <- Source code for use in this project.
    │   ├── utl           <- Scripts to process and gather chemical specific data
    │   ├── sim           <- Scripts to run each of the MetSim tools
   
 
--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
