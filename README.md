# Metabolism Simulation Framework (MetSim)
==============================

Code repository supporting manuscript: Groff L, Williams A, Shah I, Patlewicz G*. 2024. MetSim: Integrated Programmatic Access and Pathway Management for Xenobiotic Metabolism Simulators. Chem. Res. Toxicol *in press*  <https://doi.org/10.1021/acs.chemrestox.3c00398>

Associated data files for the manuscript are accessible at <https://doi.org/10.23645/epacomptox.25463926> as a gz zipped tarball. When cloning the repository, a data folder structure as shown below needs to be created. Use  tar xvzf metsim_data.tar.gz.

MetSim comprises three main components: a graph based schema to allow metabolism information to be harmonised for subsequent analysis; an application programming interface (API) for 4 metabolic simulators: BioTransformer, the OECD Toolbox, EPA's Chemical Transformation Simulator (CTS) and TIssue Metabolism Simulator (TIMES); and functions to help evaluate simulator performance.

## Prequisites

Preprocessing and standardisation of chemical information rely on the API calls from the EPA CompTox Chemicals Dashboard and the EPA Cheminformatics Modules. Further details on both these web applications can be found at <https://www.comptox.epa.gov>. Functions to retrieve chemical information on the basis of SMILES are provided in the utilities folder of the metsim package.

**CTS** predictions can be made using the web-based API at <https://qed.epa.gov/cts/rest>. For the other 3 simulators, local installations are required.

**BioTransformer 3.0** is an open access software tool available from <https://biotransformer.ca/> 
A local download of the runnable Java Archive (JAR) folder is required, and the current working directory must be set to this folder in order for the BioTransformer functions within MetSim to work. 
The latest version of the BioTransformer executable JAR directory is available at the following link: <https://bitbucket.org/wishartlab/biotransformer3.0jar/src/master/> 
The June 15, 2022 specific version of BioTransformer JAR directory used in this manuscript is available at the following link: <https://bitbucket.org/wishartlab/biotransformer3.0jar/src/cc4006a/>

A Dockerfile is supplied in the main directory of the repo from which a Docker image can be built. An example notebook showing how to create the output from the MetSim functions is also provided.

Further information about BioTransformer have been published here: 
Wishart DS, Tian S, Allen D, Oler E, Peters H, Lui VW, Gautam V, Djoumbou Feunang Y, Greiner R, Metz TO; BioTransformer 3.0 – A Web Server for Accurately Predicting Metabolic Transformation Products [Submitted in Nucleic Acids Research, Webserver Issue.Apr.2022] 

Djoumbou Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R, and Wishart DS; BioTransformer: A Comprehensive Computational Tool for Small Molecule Metabolism Prediction and Metabolite Identification; Journal of Cheminformatics; 2019; Journal of Cheminformatics 11:2; DOI: 10.1186/s13321-018-0324-5 

**The OECD QSAR Toolbox** is a freely available tool downloadable from <https://qsartoolbox.org/> upon registration. Once installed, the local web server instance needs to be started.  With Toolbox version 4.5 installed:
1. Open Command Prompt using CMD in start menu
2. cd to C:\Program Files (x86)\QSAR Toolbox\QSAR Toolbox 4.5\Toolbox Server\Bin
3. Launch the WebAPI by entering: LMC.Toolbox.WebAPI.exe --urls="http://0.0.0.0:{tb_port}"
Note: {tb_port} is any user-defined integer, as long as that port number is not already in use by another process in the system (valid example, --urls="http://0.0.0.0:16384"). After a few moments, the command line console should print several messages, one of which will say "Now Listening on: http://0.0.0.0:{tb_port}"
4. Leave this command prompt window open for the Toolbox functions to work in Jupyter-lab.


**TIMES** is a commercial tool, licenced from the Laboratory of Mathematical Chemistry, University As. Zlatarov, Burgas, Bulgaria. Further information can be found here <http://oasis-lmc.org/products/software/times.aspx>. Whilst there are means of accessing TIMES programmatically, the functions in metsim are limited to processing the report output produced from TIMES after running the rat metabolism simulators. After loading a TXT containing CASRN, Names and SMILES (or alternative a sdf V2000), load the relevant metabolic simulator and click on All in the **Predict** tab. Once completed, go to the **Report** tab and click on **Tree Report** to generate the relevant output file. This can be copied/pasted into Excel or equivalent or Saved as a Text file directly. The file will be read using the pd.read_csv(file) using the Pandas library as part of the MetSim function.
    
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
