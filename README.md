Metabolism Simulation Framework (MetSim)
==============================

Code repository supporting manuscript entitled MetSim: Integrated Programmatic Access and Pathway Management for Xenobiotic Metabolism Simulators

Project Organization
------------

   
   
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── metsim             <- Source code for use in this project.
    │   ├── __init__.py    <- Makes metsim a Python module
    │   │
    │   ├── util           <- Scripts to process and gather chemical specific data
    │   │   └── make_dataset.py
    │   │
    │   ├── sim           <- Scripts to run each of the MetSim tools
    │   │   │                 
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
   
   


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
