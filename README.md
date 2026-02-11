# Resource Management System Library Manual

The RMS Library provides a series of analysis pipeline (book) to facilitate codes sharing and reuse. Users could quickly access to all required binary and environments information. 

This manual is mainly documented for users directly using the RMS library. Advanced usage could be found in the `rmsp` manual. 



## Quick Start

A typical RMS library contains many different books (Analysis pipelines), with different editions (version of the pipelines). Within each book edition, it contains multiple selectable sections for users to run. These sections could be quickly accessed via 'bookmarks'. Download the book edition you like, and start running the analysis. 

### CLI solution

In the CLI solution, users should create a `parameters.json` file the include what commands to execute and the corresponding parameters.

```sh
# Normal RMS setup
rmstools setup_wizard -dbpath DBPATH -dbname DBNAME

# Register input files. 
rmstools register_files -dbpath DBPATH -dbname DBNAME -files full_path_input_file1 full_path_input_file2

# Download templates and resolve the template environment
rmstools download_templates -dbpath DBPATH -bookname BOOK1 -editionname EDITION1
conda env create --file="RMSLibrary/BOOK1/EDITION1/envs/conda/environment.yml"

# Execute template commands
rmstools execute_template_commands -dbpath DBPATH -dbname DBNAME -parameter_files parameters.json
```

Here is an example parameter file:

```json
{"commands":[
    {"bookname":"BOOK1",
     "editionname":"EDITION1",
     "bookmark":"BOOKMARK1",
     "parameters":{
         "conda_env":"CONDA_ENV",
         "Some_other_parameters":"Parameters",
     }
    }
]}
```



### Python solution

In the python solution, users could provide `args` and `kwargs` as parameters directly. 

```python
# Normal RMS setup
from rmsp import rmstools
rmstools.setup_wizard(dbpath, dbname)

# Download templates and resolve the template environment
rmstools.download_templates(dbpath, BOOK1, EDITION1)
import os; os.system('conda env create --file="RMSLibrary/BOOK1/EDITION1/envs/conda/environment.yml"') 
```



```python
from rmsp import rmscore, rmsbuilder, rmstemplate

# Connect to the RMS database
rms = rmscore.ResourceManagementSystem(dbpath + dbname, dbpath + "RMSResources/")
rmspool = rmsbuilder.RMSProcessWrapPool(rms, nthread)
rmsb = rmsbuilder.RMSUnrunTasksBuilder(rmspool)
rmstlib = rmstemplate.RMSTemplateLibrary(rmsb, dbpath + "RMSLibrary/")

# Register input files. 
rms.register_file("full_path_input_file1")
rms.register_file("full_path_input_file2")

# Execute the template commands
rmstlib.run(BOOK1, EDITION1, BOOKMARK1, args, kwargs)
rmstlib.execute_builder()
```



### GUI solution

This requires the installation of `rmsp-gui` package. 



## RMS Library Layout

````
RMSLibrary/
|___ Book1/              # A book is basically an analysis
|_______ Edition1/   
|___________ bin/            # Additional executable
|___________ doc/            # Documentations
|___________ envs/           # Environments
|___________ scripts/        # The actual RMS template scripts
|___________ manifest.json   # An index that organize all the analysis functions
|_______ Edition2/
|___________ ...
|___ Book2/
|_______ ...

````



## FAQ

- Why do you name the templates "book" and "edition" instead of  "program" and "version"?

  When naming the templates, we want to avoid potential confusion - the actual running programs and their version. 
  
- The program prevents me from running a certain workflow because output file already exists. 

  RMS forbids overwriting a file since this could cause significant confusion in the overall workflow. You need to delete all entries (from the database) that generating the existing output files (and all downstream files) before proceeding.

- I have problems in generating matplotlib figures using jupyter. An error occurs related to missing matplotlib_inline library.

  The problem may occur when you are running a bash function on another environment. If you are using jupyter, you may consider turning off matplotlib inline:

  ```
  %env MPLBACKEND=
  ```

  
