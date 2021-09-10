# SSA_project
Python script for Secondary Structure Assignement base on hydrogen bonds

`git clone https://github.com/thomas-bailly/SSA_project.git`
## Requirement
### Python
For windows users who don't have Python you can download it at this address: https://docs.conda.io/en/latest/miniconda.html

For MacOSX and systems Linux users if you don't have Miniconda3, I recommend you to download it with the link above.

For the install on Linux, enter in a shell the command below : `bash Miniconda3-latest-Linux-x86_64.sh` 

For the install on MacOSX, enter in a shell the command below : `bash Miniconda3-latest-MacOSX-x86_64.sh`
### Librairies
- Numpy
- Scipy
- Biopython

If you have conda, enter the command below: `conda env create -f env.yml` in the SSA_project directory, then: `conda activate ssa_project` for activate the environment with the librairies. enter `conda deactivate` for exit the environment.

If you don't have conda, enter the command: `python3 -m pip install --user numpy scipy biopython`
### Pdb file with hydrogens
The script needs the hydrogens to be in the pdb file to work properly. For add the hydrogens you can go in [Molprobity](http://molprobity.biochem.duke.edu/), upload your file or fetch the pdb code and add hydrogens. For Linux users you can donwload the reduce program wit the command :

`conda install -c ostrokach-forge reduce`

#### Reduce
If you just enter reduce in a shell, you can see the help.

here's an example of use: `reduce -BUILD -i 2obv.pdb > 2obvFH.pdb`

## Usage

Go to the SSA_project directory and enter: `python3 ./src/ssa.py help` for obtain the help. 

For run the script enter a command like this: `python3 ./src/ssa.py ./data/1bat1FH.pdb A`

You can redirect the result with > : `python3 ./src/ssa.py ./data/1bta1FH.pdb A > ./result/1bta_ssa.txt`

The result for the file 3h7hFH.pdb is already present in the subdirectory result.
