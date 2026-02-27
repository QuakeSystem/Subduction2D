# Subduction2D
## Contents
This repository contains the forward modeling developments of QuakeSystem.
It consists of the original JustRelax miniapp on 2D subduction, and the current version of the advanced model.
The advanced model is named SZU2019, after van Dinther et al. 2019 on A Secondary Zone of Uplift Due to Megathrust Earthquakes.

## How to run the original model
Run the original from the Subduction2D directory by running the command: 
julia --project Original/Subduction2D.jl

The project flag searches for the Project.toml file in the root. This contains information on all the dependencies.

## How to run the advanced model
The advanced model currently requires our QuakeSystem fork of JustRelax, instead of the default JustRelax.
You need to clone the fork, and then link julia's package manager to this clone (instead of the default).

1. Clone [our QuakeSystem fork of JustRelax](https://github.com/QuakeSystem/JustRelax.jl) to your local device, I recommend next to your Subduction2D directory.
2. Start a terminal and navigate to your local JustRelax directory
3. Enter the julia REPL and enter the package manager with the following commands:

julia (brings you to the julia REPL)
] (brings you to the package manager)
4. Now point the package manager to your new local JustRelax package version, rather than the JustRelax version on GitHub which might be not compatible anymore.
dev /full/path/to/JustRelax.jl 

The script currently needs you to switch to the bert_dev branch of JustRelax before you can run the advanced model.
1. In a terminal, navigate to your local JustRelax folder
2. Pull the latest updates from GitHub by running:
git pull
3. Change the git branch by using:
git checkout bert_dev

Now you can return to the Subduction2D repository and run the model.
* Run the advanced model locally using 'julia --project Subduction2D_SZU2019/Subduction2D_SZU2019.jl'
* Run the advanced model on eejit using 'sbatch run_eejit_sub2d.sbatch'