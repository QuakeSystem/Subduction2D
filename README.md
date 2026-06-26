# Subduction2D
## Contents
This repository contains the forward modeling developments of QuakeSystem.
It consists of the original JustRelax miniapp on 2D subduction, and the current version of the advanced model.
The advanced model is named SZU2019, after van Dinther et al. 2019 on A Secondary Zone of Uplift Due to Megathrust Earthquakes.

## How to run the original model
Run the original from the Subduction2D directory by running the command: 
julia --project Original/Subduction2D.jl

The project flag searches for the Project.toml file in the root. This contains information on all the dependencies.

## How to install the advanced model
The advanced model currently requires our QuakeSystem fork of JustRelax, instead of the default JustRelax.
You need to clone the fork, and then link julia's package manager to this clone (instead of the default).

1. Clone [our QuakeSystem fork of JustRelax](https://github.com/QuakeSystem/JustRelax.jl) to your local device, I recommend next to your Subduction2D directory.
2. Start a terminal and navigate to your local JustRelax directory
3. Enter the julia REPL and enter the package manager with the following commands:

julia (brings you to the julia REPL)
] (brings you to the package manager)
4. Now point the package manager to your new local JustRelax package version, rather than the JustRelax version on GitHub which might be not compatible anymore.
dev /full/path/to/JustRelax.jl 

Now you can return to the Subduction2D repository and run the model.

## How to run the advanced model
The advanced model is found at Subduction2D/Subduction2D_SZU2019/ and consists of:
- Subduction2D.jl (main script)
- Subduction2D_rheology.jl (material parameters)
- Subduction2D_setup.jl (initial geometry and temperature)

To simply run the model you do have to adjust two items:
1. You _have to_ configure whether you run on CUDA (GPU's) or not by configuring, in line 3 of the main. ```const isCUDA = true``` or ```const isCUDA = false```
2. The output directory is currently hardcoded in the main script, in function 'prepare_visualisation'. 
```figdir   = "Subduction2D_SZU2019/data/Subduction2D_JRv0.6.0/$version"```

* Run the advanced model locally using 'julia --project Subduction2D_SZU2019/Subduction2D.jl'
* Run the advanced model remotely on Eejit, using GPU, by running ```sbatch run_eejit_sub2d.sbatch``` . Do adjust your path-to-julia.
* Automated workflow: Run the advanced model on eejit by running 'bash new_run.sh <name_of_version>'. This stores your input files and the slurm script to the output directory, and from there _calls the slurm script_. This is done to automate the process and to prevent issues with queued jobs using the wrong input files.

## Assumed structure of output folders
Subduction2D_SZU2019/data/<current_main_version_of_justrelax>/<sub2d_version>
examples:
- Subduction2D/Subduction2D_SZU2019/data/Subduction2D_JRv0.5.1/v0.350_40myr_halfspace
- Subduction2D/Subduction2D_SZU2019/data/Subduction2D_JRv0.6.0/v0.358_old_density_40Myr_heatprod

## Useful scripts in utils/
utils/ has convenient scripts to improve the workflow.
make_all_movies.sh uses make_movies_in_directory on all sub2d versions in an output directory. 
It bundles all .png files made of a run (full domain, and zoomed in) and makes an mp4 out of them.
They're stored in /movies/

Bert has it configured currently to an alias is .bashrc. Then you can just run 'makemovies' and then 'collectmovies'
```alias makemovies='cd Subduction2D_SZU2019/data/Subduction2D_JRv0.6.0/ && bash ../../../utils/make_all_movies.sh 20 && cd ../../../' ```
Puts the mp4 movies in, for example: Subduction2D_SZU2019/data/Subduction2D_JRv0.6.0/v0.358_new_density_40Myr/movies
```alias collectmovies='bash /scratch/tectonics/bert/Subduction2D/utils/collect_movies.sh' ```
Then you will find all movies in Subduction2D_JRv0.6.0 added to Subduction2D_SZU2019/data/subduction_movies/.