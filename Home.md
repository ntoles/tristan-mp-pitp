Welcome to the tristan-mp-pu wiki!

# Getting started
## Setting up your personal computer to access Princeton's cluster

* [First time logging in to the Princeton Computer cluster](Setting-up-ssh-keys#first-time-accessing-perseus-andor-tigressdata)

  * [Using VPN](Setting-up-ssh-keys#setting-up-a-vpn)
  * [Using SSH keys](Setting-up-ssh-keys#ssh-keys)

    * [Generating your own SSH keypair](Setting-up-ssh-keys#creating-ssh-keys-on-your-local-computer)
    * [Putting your key on the remote machine](Setting-up-ssh-keys#putting-your-ssh-keys-on-the-cluster)
    * [Using your keys to log in](Setting-up-ssh-keys#logging-in-using-your-ssh-key)

## Running your first Tristan-MP Job on Princeton's cluster
* [Using Perseus](Logging-in-to-perseus)
  * [Logging into Perseus and configuring your environment](Logging-in-to-perseus#logging-in-to-perseus)
  * [Familiarizing yourself with the Filesystem](Logging-in-to-perseus#princeton-clusters-file-system)

* [Downloading and Compiling Tristan-MP](Downloading-and-Compiling-Tristan)

  * [Downloading Tristan-MP from github](Downloading-and-Compiling-Tristan#downloading-tristan-mp-from-github)
  * [Starting your own branch](Downloading-and-Compiling-Tristan#creating-your-own-tristan-mp-branch-and-user-files)
  * [Compiling Tristan-MP](Downloading-and-Compiling-Tristan#compiling)


* [Running your first Tristan-MP Simulation](running-your-first-tristan-mp-simulation)
  * [Getting ready](running-your-first-tristan-mp-simulation#getting-ready)
  * [Editing the input file](running-your-first-tristan-mp-simulation#editing-the-input-file)
  * [Submitting the job](running-your-first-tristan-mp-simulation#submitting-jobs)

## Analyzing and Visualizing your simulations

* [Tigressdata2](Logging-in-to-tigressdata2)

  * [Starting a VNC server](Logging-in-to-tigressdata2#starting-a-vnc-server-on-tigressdata2)
  * [Connecting to Tigressdata2](Logging-in-to-tigressdata2#connecting-to-tigressdata2)
  * [Closing a VNC Session](Logging-in-to-tigressdata2#closing-a-vnc-session)
  * [Enabling Tigressdata2's GPUs](Logging-in-to-tigressdata2#enabling-tigressdata2s-gpus)
  * [Accessing Jupyter Notebook using ssh](Logging-in-to-tigressdata2#jupyter-notebook-through-ssh)
* [Writing your own python analysis scripts](Writing-your-own-python-scripts)

# Example Simulations

* [Weibel instability](Sample-Weibel-run)

* [A Shock](Sample-shock-run)

# Tristan-MP Features and Code Structure

* [Overview](Code-Features#overview)

  * [General Code Structure](Code-Features#general-code-structure)

    * [Main Files](Code-Features#main-files)
    * [Main Algorithms](Code-Features#main-algorithms)
    * [Outputs](Code-Features#outputs)
    * [Initialization](Code-Features#initialization)
    * [Communication](Code-Features#communication)
    * [Auxiliary Modules](Code-Features#auxiliary-modules)
    * [User Modules](Code-Features#user-module)

  * [Input Files](Input-Files)

    * [Input File Structure](Input-Files#input-file-structure)
    * [Sample Input File](Input-Files#sample-input-file)

  * [User Files](User-Files)

# [User Manual](User-Manual)