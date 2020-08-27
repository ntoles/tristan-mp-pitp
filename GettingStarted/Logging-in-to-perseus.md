---
title: Logging into Perseus
has_children: False
parent: Getting Started
nav_order: 1
---

## First time accessing Perseus and/or Tigressdata

The first time you try to access Perseus and/or Tigressdata from you local machine you must be on Princeton’s network or use a VPN connection.

If you’re off-campus first you’ll need to install a VPN to connect to the cluster. VPN access to Princeton is well documented and you will find information about how to do it here.. Follow the instructions on how to install sonic wall connect and VPN into the local Princeton network.

## Perseus Information
[Perseus](https://researchcomputing.princeton.edu/systems-and-services/available-systems/perseus) is the supercomputer where we run most of our Tristan-MP development runs. Perseus is a 320 node Dell Beowulf cluster. All compute nodes in this cluster are connected via an Infiniband network designed for high speed and low latency to enable excellent performance for tightly coupled MPI codes. Each node has 28-core Broadwell processors and 128 GB of RAM.

### General Guidelines

The head node, perseus, should be used for interactive work only, such as compiling programs, and submitting jobs as described below. No jobs should be run on the head node, other than brief tests that last no more than a few minutes. Where practical, we ask that you entirely fill the nodes so that CPU core fragmentation is minimized. For this cluster, perseus, that means multiples of 28 cores.

*Please remember that these are shared resources for all users.*

## Logging into Perseus.
Connect to perseus. Type in terminal:
```bash
ssh netID@perseus.princeton.edu
```
If you have not already done so, it is highly recommended that you copy and source Anatoly’s `.bashrc` file into your perseus account.
```bash
cp ~anatoly/.bashrc ./
source ~/.bashrc
```    
If instead, you prefer manually loading the modules necessary to compile Tristan-MP, type:
```bash
module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16
```

## Princeton Clusters file system.
Right now you should be in your perseus home folder. While you may find it useful to have a few files here, most time you will save your data on a large filesystem that is backed up to tape & accessible from all of the Princeton supercomputers, Tigress.

Go to your tigress folder. This space is your main folder where we will keep all our simulations.

```bash
cd /tigress/netID
```

You may find it useful to create a symlink in your home directory to your tigress folder. e.g.
```bash
ln -sf /tigress/netID ~/tig
```
This enables you to `cd` to the directory by simply typing `cd ~/tig`.

There also is a scratch directory which is not backed up to tape that you can use to save development runs. Make a  symlink to the scratch disk.
```bash
ln -sf /scratch/gpfs/netID ~/scratch
```
