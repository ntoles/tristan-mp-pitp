## Tigressdata2

Tigressdata2 is a computer that has access to the files you created when running your simulations. It has a lot of RAM and is very convenient way to analyze your simulations without have to move the data anywhere.  Most of the time, you'll login to tigressdata2 via VNC, a way to have a remote desktop.

## Starting a VNC server on tigressdata2

Open terminal on your local machine and login to tigressdata2 
```bash
ssh netID@tigressdata2.princeton.edu 
```
Start a vncserver on tigressdata2. The geometry flag sets the size of the vnc window.
```bash
vncserver -geometry 1680x1050
```

The first time you start vncserver you will be prompted to create a password to access you vnc desktops. This can be anything you like. Should you choose to change it in the future simply delete the file ~/.vnc/passwd, and the next time to start a vnc session you will be prompted to recreate a password.

You may choose to have a view-only password, where others can view your desktop, but keystrokes are disabled.

You should get an output that looks like this:
```
Desktop 'TurboVNC: tigressdata2.princeton.edu:1 (NetID)' started on 
display tigressdata2.princeton.edu:1

Creating default startup script /home/NetID/.vnc/xstartup.turbovnc
Starting applications specified in /home/NetID/.vnc/xstartup.turbovnc

Log file is /home/NetID/.vnc/tigressdata2.princeton.edu:1.log
```

You will want to note the note the session number. In the above example it is :1.

You can always check what vnc session(s) you have open with the `vncserver -list` command which gives an output like:
```
TurboVNC server sessions:

X DISPLAY # PROCESS ID
:1 18966
```

## Connecting to Tigressdata2
First download and install TurboVNC on your local machine: [turboVNC](http://sourceforge.net/projects/turbovnc)

From your local machine, you will use the TurboVNC viewer to connect, but first you need to set up an ssh tunnel. 
This is necessary because TurboVNC is not encrypted. The session number from earlier is 'n' here. e.g. on a Unix machine like MacOS or Linux type
```bash 
ssh -A -L <5900+n>:localhost:<5900+n> <user>@tigressdata2.princeton.edu
```
Example for the above session (:1)
```bash
ssh -A -L 5901:localhost:5901 netID@tigressdata2.princeton.edu
```
You may see a message or dialog box similar to this: `The server's host key is not cached in the registry.` If you do, type or click `yes`.

Enter the password associated with your netID.

Now on your local computer, open the TurboVNC Viewer and in the VNC server box type: `localhost:<n>` and click connect. For example: `localhost:1`.

Enter the password you set up when you created your vncserver session above.

You should now be able to see your tigressdata2 desktop.

## Closing a VNC Session
You can close the TurboVNC viewer, reconnect at a later time, and pick up right where you left off. Just follow the procedure for connecting to a remote VNC session above.

If the tunnel is closed, the session will end, but your desktop will still be intact on tigressdata. Just reconnect.
 
When you are finished with a VNC session entirely, you should end it.

To end a vnc session
ssh to tigressdata2

If you haven't already, load the turbovnc module
```bash
module load turbovnc
```

Check the number of your vnc session with `vncserver -list` and note the display # of the session you want to end
 
Do `vncserver -kill :<n>` where n is the number (with a colon).

e.g
```bash
vncserver -kill :1
```

## Enabling tigressdata2's GPUs
To take full advantage of the local graphics hardware of tigressdata2 we use VirtualGL. Note this will only work within a VNC session. The instructions here are for a terminal open on your VNC desktop.

Check to make sure you have the virtualgl module loaded, if not do
```bash
module load virtualgl
```

Start a program as you normally would but precede the command with vglrun. Some examples:
```bash
vglrun visit
vglrun paraview 
vglrun /usr/licensed/matlab/bin/matlab
```

## Jupyter notebook through ssh

Most of the time the data is stored on the server, but it's a lot more convenient and sometimes faster to analize the results from a personal laptop. For that you can start a jupyter session on a `tigressdata2` machine and connect to it through ssh from you personal device and access all the jupyter capabilities from your local browser. 

For that you would first need a jupyter loaded on `tigressdata2`, which is contained in `anaconda` package, so run this in terminal

```bash
module load anaconda
```

Then you want to run the following command in any folder you want the jupyter root to be in

```bash
jupyter notebook --no-browser --port=8889 --ip=127.0.0.2
```

Port number and ip could be anything; if a particular value doesn't work, it may be used by some other process. You may also specify a bash function in your `~/.bashrc` or `~/.bash_profile` in the following way

```bash
jupyter_run() {
    jupyter notebook --no-browser --port=$1 --ip=127.0.0.2
}
```
which you later can use the following way `$ jupyter_run 8889`.

Now you need to connect to this server from your local laptop through ssh. For that you can simply run the following command in your local terminal

```bash
ssh -N -L localhost:8889:127.0.0.2:8889 <username>@tigressdata2.princeton.edu
```
where `<username>` is your `tigressdata2` username, and the port number `8889` should correspond to the one specified earlier. In the same fashion you can create an alias function in `~/.bashrc` or `~/.bash_profile`

```bash
ssh_jupyter() {
    ssh -N -L localhost:$1:127.0.0.2:$1 <username>@tigressdata2.princeton.edu
}
```
and call it later as `ssh_jupyter 8889`. 

You can now access the jupyter notebook via `localhost:8889` in your local browser. 

## Additional Information
See instructions [here](https://researchcomputing.princeton.edu/faq/how-do-i-use-vnc-on-tigre)
