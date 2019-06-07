## First time accessing Perseus and/or Tigressdata
The first time you try to access Perseus and/or Tigressdata from you local machine you must be on Princeton's network or use a VPN connection. If you only ever plan on accessing Princeton's cluster from campus, you do not need to follow these instructions. But I highly recommend that you use SSH keys to log in regularly and set-up a VPN just in case.

**Note:** SSH keys are a much easier way to access the cluster. If you do not follow the next steps, then you will ALWAYS have to use VPN to connect to the cluster. In any rate, SSH key access will be going away soon with two-factor authentication being implemented. :(

## Setting up a VPN
If you're off-campus first you'll need to install a VPN to connect to the cluster. VPN access to Princeton is well documented and you will find information about how to do it [here.](https://princeton.service-now.com/snap?sys_id=6023&id=kb_article). Follow the instructions on how to install sonic wall connect and VPN into the local Princeton network.
 
## SSH keys
To avoid using VPN or wired connection anytime you want to access the cluster, you will need to set up a pair of ssh keys that will be shared by your local machine and the remote server. The pair consists of a *private key* (default name: id_rsa) and a *public key* (default name: id_rsa.pub). **Never share your private key, treat it like a password.**  

For more info about ssh keys go [here](https://wiki.archlinux.org/index.php/SSH_keys)

## Creating SSH keys on your local computer.

Open a terminal on your local computer and see if you already have an ssh key in your local folder by typing
```bash
ls ~/.ssh
```

If your directory is empty, you must create a pair of keys following [these instructions](#if-you-do-not-have-any-ssh-keys-in-your-local-folder). If you have a pair of keys in the directory, [go here](#if-you-do-have-an-ssh-key-in-your-local-folder) 


### If you do not have any ssh keys in your local folder

You must generate a pair of SSH keys. To generate RSA keys, on the command line, enter:  
```bash
ssh-keygen -t rsa
```
You will receive a prompt to supply a filename (for saving the key pair) and a password (for protecting your private key. that looks something like this:
```
$ ssh-keygen -t rsa
Generating public/private rsa key pair.
Enter file in which to save the key (/home/USER/.ssh/id_rsa):
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /home/USER/.ssh/id_rsa.
Your public key has been saved in /home/USER/.ssh/id_rsa.pub.
The key fingerprint is:
ef:69:3b:9e:3b:2d:99:0d:ac:57:4e:b2:92:82:bd:9f dave@hostname
The key's randomart image is:
+--[ RSA 2048]----+
|                 |
|                 |
|                 |
|                 |
|        S.       |
|         .+ o    |
|     o   o.%     |
|    . o +oXo+    |
|      .+E=B*     |
+-----------------+
``` 

**Filename:** To accept the default filename (and location) for your key pair, press `Enter` or 
`Return` without entering a filename. Alternatively, you can enter a filename (e.g., `my_ssh_key`) at the prompt, and then press `Enter` or `Return`. However, many remote hosts  are configured to accept private keys with the default filename and path (`~/.ssh/id_rsa` for RSA keys; `~/.ssh/id_dsa` for DSA keys) by default. Consequently, to authenticate with a private key that has a different filename, or one that is not stored in the default location, you must explicitly invoke it either on the SSH command line or in an SSH client configuration file (`~/.ssh/config`). It's easiest to just accept the default filename (More [here](https://kb.iu.edu/d/aews)).

**Password:** Enter a password that contains at least five characters, and then press `Enter` or 
`Return`. If you press `Enter` or `Return` without entering a password, your private key will be generated without password-protection 

### If you do have an ssh key in your local folder 
You can choose to either generate a different pair of keys for the connection to the remote server (e.g. `id_rsa_perseus.pub`) following the instructions in the [previous section](#if-you-do-not-have-any-ssh-keys-in-your-local-folder). Otherwise, or you can use the already existing key (e.g. `id_rsa.pub`), and you can move on to the next step.

## Putting your SSH keys on the cluster

After generating the pair of keys, you must add the contents of your public key to the authorized_keys file in the remote server. Use SCP to securely copy your public key (e.g., `~/.ssh/id_rsa.pub`) onto Perseus and/or Tigressdata2. On the local network or using VPN, securely:
```bash
scp ~/.ssh/id_rsa.pub netID@perseus.princeton.edu
```
You'll be prompted for your account password. Your public key will be copied to your home directory (and saved with the same filename) on the remote system.

Log into the remote system using your account username and password.
```bash
ssh netID@perseus.princeton.edu
```
On the command line in perseus, enter the following commands that ensures that persues has a `~/.ssh/authorized_keys` file:
```bash
   mkdir -p ~/.ssh 
   touch ~/.ssh/authorized_keys
```
*Note:*If your account on Perseus already has a `~/.ssh/authorized_keys` file, executing these commands will not damage the existing directory or file.

On the remote system, add the contents of your public key file (e.g., `~/id_rsa.pub`) to a new line in your `~/.ssh/authorized_keys` file; on the command line in Perseus, enter:
```bash
cat ~/id_rsa.pub >> ~/.ssh/authorized_keys
```

You may want to check the contents of ~/.ssh/authorized_keys to make sure your public key was added properly; on the command line, enter:
```bash
less ~/.ssh/authorized_keys
```
*Note:* you can close less by hitting 'q'.

You may now safely delete the public key file (e.g., `~/id_rsa.pub`) from your account on the remote system; on the command line, enter:
```bash
rm ~/id_rsa.pub
```
Alternatively, if you prefer to keep a copy of your public key on the remote system, move it to your .ssh directory; on the command line, enter:
```bash
mv ~/id_rsa.pub ~/.ssh/
```

*You should follow these instructions both for tigressdata2 & perseus.*

## Logging in using your ssh key

If you have only one pair of keys on your local machine then you can then access the remote server from any network by typing:
```bash 
ssh netID@perseus.princeton.edu
```

If you have more than one pair of keys on your local machine, then you need to specify which key to be used for accessing the cluster. For example, let’s assume you have two keys in your local machine (`id_rsa` and `id_rsa_new`) and you have copied the contents of `id_rsa_new.pub` into the `~/.ssh/authorized_keys` file on Perseus. Then, you can access the cluster by typing: 
```bash
ssh -i ~/.ssh/id_rsa_new  netID@perseus.princeton.edu
```

## Notes

1. If you add manually the content of the public key to the file authorized_keys, be careful not to copy paste the content of the key as displayed on the terminal to the file. This may create broken lines and the key will not be recognized! Use `cat` as we showed above. For more details about the authorized_keys file see [here](https://wiki.mcs.anl.gov/IT/index.php/SSH_Keys:authorized_keys)

2. If you want to avoid typing long commands on the terminal you can either create an alias 
to your .bashrc file on your own computer, e.g.,

```
alias sshpers='ssh -i $HOME/.ssh/id_rsa_new netID@perseus.princeton.edu
```

Or you canmodify your ssh_config file. This can be found in `~/.ssh/config` or `/etc/ssh/ssh_config` folders. This is possible for your local machine, but usually one cannot modify the ssh_config file of a remote server, so using an alias is likely the better solution. For more info 
check [here](https://kb.iu.edu/d/aews)

3. If you get the following error:
``` Permission denied (publickey,gssapi-keyex,gssapi-with-mic)```

you should read the log by typing
```
ssh -vvv -i  ~/.ssh/id_rsa_new username@perseus.princeton.edu
```
Where you should change ```~/.ssh/id_rsa_new``` to the location of your private ssh key.

Useful information can be found at:
* https://askubuntu.com/questions/692480/ssh-authentication-with-id-rsa-key-not-working
* https://www.experts-exchange.com/questions/27464003/SSH-Pub-Key-Exchange-Failure.html