#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#define MACSTYLE

#define CONCAT(prefix, name) prefix ## name
#ifdef MACSTYLE
	#define F90(name) CONCAT(name,_)
#else
	#define F90(name) name
#endif


/*change umask in unix type systems to something else than the default*/

void F90(umask_f)(void)
{
//	mode_t umasknumber= (S_IRGRP | S_IROTH | S_IRUSR | S_IWUSR);

	umask(18);
}


/* Remove a file */

void F90(remove_f)(const char *path, int *ierr)
{
 if (*ierr = remove(path)) *ierr = errno;
}

/* move a file */

void F90(move_f)(const char *pathprev, const char *pathfinal, int *ierr)
{
 if (*ierr = rename(pathprev, pathfinal)) *ierr = errno;
}


/* Create a directory */

void F90(mkdir_f)(const char *path, int *ierr)
{
 char uppath[256], *p;
    
 if (mkdir(path,S_IRWXU | (S_IRGRP | S_IXGRP ) | (S_IROTH | S_IXOTH) ))
   switch (errno) {
     case ENOENT : // A component of the path does not exist
        
        // get upper path
        
        strcpy(uppath, path);
        p = uppath + strlen(uppath);
        while(*p!='/') p--;
        *p=0;
        
        // recursively build the path
        
        F90(mkdir_f)(uppath, ierr);
        if (!*ierr) F90(mkdir_f)(path, ierr);
     
     case EEXIST : // if directory already exists ignore the error
        *ierr = 0;
        break;
     default: *ierr = errno;
   }
 else *ierr = 0;
 
}

/* change working a directory */

void F90(chdir_f)(const char *path, int *ierr)
{
  if (*ierr = chdir(path)) *ierr=errno;
}

/* get current working a directory */

void F90(getcwd_f)(char *path, int *len, int *ierr)
{
  char *buf = (char*) malloc( (*len+1) );

  if (getcwd(buf, *len)) {
     strcpy(path, buf);
     *ierr = 0;
   } else {
    strcpy(path, "\0");
    *ierr=errno;
   }
  
  free(buf);
}

/* Create a symbolic link */

void F90(symlink_f)(const char *name1, const char *name2, int *ierr)
{
  if (*ierr = symlink(name1,name2)) *ierr=errno;
}

/* Get the environment variable value */

void F90(getenv_f)(const char *name, char *value, int *ierr)
{
 char *lvalue;
  
 lvalue = (char *) getenv(name);
  
 if (lvalue == NULL) *ierr = -1;
 else 
 {
   *ierr = 0;
   strcpy(value,lvalue);
 }
 
}

/* Get the hostname */

void F90(gethostname_f)(char *hostname, int *ierr)
{
 char lhostname[256];
  
 if (gethostname(lhostname, 255)) {
   strcpy(hostname,"\0");
   *ierr = errno; 
 } else {
   strcpy(hostname,lhostname);
   *ierr = 0;
 }
}

