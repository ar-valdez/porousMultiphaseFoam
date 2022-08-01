---
		 _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    

**** PorousMultiphaseFoam (PMF) for OpenFOAM ****

---

# Branches

## The current branches and source codes are being used

- OpenFOAM v9 from openfoam.org
- branch **openfoam-v9**

# Stuff done

- Corey 1954 model implemented
- Chierici 1984 model implemented
- LET 2005 model implemented

# Stuff to do

- Make a simple test following Gabriel's advice

# Working tests

- Make a rectangular domain, taking core params. and run the injection to remove
oil from the reservoir.

# Installation instructions :

- First, source the OpenFOAM configuration file, i.e. (example for ubuntu
  version) :

> source /opt/openfoamv6/etc/bashrc

- then in the "porousMultiphaseFoam" directory, run :

> ./Allwmake -jX

  to install the package (with X the number of processors).

- Dynamic libraries are compiled in the standard OpenFOAM user directory :

> $FOAM_USER_LIBBIN

- The executable solver "impesFoam" is placed in the standard OpenFOAM user
  directory $FOAM_USER_APPBIN.

- Each tutorial directory contains "run" and "clean" files to test installation
  and validate the solver.

- A python script runTutorials.py can be used to test all components.

- To remove compilation and temporary files, run :

> ./Allwclean

- see the ReleaseNotes.txt file for detailed information about the toolbox.

---
