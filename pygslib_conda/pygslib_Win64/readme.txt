This package requires manual fix to make it work
a) add custom channel 
conda config --add channels https://conda.binstar.org/ioos
b) run build
conda build pygslib_W64
c) if numpy <1.10 installed you may need to fix a bug in distutils
go to the build environment and around the line 337 you will see:

  File "C:\Users\yourusername\Anaconda\envs\_build\lib\site-packages\numpy\distutils\fcompiler\gnu.py", line 337

if is_win64():
    c_compiler = self.c_compiler
    if c_compiler and c_compiler.compiler_type == "msvc":
        return []
    else:
        raise NotImplementedError("Only MS compiler supported with gfortran on win64")
		
rewrite the code like this:

if is_win64():
        c_compiler = self.c_compiler
        if c_compiler and c_compiler.compiler_type == "msvc":
            return []
        else:
            return [] #raise NotImplementedError("Only MS compiler supported with gfortran on win64")
			
then run again 

conda build pygslib_W64

#issue with MSVC compiler
This may not work, instead use MINGW 
conda install mingw
conda install libpython  # this may prevent Anaconda to use MINGW instead MSVC