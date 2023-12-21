efros_leung - Efros and Leung texture synthesis

===========================

Version 1.0 - Janvier 11, 2013
by Cecilia Aguerrebere <aguerreb@telecom-paristech.fr>

-------------------------------------------------------------------------------
Introduction
------------

efros_leung is an implementation of the texture synthesis algorithm by 
Efros and Leung described in the paper:


  "Texture synthesis by non-parametric sampling"
  by Alexei A. Efros and Thomas Leung, Proc. Int. Conf. on 
  Comp. Vision, vol 2, pages 1033–1038, Kerkyra, Greece, 1999.

-------------------------------------------------------------------------------
Files
-----

README.txt          - This file.

-- src folder:
Makefile            - Compilation instructions for 'make'.
efros_leung_cmd.c   - command line interface for efros_leung, ANSI C code.
efros_leung.c  	    - efros_leung module for color images ANSI C code
synthesis.c	    - basic synthesis sub-functions ANSI C code
auxiliary.c         - auxiliary functions ANSI C code
io_png.c  	    - png input/ouput functions ANSI C code
mt.c  		    - random number generation

-- inc folder:
efros_leung.h       - efros_leung module ANSI C header	
synthesis.h	    - basic synthesis module ANSI C header	
auxiliary.h         - auxiliary module ANSI C header	
io_png.h            - png input/output module ANSI C header	
mt.h                - random number generator module ANSI C header	

-- doc folder:
Html code documentation. 

-------------------------------------------------------------------------------
Compiling
---------

Requirements:

- libpng - Portable Network Graphics (PNG) Reference Library (tested with version 1.2.46)
- OpenMP (can be disabled removing the -fopenmp flag on the Makefile file) 

To compile run

  make

from the directory where the source codes and the Makefile are located.

To verify a correct compilation you can apply efros_leung to the test
image "texture.png" and compare the result to the provided one "output_example.png".
For that purpouse, from the src folder, run the command

./efros_leung texture.png

The result is saved in the "output.png" file.

-------------------------------------------------------------------------------
Running efros_leung command
-------------------

Run ./efros_leung without parameters to see the program usage.

Usage: ./efros_leung [options] <input file> <output file> <synthesis map sample file> <synthesis map output file>

Only  PNG 8 bits images are supported.

Options:
  -p <number>  Half patch size  (default 4). Patches are of width 2p+1.
  -t <number>  Tolerance parameter (default 0.1).
  -s <number>  Output image size (default 128) 

Usage example:

./efros_leung -p 5 texture.png 

This example uses the default parameters and a patch of size 11x11. The output is saved 
in the "output.png" file. The synthesis maps are saved as "map.png" and "copy_map.png".

-------------------------------------------------------------------------------
Code Documentation
------------------

There is a HTML documentation of the code on the directory 'doc'. The
entry point is the file 'doc/index.html' that should be opened with a
web browser. The documentation was automatically generated from the
source code files using the Doxygen documentation system, see
http://www.stack.nl/~dimitri/doxygen/.


-------------------------------------------------------------------------------
Copyright and License
---------------------

Copyright (c)

efros_leung is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

efros_leung is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------
Feedback
------

Do not hesitate to contact me if you have any comments, especially about errors,
bugs, or strange results.
