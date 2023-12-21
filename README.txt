Efros and Leung texture synthesis

===========================

Version 1.0 - Janvier 11, 2013
by Cecilia Aguerrebere <aguerreb@telecom-paristech.fr>


Introduction
------------

The file efros_leung_1.0.zip contains the source code of the implementation of the Efros and Leung texture synthesis algorithm described in the IPOL article of which this file is part:

  "Exemplar-based texture synthesis : the Efros-Leung algorithm" 
  by Cecilia Aguerrebere, Yann Gousseau and Guillaume Tartavel,
  Image Processing On Line, 2013

Two versions of the texture synthesis algorithm by Efros and Leung are provided in this file. The first one follows exactly the algorithm described in the article 

  "Texture synthesis by non-parametric sampling"
  by Alexei A. Efros and Thomas Leung, Proc. Int. Conf. on 
  Comp. Vision, vol 2, pages 1033–1038, Kerkyra, Greece, 1999.

The source code of this implementation is found in the sub-folder "efros_leung_orig". Please refer to the README file in that folder for details on the content, compilation and utilisation of that version.

The second one is an accelerated version of the texture synthesis algorithm by Efros and Leung. The details of this implementation are described in the IPOL article previously mentioned.

The source code of this implementation is found in the sub-folder "efros_leung_accel". Please refer to the README file in that folder for details on the content, compilation and utilisation of that version.

Content
------------

efros_leung_orig   - Folder containing the classical implementation of the Efros and Leung texture synthesis algorithm.

efros_leung_accel  - Folder containing the accelerated implementation of the Efros and Leung texture synthesis algorithm.



