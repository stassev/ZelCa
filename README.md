    Author: Svetlin V. Tassev (Princeton U)
    Initial public release date: Nov 25, 2013

    This is the Zel'dovich Calculator (ZelCa): A code calculating the 
     real space and redshift space matter 2-pt functions, as well as the
     3-point function in the Zel'dovich Approximation. The code is 
     based on the paper "N-point Statistics of Large-Scale Structure in 
     the Zel'dovich Approximation", S. Tassev (2013).
     Check that paper for the details.
     
    If you use the ZelCa for scientific work, I kindly ask you 
    to reference the paper above.
    
    
    * ZelCa is open-source software, distributed under the GPLv3 license.
    * To compile the code, one needs the GSL library
      (http://www.gnu.org/software/gsl/).
    * It also needs Cuba (http://www.feynarts.de/cuba/)  
      and Eigen (http://eigen.tuxfamily.org), which in turn
      depends on cmake (http://www.cmake.org/). Those are usually not 
      readily available. So, if any of those three are missing from your
      system, ZelCa will try to download and build them from source
      for you.
    * ZelCa comes with a building script which tries to take care of the
      above dependencies and sets some environment variables. To build
      ZelCa initially one needs to run
        ./build.sh initial-zelca
      and it will guide the user through the process of building  
      the needed dependencies if missing. Just follow the instructions.
      If dependencies have to be built, they will appear 
      in "./DEPENDENCIES/". Do not delete that folder if you plan to 
      rebuild ZelCa at a later time. Otherwise you'll have to rebuild
      the dependencies. Also, the initial ZelCa build generates a file
      "install.env". Messing with it may result in the need to re-issue
      the above command.
    * Once built, the executable will be in "./executable/ZelCa"
    * For all subsequent builds of ZelCa, one should use:
        ./build.sh zelca
    * Type "./build.sh" for further options.
    * Example input and output files are provided in "./executable/"
    * For the Zel'dovich calculations, one needs to precompute Chi and
      Gamma. The code offers to do that once started. ZelCa also needs 
      a matter transfer function. It assumes it is generated with
      CAMB (http://camb.info/).
    * After each calculation (precopmuting Chi and Gamma; calculating
      linear theory; calculating 2pt in ZA; etc...), the code will exit.
      If you want a different calculation than the chosen one, rerun the 
      code. 
    * To understand the output, the easiest way is through example. 
      The directory "./executable/SAGE_analsysis/" contains a 
      Sage (http://www.sagemath.org/) worksheet (*.sws) showing you how
      how to process the data. If you do not have Sage at hand and do
      not want to bother installing it, a printout of the worksheet is 
      available in the folder as a *.pdf. It shows how I generated the 
      plots for the paper cited above, as well as for "Lagrangian or 
      Eulerian...", S. Tassev (2013).
      
      
    
    Warnings and more:
    * The growth factor and rate of growth calculations assume spatially
      flat cosmology.
    * Note that I never check whether input files exist, so if you get 
      segfaults, check the filenames you input.
    * The linear calculations were never intended to be fast, and the 
      code for them is a bit of a mess. The pieces of the code doing the 
      Zel'dovich approximation (ZA), however, are fairly clean and 
      organized. So, should be straightfoward to understand and modify
      as needed. The cleanliness of the ZA pieces of the code are 
      entirely thanks to the Eigen library, which is also blazingly
      fast at numerical linear algebra. I tested with Blas/Lapack, as 
      well as with building an optimized Atlas library. None of them 
      worked nearly as fast.
    * As a rule of thumb, the code will suffer severely from roundoff 
      errors, when the quantity one needs to calculate is <<<1, as one 
      has to subtract the disconnected pieces, one of them being 
      always 1. 
    * The code will also suffer when one has to sample large regions of 
      Q space for convergence as volume factors will explode the 
      integrand, leading once again to roundoff issues. The code tries 
      to minimize this effect, recasting the disconnected pieces and 
      pushing those inside the integrand of the n-pt, thus cancelling a
      lot of contributions from |Q|>>|X|. This is also described in the 
      Numerical implementation section of the paper. That trick helps, 
      but only to some extent.
    * As an example, the code will suffer from roundoff issues at 
      redshifts (z>~4) for the 2-pts. This particular problem has a 
      somewhat trivial solution: expand the ZA 2-pt, code it up, and 
      do it again. But that is beyond ZelCa.
    * The 3pt piece was tested only for certain triangle configurations. 
    * When one encounters convergence issues, better read the Cuba 
      documentation, and see whether switching integrators and changing 
      integrator parameters will help. Cuba integrates over the unit 
      hypercube. Thus choosing a good map between the hypercube and 
      Q-space is crucial. Unfortunately, what works is more
      of an art, but the "partview" utility of Cuba can be of great
      help. Read the Cuba documentation about it, and note that one 
      may be missing "partview" from their Cuba build if they do not 
      have the qmake (with a "q"). If you use the build script coming 
      with this code to build Cuba and you have qmake, but you are still
      missing partview (and you want to use it), try commenting out
      the line:
        sed -e "s/qmake/qmake-qt4/g" -i makefile.in
      in ./scripts/cuba-build.sh
    * One can do Cuba debugging by setting the Cuba verbosity level 
      between 0 and 3 by issuing the following before running ZelCa 
      (for verbosity=3):
        export CUBAVERBOSE=3
    * One can change the number of threads that Cuba uses by issuing the
      following before running ZelCa (for 8 cores):
        export CUBACORES=8
      

