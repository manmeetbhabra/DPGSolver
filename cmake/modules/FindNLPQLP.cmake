# - Try to find NLPQLP libraries needed
# Input variables:
#  NLPQLP_LIBRARY       - The NLPQLP library directory
#  GFORTRAN_LIBRARY     - The GFORTRAN library directory
# Output variables:
#  NLPQLP_FOUND         - System has MKL
#  NLPQLP_LIBRARIES     - The libraries needed to use NLPQLP

if(NOT NLPQLP_FOUND)
  message("Locating NLPQLP Libraries")
  include(FindPackageHandleStandardArgs)

  # If the user did not specify an MKL root directory, check for the LIBGFORTRAN_DIR and LIBNLPQLP_DIR environment variable 
  if(NOT NLPQLP_LIBRARY OR NOT DEFINED NLPQLP_LIBRARY)
    if(EXISTS $ENV{LIBNLPQLP_DIR})
      set(NLPQLP_LIBRARY "$ENV{LIBNLPQLP_DIR}")
    else()
      message("NLPQLP_LIBRARY could not be set")
    endif()
  endif()

  if(NOT GFORTRAN_LIBRARY OR NOT DEFINED GFORTRAN_LIBRARY)
    if(EXISTS $ENV{LIBGFORTRAN_DIR})
      set(GFORTRAN_LIBRARY "$ENV{LIBGFORTRAN_DIR}")
    else()
      message("GFORTRAN_LIBRARY could not be set")
    endif()
  endif()

  
  # Search for NLPQLP libraries
  find_library(LIBNLPQLP nlpqlp ${NLPQLP_LIBRARY})
  find_library(LIBGFORTRAN gfortran ${GFORTRAN_LIBRARY})

  # Compile libraries together
  set(NLPQLP_LIBRARIES ${LIBNLPQLP} ${LIBGFORTRAN})

  set(NLPQLP_FOUND TRUE)
  message("NLPQLP Libraries set")

endif()
