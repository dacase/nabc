# Script which determines the correct compile flags, link flags, link libraries, and link directories for the NAB compiler wrapper to use.

# currently, we link with the fortran compiler since we need to have libgfortran or the equivalent runtime linked in
set(NAB_C_COMPILER ${CMAKE_C_COMPILER})
set(NAB_LINKER ${CMAKE_Fortran_COMPILER})

# Compile Flags
# --------------------------------------------------------------------

# Just use the bog-standard Amber CFLAGS, but disable warnings because NAB-generated code contains lots of syntax warnings
separate_arguments(NAB_COMPILE_FLAGS UNIX_COMMAND ${CMAKE_C_FLAGS})
list(APPEND NAB_COMPILE_FLAGS ${OPT_CFLAGS} -w)

# NOTE: nab.cpp adds "-I$AMBERHOME/include"

# Link Libraries & Directories
# --------------------------------------------------------------------

# prevent linking of FORTRAN main
set(NAB_LINK_OPTIONS "")
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
	list(APPEND NAB_LINK_OPTIONS -nofor-main)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
	list(APPEND NAB_LINK_OPTIONS -Mnomain)
endif()

# Amber libraries that NAB needs to link to programs it builds
set(NAB_BUILTIN_LIBRARIES libnab sff libpbsa rism cifparse fftw netcdf)

resolve_cmake_library_list(NAB_LIBS_PATHS ${NAB_BUILTIN_LIBRARIES} ${NETLIB_LIBRARIES})
resolved_lib_list_to_link_line(NAB_LINK_LIBRARIES NAB_LINK_DIRECTORIES ${NAB_LIBS_PATHS})
list(APPEND NAB_LINK_LIBRARIES ${NAB_LINK_OPTIONS})

#printvar(NAB_COMPILE_FLAGS)
#printvar(NAB_LINK_LIBRARIES)
#printvar(NAB_LINK_DIRECTORIES)

if(MPI)

	set(MPINAB_COMPILE_FLAGS "")	
	set(MPINAB_LINK_LIBRARIES )
	set(MPINAB_LINK_DIRECTORIES "")
		
	if(EXISTS "${MPI_C_COMPILER}" AND EXISTS "${MPI_Fortran_COMPILER}")
	
		# if we have compiler wrappers available, just use that.
		set(MPINAB_C_COMPILER ${MPI_C_COMPILER})
		set(MPINAB_LINKER ${MPI_Fortran_COMPILER})
		
		# link Amber libraries, but not MPI itself
		set(MPINAB_BUILTIN_LIBRARIES nab_mpi sff_mpi pbsa rism_mpi cifparse fftw fftw_mpi netcdf ${NETLIB_LIBRARIES})
	else()
		
		# no wrappers, so try to populate the nab wrapper with the MPI flags we are using
		set(MPINAB_C_COMPILER ${CMAKE_C_COMPILER})
		set(MPINAB_LINKER ${CMAKE_Fortran_COMPILER})
		
		# MPI Compile Flags
		# --------------------------------------------------------------------
	
		set(MPINAB_COMPILE_FLAGS ${MPI_${LANG}_COMPILE_FLAGS})
	
		# MPI Link Libraries & Directories
		# --------------------------------------------------------------------
		
		set(MPINAB_BUILTIN_LIBRARIES nab_mpi sff_mpi pbsa rism_mpi cifparse fftw fftw_mpi netcdf ${NETLIB_LIBRARIES} mpi_c)
	endif()
	
	resolve_cmake_library_list(MPINAB_LIBS_PATHS ${MPINAB_BUILTIN_LIBRARIES})
	resolved_lib_list_to_link_line(MPINAB_LINK_LIBRARIES MPINAB_LINK_DIRECTORIES ${MPINAB_LIBS_PATHS})
	list(APPEND MPINAB_LINK_LIBRARIES ${NAB_LINK_OPTIONS})

			
	#printvar(MPINAB_COMPILE_FLAGS)
	#printvar(MPINAB_LINK_LIBRARIES)
	#printvar(MPINAB_LINK_DIRECTORIES)
endif()