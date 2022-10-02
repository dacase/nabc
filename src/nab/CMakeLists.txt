
project(NAB C)

#full path to the C file built by rule_parse
set(CHECKEXPR_FILE ${CMAKE_CURRENT_BINARY_DIR}/checkexpr.c)

set(NAB2C_SOURCES
	nab2c.c	
	cgen.c	
	dumpnode.c
	errormsg.c	
	fixexpr.c	
	node.c	
	symbol.c	
	y.tab.c)

#NOTE: to maintain compatibility, the libsym sources in nss are added to libnab here.

set(LIBNAB_C_SOURCES
	axis2frame.c	
	bonded_atoms.c 
	boundsutil.c	
	chirvol.c	
	database.c  
	deform_ener.c 
	dump.c		
	embed.c		
	frameutil.c	
	getdist.c	
	gsub.c	
	hashutil.c	
	jacobi.c	
	match.c		
	matop.c		
	molio.c		
	molsurf.c   
	molutil.c	
	pointutil.c	
	prm.c 
	regexp.c	
	reslib.c	
	select_atoms.c	
    time.c
	traceback.c 
	mme_init.c 
	stringutil.c	
	temputil.c	
	tm_malloc.c
	../nss/tss.c)
	
set(LIBNAB_NAB_SOURCES
	bdna.nab
	dg_helix.nab 
	fd_helix.nab
	getseq_from_pdb.nab
	getxyz_from_pdb.nab
	getchivol.nab
	getpdb_prm.nab 
	linkprot.nab	
	link_na.nab	
	na_anal.nab 	
	setboundsfromdb.nab
	setchiplane.nab
	wc_basepair.nab	
	wc_complement.nab	
	wc_helix.nab 
	../nss/symop.nab)
	
#headers for NAB to install
set(NAB_INSTALL_HEADERS
	nabcode.h 
	nabtypes.h
	defreal.h 
	nab.h 
	xmin_opt.h
	lmod_opt.h
	nab_netcdf.h)

#source files compiled with optimizations (all others are unoptimized)
set(NAB_OPTIMIZED_SOURCES molio.c me_init.c boundsutil.c embed.c bonded_atoms.c ms.c sasad.c)
	
#Rule files in semantics for checkexpr.c
set(RULES SYM_ALLOCATE.rules 
	SYM_AND.rules 
	SYM_ASSERT.rules 
	SYM_ASSIGN.rules 
	SYM_ATSIGN.rules 
	SYM_CALL.rules 
	SYM_DEALLOCATE.rules 
	SYM_DELETE.rules 
	SYM_DONT_MATCH.rules 
	SYM_EQUAL.rules 
	SYM_FOREACH.rules 
	SYM_GREATER.rules 
	SYM_GREATER_EQUAL.rules 
	SYM_IN.rules 
	SYM_INDEX.rules 
	SYM_LBRACK.rules 
	SYM_LESS.rules 
	SYM_LESS_EQUAL.rules 
	SYM_LPAREN.rules 
	SYM_MATCH.rules 
	SYM_MINUS.rules 
	SYM_MINUS_ASSIGN.rules 
	SYM_MINUS_MINUS.rules 
	SYM_MODULUS.rules 
	SYM_MODULUS_ASSIGN.rules 
	SYM_NEGATE.rules 
	SYM_NOT.rules 
	SYM_NOT_EQUAL.rules 
	SYM_OR.rules 
	SYM_PARM.rules 
	SYM_PLUS.rules 
	SYM_PLUS_ASSIGN.rules 
	SYM_PLUS_PLUS.rules 
	SYM_SLASH.rules 
	SYM_SLASH_ASSIGN.rules 
	SYM_STAR.rules 
	SYM_STAR_ASSIGN.rules 
	SYM_TEST.rules 
	SYM_UPARROW.rules 
	SYM_UPARROW_ASSIGN.rules)

#attribute.tab file for rule_parse
set(ATTRIBUTE_TAB ${CMAKE_SOURCE_DIR}/AmberTools/src/nab/semantics/attribute.tab)

#parser files
set(NAB_YACC nabgrm.y)
set(NAB_COMPILED_YACC_C ${CMAKE_CURRENT_BINARY_DIR}/y.tab.c)
set(NAB_COMPILED_YACC_H ${CMAKE_CURRENT_BINARY_DIR}/y.tab.h)

set(NAB_LEX_SOURCE nablex.l)
set(NAB_LEX_LEXED ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.c)

#set up directory for intermediates
set(INTERMED_DIR ${CMAKE_CURRENT_BINARY_DIR}/lex-intermediates)
file(MAKE_DIRECTORY ${INTERMED_DIR})

set(DG_OPTIONS_SOURCE dg_options.l)
set(DG_OPTIONS_INTERMEDIATE ${INTERMED_DIR}/lex.dgo.c)
set(DG_OPTIONS_LEXED ${CMAKE_CURRENT_BINARY_DIR}/lex.dgo.c)

#default compile flags
set_property(DIRECTORY PROPERTY COMPILE_OPTIONS ${NO_OPT_CFLAGS})

#set specific source file properties
set_property(SOURCE ${NAB_OPTIMIZED_SOURCES} PROPERTY COMPILE_FLAGS ${OPT_CFLAGS_SPC})
set_property(SOURCE mme_init.c embed.c PROPERTY COMPILE_DEFINITIONS flex)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#make the lex/yacc parser (needed by checkexpr)
flex_target(NabLex ${NAB_LEX_SOURCE} ${NAB_LEX_LEXED})
flex_target(DgOptionsLex ${DG_OPTIONS_SOURCE} ${DG_OPTIONS_INTERMEDIATE} COMPILE_FLAGS -Pdgo)
bison_target(NabYacc ${NAB_YACC} ${NAB_COMPILED_YACC_C})
add_flex_bison_dependency(NabLex NabYacc)

#do the textual replaces on dg_options
buildtime_file_replace(${DG_OPTIONS_INTERMEDIATE} ${DG_OPTIONS_LEXED} TO_REPLACE stdout REPLACEMENT nabout)

add_custom_target(nab_lex DEPENDS ${NAB_LEX_LEXED} ${DG_OPTIONS_LEXED})


# --------------------------------------------------------------------
# build rule_parse
if(NOT USE_HOST_TOOLS)
	add_subdirectory(semantics)
endif()

# build checkexpr.c using rule_parse
# --------------------------------------------------------------------

# ensure that we have a place to put the HTML files that are generated from rule_parse
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/checkexpr-html)

# convert rules to full paths for the dependency list
set(RULES_FULLPATH "")
foreach(RULE ${RULES})
	list(APPEND RULES_FULLPATH "${CMAKE_CURRENT_SOURCE_DIR}/semantics/${RULE}")
endforeach()

add_custom_command(OUTPUT ${CHECKEXPR_FILE} 
	COMMAND ${RUNNABLE_rule_parse} ARGS 
		--amber-data-dir "${CMAKE_SOURCE_DIR}/dat" 
		--y-tab-h "${NAB_COMPILED_YACC_H}"
		--attributes-tab "${ATTRIBUTE_TAB}"
 		--html-output-dir "${CMAKE_CURRENT_BINARY_DIR}/checkexpr-html"
 		--nab-h-path "${CMAKE_CURRENT_SOURCE_DIR}/nab.h" 
 		-o "${CHECKEXPR_FILE}"
 		${RULES} 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/semantics 
	DEPENDS ${NAB_COMPILED_YACC_H} ${ATTRIBUTE_TAB} ${RUNNABLE_rule_parse} "${CMAKE_CURRENT_SOURCE_DIR}/nab.h"
	COMMENT "Parsing nab rules and generating checkexpr.c")
	
add_library(nab_checkexpr_obj OBJECT ${CHECKEXPR_FILE})
	

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#build nab2c

# flex-generated files include unistd.h unless you tell them not to
if(TARGET_WINDOWS)
	add_definitions(-DYY_NO_UNISTD_H)
endif()

add_executable(nab2c ${NAB2C_SOURCES} $<TARGET_OBJECTS:nab_checkexpr_obj>)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

#the lexed file is included into another source file, so make sure it gets built first
add_dependencies(nab2c nab_lex)
install(TARGETS nab2c DESTINATION ${BINDIR})

# also build the MPI version.
# The MPI version doesn't actually use MPI; it just generates code that uses MPI, so we just leave it enabled all the time.
copy_target(nab2c mpinab2c INSTALL)
target_compile_definitions(mpinab2c PRIVATE MPI)
add_dependencies(mpinab2c nab_lex)

install(TARGETS mpinab2c DESTINATION ${BINDIR})


if(NOT BUILD_HOST_TOOLS)
	# only build supporting libraries if not building host tools

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#libnab
	set(NAB_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})
	nab_compile(LIBNAB_COMPILED_FILES ${LIBNAB_NAB_SOURCES})
	add_library(libnab STATIC ${LIBNAB_COMPILED_FILES} ${LIBNAB_C_SOURCES})
	remove_prefix(libnab)
	make_pic_if_needed(libnab)
	
	#add interface link libraries
	target_link_libraries(libnab cifparse sff libpbsa netlib)
	
	target_include_directories(libnab PUBLIC .)
	
	install_libraries(libnab)
	
	# MPI version
	if(MPI)
		nab_compile(LIBNABMPI_COMPILED_FILES ${LIBNAB_NAB_SOURCES} MPI)
		add_library(libnab_mpi STATIC ${LIBNABMPI_COMPILED_FILES} ${LIBNAB_C_SOURCES})
		remove_prefix(libnab_mpi)
		make_pic_if_needed(libnab_mpi)
		
		#add interface link libraries
		target_link_libraries(libnab_mpi cifparse sff_mpi libpbsa netlib)
		
		target_include_directories(libnab_mpi PUBLIC .)
		
		install_libraries(libnab_mpi)
	endif()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#header files
	
	install(FILES ${NAB_INSTALL_HEADERS} DESTINATION ${INCDIR})
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# finally, build the nab wrapper program

if(NOT BUILD_HOST_TOOLS)
	
	include(${CMAKE_CURRENT_LIST_DIR}/NABCompileSettings.cmake)

	#figure out parameters that get compiled into the nab wrapper
	set(REQUIRED_HEADERS "string")
	
	create_list_source_file(TARGETNAME nab 
		OUTPUTFILEVAR NAB_VECTOR_FILE 
		INCLUDE "string" 
		QUOTE 
		builtinLibraries std::string NAB_LINK_LIBRARIES 
		builtinCompilerFlags std::string NAB_COMPILE_FLAGS
		builtinLinkDirs std::string NAB_LINK_DIRECTORIES)
		
	add_executable(nab nab.cpp "${NAB_VECTOR_FILE}")

	target_compile_definitions(nab PRIVATE CC="${NAB_C_COMPILER}" LD="${NAB_LINKER}")
	 	
	install_executables_wrapped(TARGETS nab)
	
	if(MPI) # make MPI version
		
		create_list_source_file(TARGETNAME mpinab 
			OUTPUTFILEVAR MPINAB_VECTOR_FILE 
			INCLUDE "string" 
			QUOTE 
			builtinLibraries std::string MPINAB_LINK_LIBRARIES 
			builtinCompilerFlags std::string MPINAB_COMPILE_FLAGS
			builtinLinkDirs std::string MPINAB_LINK_DIRECTORIES)
		
		add_executable(mpinab nab.cpp "${MPINAB_VECTOR_FILE}")
		target_compile_definitions(mpinab PRIVATE MPI CC="${MPINAB_C_COMPILER}" LD="${MPINAB_LINKER}")
		install_executables_wrapped(TARGETS mpinab COMPONENT MPI)
	endif()
		
endif()
