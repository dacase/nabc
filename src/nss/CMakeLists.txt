project(nss)

# NOTE: the original makefile compiles tss.c and symop.nab in this folder, then adds them with ar to the preexisting libnab archive.
# So they're locically nab source, but they're in this folder and Makefile.
# I'm not exactly sure why.
# I just moved them into the libnab source list in the nab CMakeLists.txt.

set(NAB_EXECUTABLES matextract matgen matmerge matmul transform tss_init tss_main tss_next)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#not sure if this stuff is supposed to be optimized.  I'm guessing yes.
add_compile_options(${OPT_CFLAGS})

#process one-file nab executables
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
foreach(NAB_EXEC ${NAB_EXECUTABLES})
    add_nab_executable(${NAB_EXEC} NAB_SOURCES ${NAB_EXEC}.nab)
endforeach()

install(TARGETS ${NAB_EXECUTABLES} DESTINATION ${BINDIR})