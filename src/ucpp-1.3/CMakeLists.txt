# AMBER CMake buildfile for ucpp
project(ucpp C)

#we are buiding a standalone executable
add_definitions(-DSTAND_ALONE)

add_executable(ucpp mem.c nhash.c cpp.c lexer.c assert.c macro.c eval.c)

install(TARGETS ucpp DESTINATION ${BINDIR})

