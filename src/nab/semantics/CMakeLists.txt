#This file builds rule_parse, a component of the NAB compiler.
#NOTE: it'snot included at all if USE_HOST_TOOLS is enabled.

set(RULEPARSE_SOURCE parse2.c)

add_executable(rule_parse ${RULEPARSE_SOURCE})
set_property(TARGET rule_parse PROPERTY COMPILE_OPTIONS ${NO_OPT_CFLAGS})

#install rule_parse in host tool mode
if(BUILD_HOST_TOOLS)
	install(TARGETS rule_parse DESTINATION ${BINDIR})
endif()