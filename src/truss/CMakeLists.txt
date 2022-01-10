set(programs
    Example01
    Example02
    Example03
    Example04)
		
# include OpenMP
find_package(OpenMP)

foreach(_program ${programs})
	add_executable(${_program} ${_program}.cc)
	target_link_libraries(${_program} PRIVATE OpenMP::OpenMP_CXX)
endforeach()
