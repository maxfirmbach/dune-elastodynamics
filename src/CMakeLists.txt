set(programs
		Example01
		Example02
		Example03)

foreach(_program ${programs})
	add_executable(${_program} ${_program}.cc)
endforeach()
