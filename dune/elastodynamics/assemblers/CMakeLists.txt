install(FILES
	consistentmassassembler.hh
	hooketensor.hh
	hrzlumpedmassassembler.hh
	lobattolumpedmassassembler.hh
	operatorassembler.hh
  operatorassembler_parallel.hh
	stiffnessassembler.hh
	symmetrictensor.hh
	stiffnessassembler_truss.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/elastodynamics/assemblers)
