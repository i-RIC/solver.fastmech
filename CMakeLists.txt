cmake_minimum_required(VERSION 3.12)

project(Fastmech Fortran)

include(GNUInstallDirs)

# set CMAKE_PREFIX_PATH to the directory containing iriclibConfig.cmake
# ie -DCMAKE_PREFIX_PATH=d:/Downloads/iriclib-artifacts/iricsdk/iriclib-4.0.0/lib/cmake/iriclib
find_package(iriclib 4 REQUIRED)

# Change this to ON if solver requires OpenMP
option(IRIC_SOLVER_USE_OPENMP "Solver requires OpenMP" OFF)
if(IRIC_SOLVER_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

add_executable(Fastmech
  src/betaCDFandINV.f90
  src/betavol.f90
  src/elevdt.f90
  src/functionrough.f90
  src/initdt.f90
  src/initHdt.f90
  src/initimeta.f90
  src/rough.f90
  src/roughnessinv.f90
  src/volfrac.f90
  src/CalcCondMod.F90
  src/GridCondMod.F90
  src/GridCoordMod.F90
  src/ReadCGNS.F90
  src/WriteCGNS.F90
  src/ParmSet.f90
  src/StbBkwOutput.f90
  src/StpBkwInput.f90
  src/Support.f90
  src/xAnlzFlow.f90
  src/xFArea2.f90
  src/xFChPreis.f90
  src/xFInput.f90
  src/xFlowSupport.F90
  src/xFMatInv.f90
  src/xFNetInitConds.f90
  src/xFNetPreisMain.f90
  src/xFTSBoundCond.f90
  src/CsedMod.f90
  src/CsedMod_DT.f90
  src/Einit.f90
  src/GeometryMod.F90
  src/Median.f90
  src/RivCalcInitCond.F90
  src/RivConnectivityMod.F90
  src/RivRoughnessMod.F90
  src/RivStagrMod_jmn.F90
  src/RivVarMod.f90
  src/RivVarTimeMod.f90
  src/RivVarVertMod.f90
  src/RivVarWMod.f90
  src/RivVertMod2.f90
  src/Stagr_Console.f90
  src/stressDivMod.f90
  src/tridag.f90
  src/Uinit.f90
  src/Welcome.f90
)

# default real kind 8
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(WIN32)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /real-size:64")
  else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

target_link_libraries(Fastmech PRIVATE iriclib::iriclib)

if(IRIC_SOLVER_USE_OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(Fastmech PRIVATE OpenMP::OpenMP_Fortran)
endif()

# set stack on windows
if(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /HEAP:100000000,100000000 /STACK:1000000000,1000000000")
#else()
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fmax-stack-var-size=100000000")
endif()

install(TARGETS Fastmech DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(FILES definition.xml DESTINATION "${CMAKE_INSTALL_BINDIR}")

include(CTest)

add_subdirectory(tests)

# POST_BUILD
if (WIN32)
  # ie "-DCMAKE_PREFIX_PATH=d:/Downloads/iriclib-artifacts/iricsdk-4.0.9/iriclib-4.0.9/lib/cmake/iriclib;d:/Downloads/iriclib-artifacts/iricsdk-4.0.9/hdf5-1.8.21/cmake/hdf5;d:/Downloads/iriclib-artifacts/iricsdk-4.0.9/poco-1.9.4/lib/cmake/Poco"
  find_package(HDF5 1.8 REQUIRED)
  find_package(PocoFoundation REQUIRED)
  # copy reqd dlls to test1 directory
  add_custom_command(TARGET Fastmech POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:iriclib::iriclib>"   "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:Poco::Foundation>"   "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:hdf5::hdf5-shared>"  "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:zlib-shared>"        "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:szip-shared>"        "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\2021.2.0\\windows\\redist\\intel64_win\\compiler\\$<IF:$<CONFIG:Release>,libifcoremd.dll,libifcoremdd.dll>" "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\2021.2.0\\windows\\redist\\intel64_win\\compiler\\$<IF:$<CONFIG:Release>,libmmd.dll,libmmdd.dll>"           "$<TARGET_FILE_DIR:Fastmech>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\2021.2.0\\windows\\redist\\intel64_win\\compiler\\svml_dispmd.dll"                                          "$<TARGET_FILE_DIR:Fastmech>"
    )
endif()

# message(STATUS "HDF5_Fortran_INCLUDE_DIR = ${HDF5_Fortran_INCLUDE_DIR}")
# message(STATUS "HDF5_Fortran_LIBRARIES = ${HDF5_Fortran_LIBRARIES}")
# message(STATUS "HDF5_Fortran_HL_LIBRARIES = ${HDF5_Fortran_HL_LIBRARIES}")
# message(STATUS "CMAKE_Fortran_MODULE_DIRECTORY = ${CMAKE_Fortran_MODULE_DIRECTORY}")
# message(STATUS "HDF5_INCLUDE_DIRS = ${HDF5_INCLUDE_DIRS}")
# message(STATUS "hdf5::h5diff = $<TARGET_FILE_DIR:hdf5::h5diff>")
