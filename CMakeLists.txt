cmake_minimum_required(VERSION "3.12")
project(photochem Fortran C)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(BUILD_STATIC_PHOTOCHEM "if ON, then will build the
static library libphotochem.a" ON)

option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

option(CMAKE_POSITION_INDEPENDENT_CODE "Makes code position independent." ON)

option(BUILD_F2PY_PHOTOCHEM "if ON, then will build a python
version with f2py" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" ON)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

# adding src build dependencies
add_subdirectory(src)

set(PHOTOCHEM_SOURCES ${CMAKE_SOURCE_DIR}/src/photochem_data.f90 
                      ${CMAKE_SOURCE_DIR}/src/photochem_vars.f90 
                      ${CMAKE_SOURCE_DIR}/src/photochem_wrk.f90 
                      ${CMAKE_SOURCE_DIR}/src/photochem_lightning.f90 
                      ${CMAKE_SOURCE_DIR}/src/photochem_clima.f90
                      ${CMAKE_SOURCE_DIR}/src/photochem.f90
                      ${CMAKE_SOURCE_DIR}/src/lin_alg.f)

set(PHOTOCHEM_DEPS ${CMAKE_SOURCE_DIR}/src/determine_dimensions.f90
                   ${CMAKE_SOURCE_DIR}/src/read_species.f90
                   ${CMAKE_SOURCE_DIR}/src/read_reactions.f90
                   ${CMAKE_SOURCE_DIR}/src/read_settings.f90
                   ${CMAKE_SOURCE_DIR}/src/photgrid.f90
                   ${CMAKE_SOURCE_DIR}/src/Initphoto.f90.f90
                   ${CMAKE_SOURCE_DIR}/src/Initmie.f90
                   ${CMAKE_SOURCE_DIR}/src/Rates.f90
                   ${CMAKE_SOURCE_DIR}/src/Xsections.f90
                   ${CMAKE_SOURCE_DIR}/src/Rainout.f90
                   ${CMAKE_SOURCE_DIR}/src/Aqueous.f90
                   ${CMAKE_SOURCE_DIR}/src/Ltning.f90
                   ${CMAKE_SOURCE_DIR}/src/Aertab.f90
                   ${CMAKE_SOURCE_DIR}/src/Densty.f90
                   ${CMAKE_SOURCE_DIR}/src/Aercon.f90
                   ${CMAKE_SOURCE_DIR}/src/PhotSatrat.f90
                   ${CMAKE_SOURCE_DIR}/src/Difco.f90
                   ${CMAKE_SOURCE_DIR}/src/Sedmnt.f90
                   ${CMAKE_SOURCE_DIR}/src/Dochem.f90
                   ${CMAKE_SOURCE_DIR}/src/Chempl.f90
                   ${CMAKE_SOURCE_DIR}/src/Photo.f90
                   ${CMAKE_SOURCE_DIR}/src/Rayleigh.f90
                   ${CMAKE_SOURCE_DIR}/src/Twostr.f90
                   ${CMAKE_SOURCE_DIR}/src/setup.f90
                   ${CMAKE_SOURCE_DIR}/src/integrate.f90
                   ${CMAKE_SOURCE_DIR}/src/right_hand_side.f90
                   ${CMAKE_SOURCE_DIR}/src/jacobian.f90
                   ${CMAKE_SOURCE_DIR}/src/cvode.f90
                   ${CMAKE_SOURCE_DIR}/src/redox_conservation.f90
                  )
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/photochem.f90 
                            PROPERTIES OBJECT_DEPENDS ${PHOTOCHEM_DEPS})

set_source_files_properties(${CMAKE_SOURCE_DIR}/src/photochem.f90
PROPERTIES OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/cvode_funcs.f90)

if (BUILD_STATIC_PHOTOCHEM)  
  # libphotochem.a
  add_library(photochem 
  ${PHOTOCHEM_SOURCES})
  
  target_link_libraries(photochem 
                        yaml 
                        sundials_fcvode_mod
                        sundials_cvode
                        yaml 
                        minpack)
  # target_compile_options(photochem PRIVATE -Wunused-variable)
  # target_include_directories(photochem PUBLIC
  #                           "${CMAKE_BINARY_DIR}/modules"
  #                           )

  # the main executable
  add_executable(photo.run ${CMAKE_SOURCE_DIR}/src/main.f90)
  target_link_libraries(photo.run photochem )
  # target_include_directories(photo.run PRIVATE ${CMAKE_BINARY_DIR}/modules)

  #  tests
  # add_executable(test.run ${CMAKE_SOURCE_DIR}/tests/test.f90)
  # target_link_libraries(test.run photochem)
  # target_include_directories(test.run PRIVATE ${CMAKE_BINARY_DIR}/modules)

  # enable_testing()
  # add_test(test_photochem test.run)
endif()

if (BUILD_F2PY_PHOTOCHEM)
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(F2PY REQUIRED)
  
  message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")
  
  set(f2py_module_name "_photochem")
  set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
  
  add_custom_target(${f2py_module_name} ALL
  DEPENDS ${PHOTOCHEM_SOURCES}
  )
  
  add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${PHOTOCHEM_SOURCES}
    --lower
    only: setup right_hand_side jacobian 
          integrate cvode cvode_save cvode_equilibrium
          steam2photochem :
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${PHOTOCHEM_SOURCES}
  )
  
  add_library(${generated_module_file} MODULE
              "${f2py_module_name}module.c"
              "${F2PY_INCLUDE_DIR}/fortranobject.c"
              "${f2py_module_name}-f2pywrappers2.f90"
               ${PHOTOCHEM_SOURCES})
    
  target_link_libraries(${generated_module_file} 
                        yaml 
                        sundials_fcvode_mod
                        sundials_cvode
                        yaml 
                        minpack)
  target_include_directories(${generated_module_file} PUBLIC
                            "${CMAKE_BINARY_DIR}/modules"
                            ${F2PY_INCLUDE_DIRS}
                            ${PYTHON_INCLUDE_DIRS})
                            
  set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
  set_target_properties(${generated_module_file} PROPERTIES PREFIX "")    
  
  if (UNIX)
    if (APPLE)
      set_target_properties(${generated_module_file} PROPERTIES
      LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
    else()
      set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,--allow-shlib-undefined')
    endif()
  endif()
  
  if (SKBUILD)
    install(TARGETS ${generated_module_file} DESTINATION PhotochemPy)
  else()
    install(TARGETS ${generated_module_file} DESTINATION ${CMAKE_SOURCE_DIR}/PhotochemPy)
  endif()
                    
endif()



