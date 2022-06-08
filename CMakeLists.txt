# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

  cmake_minimum_required(VERSION 3.0.2)
  project(LANGUAGES CXX)

# Adding IBEX

  # In case you installed IBEX in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/ibex-lib/build_install")

  find_package(IBEX REQUIRED)
  ibex_init_common() # IBEX should have installed this function
  message(STATUS "Found IBEX version ${IBEX_VERSION}")
  message(STATUS "Found IBEX inclusion ${IBEX_INCLUDE_DIRS}")

# Adding Eigen3

  # In case you installed Eigen3 in a local directory, you need
  # to specify its path with the CMAKE_PREFIX_PATH option, e.g.
  # set(CMAKE_PREFIX_PATH "~/eigen/build_install")

  find_package(Eigen3 REQUIRED NO_MODULE)
  message(STATUS "Found Eigen3 version ${EIGEN3_VERSION}")

# Adding Codac

  # In case you installed Codac in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/codac/build_install")

  find_package(CODAC REQUIRED codac)
  message(STATUS "Found Codac version ${CODAC_VERSION}")
  message(STATUS "Found Codac inclusion ${CODAC_INCLUDE_DIRS}")
  message(STATUS "Found Codac flags ${CODAC_CFLAGS}")

# Adding Codac

  # In case you installed Codac in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/codac/build_install")
  include(FindPkgConfig)
  pkg_search_module(PKG_CAPD REQUIRED capd capd-gui mpcapd mpcapd-gui)
  message(STATUS "[capd2codac] PKG_CAPD_INCLUDE_DIRS = ${PKG_CAPD_INCLUDE_DIRS}")
  message(STATUS "[capd2codac] PKG_CAPD_CFLAGS = ${PKG_CAPD_CFLAGS}")
  message(STATUS "[capd2codac] PKG_CAPD_CFLAGS_OTHER = ${PKG_CAPD_CFLAGS_OTHER}")
  message(STATUS "[capd2codac] PKG_CAPD_LDFLAGS = ${PKG_CAPD_LDFLAGS}")
  message(STATUS "[capd2codac] PKG_CAPD_LDFLAGS_OTHER = ${PKG_CAPD_LDFLAGS_OTHER}")

# Executables
set(LIBS diff_inclusion)
set(EXECS test_scenar test_equadiff test_expmat)

# Compilation

foreach(LIB1 ${LIBS})
  add_library(${LIB1} STATIC expIMat.cpp IVdouble.cpp IVparals.cpp diffincl.cpp)
  target_include_directories(${LIB1} PUBLIC ${CODAC_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
  target_compile_options(${LIB1} PUBLIC ${CODAC_CFLAGS} ${PKG_CPAD_CFLAGS})
  target_link_libraries(${LIB1} ${CODAC_LIBRARIES} Ibex::ibex ${PKG_CAPD_LDFLAGS} /home/damien/Recherche/invariant-lib/ibex-lib/filibsrc/libprim.a)
endforeach(LIB1 ${LIBS})
foreach(EXECUTABLE ${EXECS})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cpp)
  target_include_directories(${EXECUTABLE} PUBLIC ${CODAC_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
  target_link_libraries(${EXECUTABLE} ${CODAC_LIBRARIES} Ibex::ibex ${PKG_CAPD_LDFLAGS} ${LIBS} /home/damien/Recherche/invariant-lib/ibex-lib/filibsrc/libprim.a)
endforeach(EXECUTABLE ${EXECUTABLES})
