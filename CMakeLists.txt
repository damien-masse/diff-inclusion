# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

  cmake_minimum_required(VERSION 3.0.2)
  project(test_equadiff LANGUAGES CXX)

# Adding IBEX

  # In case you installed IBEX in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/ibex-lib/build_install")

  find_package(IBEX REQUIRED)
  ibex_init_common() # IBEX should have installed this function
  message(STATUS "Found IBEX version ${IBEX_VERSION}")

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

  find_package(CODAC REQUIRED)
  message(STATUS "Found Codac version ${CODAC_VERSION}")

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

# Compilation

  add_executable(${PROJECT_NAME} expIMat.cpp IVdouble.cpp IVparals.cpp diffincl.cpp test_equadiff.cpp)
  target_compile_options(${PROJECT_NAME} PUBLIC ${CODAC_CXX_FLAGS} ${PKG_CPAD_CFLAGS})
  target_include_directories(${PROJECT_NAME} PUBLIC ${CODAC_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${CODAC_LIBRARIES} Ibex::ibex ${PKG_CAPD_LDFLAGS} /home/damien/Recherche/invariant-lib/ibex-lib/filibsrc/libprim.a)
