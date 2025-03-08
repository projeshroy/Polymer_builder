# Install script for directory: /usr/include/Eigen_3.3.1/Eigen/src/SparseLU

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen/src/SparseLU" TYPE FILE FILES
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLUImpl.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_Memory.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_Structs.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_Utils.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_column_bmod.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_column_dfs.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_gemm_kernel.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_kernel_bmod.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_panel_bmod.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_panel_dfs.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_pivotL.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_pruneL.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseLU/SparseLU_relax_snode.h"
    )
endif()

