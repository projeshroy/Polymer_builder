# Install script for directory: /usr/include/Eigen_3.3.1/Eigen/src/SparseCore

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen/src/SparseCore" TYPE FILE FILES
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/AmbiVector.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/CompressedStorage.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/MappedSparseMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseAssign.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseBlock.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseColEtree.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseCompressedBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseCwiseBinaryOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseCwiseUnaryOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseDenseProduct.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseDiagonalProduct.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseDot.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseFuzzy.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseMap.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseMatrixBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparsePermutation.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseProduct.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseRedux.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseRef.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseSelfAdjointView.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseSolverBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseSparseProductWithPruning.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseTranspose.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseTriangularView.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseUtil.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseVector.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/SparseView.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/SparseCore/TriangularSolver.h"
    )
endif()

