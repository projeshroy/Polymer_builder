# Install script for directory: /usr/include/Eigen_3.3.1/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/usr/include/Eigen_3.3.1/Eigen/Cholesky"
    "/usr/include/Eigen_3.3.1/Eigen/CholmodSupport"
    "/usr/include/Eigen_3.3.1/Eigen/Core"
    "/usr/include/Eigen_3.3.1/Eigen/Dense"
    "/usr/include/Eigen_3.3.1/Eigen/Eigen"
    "/usr/include/Eigen_3.3.1/Eigen/Eigenvalues"
    "/usr/include/Eigen_3.3.1/Eigen/Geometry"
    "/usr/include/Eigen_3.3.1/Eigen/Householder"
    "/usr/include/Eigen_3.3.1/Eigen/IterativeLinearSolvers"
    "/usr/include/Eigen_3.3.1/Eigen/Jacobi"
    "/usr/include/Eigen_3.3.1/Eigen/LU"
    "/usr/include/Eigen_3.3.1/Eigen/MetisSupport"
    "/usr/include/Eigen_3.3.1/Eigen/OrderingMethods"
    "/usr/include/Eigen_3.3.1/Eigen/PaStiXSupport"
    "/usr/include/Eigen_3.3.1/Eigen/PardisoSupport"
    "/usr/include/Eigen_3.3.1/Eigen/QR"
    "/usr/include/Eigen_3.3.1/Eigen/QtAlignedMalloc"
    "/usr/include/Eigen_3.3.1/Eigen/SPQRSupport"
    "/usr/include/Eigen_3.3.1/Eigen/SVD"
    "/usr/include/Eigen_3.3.1/Eigen/Sparse"
    "/usr/include/Eigen_3.3.1/Eigen/SparseCholesky"
    "/usr/include/Eigen_3.3.1/Eigen/SparseCore"
    "/usr/include/Eigen_3.3.1/Eigen/SparseLU"
    "/usr/include/Eigen_3.3.1/Eigen/SparseQR"
    "/usr/include/Eigen_3.3.1/Eigen/StdDeque"
    "/usr/include/Eigen_3.3.1/Eigen/StdList"
    "/usr/include/Eigen_3.3.1/Eigen/StdVector"
    "/usr/include/Eigen_3.3.1/Eigen/SuperLUSupport"
    "/usr/include/Eigen_3.3.1/Eigen/UmfPackSupport"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/usr/include/Eigen_3.3.1/build/Eigen/src/cmake_install.cmake")

endif()

