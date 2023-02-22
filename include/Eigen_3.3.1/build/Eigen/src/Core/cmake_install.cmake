# Install script for directory: /usr/include/Eigen_3.3.1/Eigen/src/Core

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen/src/Core" TYPE FILE FILES
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Array.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/ArrayBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/ArrayWrapper.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Assign.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/AssignEvaluator.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Assign_MKL.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/BandMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Block.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/BooleanRedux.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CommaInitializer.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CoreEvaluators.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CoreIterators.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CwiseBinaryOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CwiseNullaryOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CwiseUnaryOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/CwiseUnaryView.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/DenseBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/DenseCoeffsBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/DenseStorage.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Diagonal.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/DiagonalMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/DiagonalProduct.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Dot.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/EigenBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/ForceAlignedAccess.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Fuzzy.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/GeneralProduct.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/GenericPacketMath.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/GlobalFunctions.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/IO.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Inverse.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Map.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/MapBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/MathFunctions.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Matrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/MatrixBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/NestByValue.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/NoAlias.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/NumTraits.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/PermutationMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/PlainObjectBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Product.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/ProductEvaluators.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Random.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Redux.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Ref.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Replicate.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/ReturnByValue.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Reverse.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Select.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/SelfAdjointView.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/SelfCwiseBinaryOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Solve.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/SolveTriangular.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/SolverBase.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/SpecialFunctions.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/StableNorm.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Stride.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Swap.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Transpose.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Transpositions.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/TriangularMatrix.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/VectorBlock.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/VectorwiseOp.h"
    "/usr/include/Eigen_3.3.1/Eigen/src/Core/Visitor.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/usr/include/Eigen_3.3.1/build/Eigen/src/Core/products/cmake_install.cmake")
  include("/usr/include/Eigen_3.3.1/build/Eigen/src/Core/util/cmake_install.cmake")
  include("/usr/include/Eigen_3.3.1/build/Eigen/src/Core/arch/cmake_install.cmake")
  include("/usr/include/Eigen_3.3.1/build/Eigen/src/Core/functors/cmake_install.cmake")

endif()

