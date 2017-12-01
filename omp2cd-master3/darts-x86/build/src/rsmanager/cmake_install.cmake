# Install script for directory: /home/tgeng/omp2cd/omp2cd-master1/darts-x86/src/rsmanager

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/tgeng")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
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
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "bin" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/darts/lib" TYPE STATIC_LIBRARY FILES "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/build/src/rsmanager/librsmanager.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/darts/include" TYPE FILE FILES
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/Allocator.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/Defines.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/HardwareResource.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/HardwareSpecifier.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/HeapAllocator.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/HeapResource.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/InformationService.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/LimitSpecifier.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/MasterResource.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/MemoryAllocator.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/MemoryResource.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/NullAllocator.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/NullResource.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/Resource.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/ResourceManager.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/ResourceNeed.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/ShareSpecifier.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/Specifier.h"
    "/home/tgeng/omp2cd/omp2cd-master1/darts-x86/include/resourcemanager/Typing.h"
    )
endif()

