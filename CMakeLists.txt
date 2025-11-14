#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.28)
project(DTSim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Copy all macro scripts and geometry files to the build directory

# Copy only .mac files from macros directory
file(GLOB MACRO_FILES "${PROJECT_SOURCE_DIR}/macros/*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${CMAKE_BINARY_DIR}/)

# Copy geometry directory preserving structure, excluding unwanted files
file(
  COPY ${PROJECT_SOURCE_DIR}/geometry/ 
  DESTINATION ${CMAKE_BINARY_DIR}/geometry
  PATTERN "*.py" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
  PATTERN "*template*" EXCLUDE
  PATTERN "*test*" EXCLUDE
  PATTERN "__pycache__" EXCLUDE
)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleDTSim exampleDTSim.cc ${sources} ${headers})
target_link_libraries(exampleDTSim ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(DTSim DEPENDS exampleDTSim)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
# install(TARGETS exampleDTSim DESTINATION bin)

