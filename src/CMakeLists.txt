
set(SOURCES sympife.cpp parameters.cpp field_functions.cpp particle_functions.cpp
            poisson_solver.cpp vm_solver.cpp integrator.cpp
            pfem_extras.cpp mesh_extras.cpp)


# add the executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADER_FILES})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

#MESSAGE(STATUS "${MFEM_HOME}/${MFEM_BUILD}")
find_package(MFEM REQUIRED PATHS "${MFEM_HOME}/${MFEM_BUILD}")
if(MFEM_FOUND)
    MESSAGE(STATUS "MFEM is found.")
    MESSAGE(STATUS "MFEM include directory ${MFEM_INCLUDE_DIRS}")
    MESSAGE(STATUS "MFEM library directory ${MFEM_LIBRARIES}")

    # Add MFEM include directories to the current ones
    include_directories(SYSTEM ${MFEM_INCLUDE_DIRS})

    # MyLib links with MFEM
    #target_link_libraries(MyLib ${EXTRA_LIBS})

    # Add MFEM bin directories to the current ones
    link_directories(${MFEM_LIBRARY_DIR})


      if (MFEM_USE_OPENMP)
          message(STATUS "MFEM uses OpenMP.")
          find_package(OpenMP REQUIRED)
          target_link_libraries(${PROJECT_NAME} PRIVATE mfem OpenMP::OpenMP_CXX)
      else(MFEM_USE_OPENMP)
          message(STATUS "MFEM does not use OpenMP.")
          target_link_libraries(${PROJECT_NAME} PRIVATE mfem)
      endif(MFEM_USE_OPENMP)

else(MFEM_FOUND)
    message(FATAL_ERROR "You have to set a path for MFEM.")
endif(MFEM_FOUND)

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
