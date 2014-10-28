project(MatModelOfEnv)
cmake_minimum_required(VERSION 2.8)

# Set Global Constant #
set(MatModelOfEnv_BUILD_DIR ${PROJECT_SOURCE_DIR}/zbuild)

# Find Source #
file(GLOB_RECURSE HEADERS  ${PROJECT_SOURCE_DIR}/lib/*.h ${PROJECT_SOURCE_DIR}/lib/*.hpp ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE SRC_LIST ${PROJECT_SOURCE_DIR}/source/*.cpp ${PROJECT_SOURCE_DIR}/main.cpp)

# Copy Resource file to build directory #
file(COPY ${PROJECT_SOURCE_DIR}/resource/Task.xls DESTINATION ${MatModelOfEnv_BUILD_DIR}/resource)
file(COPY ${PROJECT_SOURCE_DIR}/resource/Task.txt DESTINATION ${MatModelOfEnv_BUILD_DIR}/resource)

# Set Path to shared library #
link_directories("${PROJECT_SOURCE_DIR}/lib/libxl/")
link_directories("${PROJECT_SOURCE_DIR}/lib/fftw/")

# Create binary file #
add_executable(${PROJECT_NAME}
    ${SRC_LIST}
    ${HEADERS}
)

# Linking binary file wit shared library #
target_link_libraries(${PROJECT_NAME} xl fftw3)
