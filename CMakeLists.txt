cmake_minimum_required(VERSION 3.10)
project(ImportanceSampling)
find_package(OpenCV REQUIRED)
#find_package(Eigen3)

include_directories(
        ${OpenCV_INCLUDE_DIRS}
)

set(CMAKE_CXX_STANDARD 14)

add_executable(ImportanceSampling example3.cpp)
target_link_libraries(ImportanceSampling ${OpenCV_LIBS})
#target_link_libraries(ImportanceSampling Eigen3::Eigen)